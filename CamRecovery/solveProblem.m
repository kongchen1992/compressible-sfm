% Function of Non-Rigid Structure from Motion using Group Sparsity
%   algorithm. This is the main algorithm, whose inputs are decomposition
%   results and outputs are camera matrix R, shape matrix S.
%
%   [R, S] = solveProblem(Pi_hat, B_hat, piStruct, isviz, GT)
%
%   Inputs:  
%       -Pi_hat:    a 2F-by-3L matrix from decomposition,
%       -B_hat:     a 3L-by-P matrix from decomposition,
%       -piStruct:  a matrix describing the sparse structure of Pi_hat,
%                   where each elements indicates number of active elements
%                   in Pi_hat's corresponding group. It's computed by 
%                   function gs_Struct.
%
%   Outputs:
%       -R:         a 2F-by-3F block diagonal matrix, whose each block is 
%                   camera matrix for corresponding frame.
%       -S:         a 3F-by-P matrix.  Each triple rows indicates 3D
%                   coordinates of key points in corresponding frame.
%       -Ri:        a F-by-1 cell array. Each cell contains camera matrix
%                   of corresponding frames.
%       -R0:        a rotation matrix that rotate Ri{1} to GT.gt.Ri{1}.

function [R, S, Ri, debugDetail] = solveProblem(Pi_hat, B_hat, GT, isviz, solveDAlg)

if nargin < 4
    isviz = 0;
    solveDAlg = 'eig';
end
if nargin < 5
    solveDAlg = 'eig';
end

L = size(B_hat, 1)/3;
F = size(Pi_hat, 1)/2;
P = size(B_hat, 2);

tic;

piStruct = gs_Struct(Pi_hat, [2,3]);

%% solve D
fprintf('\nSplit the problem into %d subproblems.\n', L);
fprintf('\tSolving corrective matrix...\n');

Dj = cell(L, 1);
C_tilde = zeros(F, L);
R_hat = cell(F, L);

% numCores = 4;
% fprintf('--------------------------------------------------------------------\n')
% fprintf('Starting a pool of workers with %d cores\n', numCores);
% parpool('local',numCores);
% 
% parfor j = 1:numel(Dj)
for j = 1:numel(Dj)
    
    piCol = Pi_hat(:, 3*j-2:3*j);    
    idx = kron(piStruct(:, j), [1;1]) ~= 0;
    delta  = piCol(idx, :);    
    if isempty(delta)
        Dj{j} = eye(3);
        continue;
    end
    Aij = cell(size(delta, 1)/2, 1);
    
    for i = 1:size(delta, 1)/2
        
        delta_ij = delta(2*i-1:2*i, :);        
        kronD = kron(delta_ij, delta_ij);        
        Aij{i} = [kronD(1,:) - kronD(4, :); kronD(2, :)];
    end
    
    Dj{j} = solveDj(Aij, solveDAlg);%, Qj_gt);
end

% fprintf('Closing the pool\n');
% delete(gcp);
% toc;
% fprintf('--------------------------------------------------------------------\n')
%

fprintf('\tSolving coefficients in Sparse Coding...\n\t')
for j = 1:numel(Dj)
    piCol = Pi_hat(:, 3*j-2:3*j);    
    idx = kron(piStruct(:, j), [1;1]) ~= 0;    
    delta  = piCol(idx, :);  
    % solve C_tilde and sign ambiguity
    tau = find(piStruct(:, j) ~= 0);
    for i = 1:numel(tau)
        delta_ij = delta(2*i-1:2*i, :);
        c_tilde = mean(sqrt(diag(delta_ij*(Dj{j})*Dj{j}'*delta_ij')));
        r_hat = delta_ij*Dj{j}/c_tilde;
        R_hat{tau(i), j} = r_hat;
        
        C_tilde(tau(i), j) = c_tilde;
    end
end
toc;

debugDetail.R_hat = R_hat;

%% solve R
fprintf('\nSolving camera matrix of each frames and shape matrix...\n\t');

existR = (piStruct ~= 0);

[~, indx] = max(sum(double(existR)));

done = false(1, L);
gottenR = false(F, 1);
Ri = cell(F, 1);

%temp_begin
csign = zeros(size(C_tilde, 1), 1);
%temp_end

for i = 1:L
    
    overlap = gottenR & existR(:, indx);
    getR = (~gottenR) & existR(:, indx);
    gottenR = gottenR | getR;
    
    if i == 1
        Rot = eye(3);
        %temp_begin
        csign(getR) = sign(C_tilde(getR, indx));
        %temp_end
    else
        Rot = cat(1, R_hat{overlap, indx}) \ cat(1, Ri{overlap});
        %temp_begin
        OverInds = find(overlap);
        useInds = randsample(OverInds, 2);
        subR_hat = R_hat(useInds, indx);

        Rot = cat(1, subR_hat{:}) \ cat(1, Ri{useInds});
        [U, S, V] = svd(Rot);
        Rot = U*V';
        if det(Rot) < 0
            Rot = -Rot;
        end
        errplus = zeros(numel(OverInds), 1);
        errminus = zeros(numel(OverInds), 1);
        for j = 1:numel(OverInds)
            errplus(j) = norm(R_hat{OverInds(j), indx}*Rot + Ri{OverInds(j)}, 'fro');
            errminus(j) = norm(R_hat{OverInds(j), indx}*Rot - Ri{OverInds(j)}, 'fro');
        end

        if max(min([errplus, errminus], [], 2)) > 1e-1
            subR_hat{1} = -subR_hat{1};
            Rot = cat(1, subR_hat{:}) \ cat(1, Ri{useInds});
            [U, S, V] = svd(Rot);
            Rot = U*V';
            if det(Rot) < 0
                Rot = -Rot;
            end
            errplus = zeros(numel(OverInds), 1);
            errminus = zeros(numel(OverInds), 1);
            for j = 1:numel(OverInds)
                errplus(j) = norm(R_hat{OverInds(j), indx}*Rot + Ri{OverInds(j)}, 'fro');
                errminus(j) = norm(R_hat{OverInds(j), indx}*Rot - Ri{OverInds(j)}, 'fro');
            end
        end
        %assert(max(min([errplus, errminus], [], 2)) < 1e-1)
        for j = 1:numel(OverInds)
            if errplus(j) < errminus(j)
                R_hat{OverInds(j), indx} = -R_hat{OverInds(j), indx};
                C_tilde(OverInds(j), indx) = -C_tilde(OverInds(j), indx);
            end
        end
        %temp_end
    end
    
    R_hat(existR(:, indx), indx) = cellfun(@(x) x*Rot, R_hat(existR(:, indx), indx), 'UniformOutput', false);
    
    Ri(getR) = R_hat(getR, indx);
    %temp_begin
    csign(getR) = sign(C_tilde(getR, indx));
    %temp_end
    
    Dj{indx} = Dj{indx}*Rot;
    
    done(indx) = true;
    [Y, indx] = max(sum(double(repmat(gottenR, 1, L) & existR & repmat(~done, F, 1))));
    if Y == 0
        if isempty(find(~gottenR, 1))
            break;
        else
            error('Some camera matrices are missed...');
        end
    end
end
toc;

%% get results
D = blkdiag(Dj{:});

R = blkdiag(Ri{:});

S = kron(C_tilde, eye(3))/D*B_hat; 


%% show result
if isviz
    fprintf('\nVisualizing results...\n');
    S_gt = GT.gt.S;
    viz_S(S, S_gt);
end

