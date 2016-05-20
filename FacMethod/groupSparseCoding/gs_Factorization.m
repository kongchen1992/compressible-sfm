% [A, B, cost] = gs_Factorization(X, sizeG, option) performs group sparse 
%   coding using Alternating Direction Method of Multipliers.  X is Data 
%   matrix. sizeG indicates the size of each group. 
%   
%   Specificially, the factorization problem is
%   arg min_{A,B} 0.5*||X - A*B||_{F}^{2} + gamma*||B||_{21}
%   subject to a_{i}^T * a_{i} <= 1, i = 1,...,size(A,2).
%
%   [A, B, cost] = gs_Factorization(X, sizeG, facOpt) 
%   specifies optional parameters from facOpt.  Check facOption for detail.
%
%   Example:
%       X = randn(2000, 40);
%       opt = facOption('gamma', 1, 'maxIter', 1e6, 'L', 100);
%       [A, B, cost] = gs_Fatorization(X, [3, 2], opt);
%
% Written by Chen Kong 2014

function [A, B,cost, C, D] = gs_Factorization(X, L, sizeG, option)
% Get dimensions
p = sizeG(1); q = sizeG(2);
[m, qn] = size(X);
n = qn/q;

gamma = option.gamma;
maxIter = option.maxIter;
isQuiet = option.isQuiet;
mu = option.initMu;
rho = option.initRho;
maxRho = option.maxRho;
maxMu = option.maxMu;
tauRho = option.tauRho;
tauMu = option.tauMu;
stopCritieria = option.stopCritieria;
stUpdate = option.stUpdate;
vizIter = option.vizIter;
tauGamma = option.tauGamma;
minGamma = option.minGamma;

if isempty(option.initA)
    A = randn(m, p*L);
else
    A = option.initA;
end
if isempty(option.initB)
    B = randn(p*L, q*n);
else
    B = option.initB;
end

% Initialize all
L_A = zeros(m, L*p); % Initialize Lagragian to be nothing (seems to work well)
L_B = zeros(L*p, n*q);
I = speye(L*p); % Set the sparse identity matrix

C = randn(L*p, n*q);
D = randn(m, L*p);
g = getGroupIdx(size(B), sizeG);

% Set the norm functions
norm2 = @(x) x(:)'*x(:); 
norm21 = @(x) sum(sqrt(sum(x.^2, 1)));

cost = zeros(maxIter/vizIter, 1);
% costL = zeros(maxIter, 1);

if ~isQuiet
    fprintf('\nUsing Alternating Direction Method of Multipliers.');
    fprintf('\n--------------------------------------------------------------------------------------------')
    fprintf('\nDim of X:\t%d-by-%d,\tDim of block:\t%d-by-%d,', ...
        size(X, 1), size(X, 2), sizeG(1), sizeG(2));
    fprintf('\nDim of A:\t%d-by-%d,\tDim of B:\t%d-by-%d,', ...
        m, L*p, size(B, 1), size(B, 2));
    fprintf('\ngamma = \t%0.4f,\tmaxIter = \t%d.', gamma, maxIter);
    fprintf('\n******************************************************************************************\n');
    fprintf('   it   gamma    rho      mu   |X-AB|_2   |C|_g0    |B-C|   |A-D|    cost   d(cost) time\n');
    fprintf('---------------------------------------------------------------------------------------------');
end

numStill = 0; tic;

for i = 1:maxIter

    % Solve sub-problem to solve C
    C(g) = fast_sthresh(B(g) + 1/rho*L_B(g), gamma/rho);

    % Solve sub-problem to solve D
    Dtilde = A + L_A/mu;
%    D_sharp = reshape(Dtilde, m*p, L);
%    normD = max(sqrt(sum(D_sharp.^2, 1)), 1);
%    D = Dtilde ./ kron(normD, ones(m, p));
    normD = max(sqrt(sum(Dtilde.^2, 1)), 1);
    D = Dtilde / diag(normD);

    % Solve sub-problem to solve A
    A = (X*B' + mu*D - L_A)/(B*B' + mu* I); 

    % Solve sub-problem to solve B
    B = (A'*A + rho*I)\(A'*X + rho*C - L_B);

    % Update the Lagrangian
    L_A = L_A + mu*(A - D);
    L_B = L_B + rho*(B - C);

    % get the current cost
    % since it's computational expensive to compute cost per iteration, we
    % compute cost per 100 iteration.
%     cost(i) = 0.5*norm2(X - A*B) + gamma*norm21(B(g)); 

%     if i > 1
%         if abs(cost(i) - cost(i-1)) < 1e-5*cost(i)
    if i > stUpdate
        rho = min(maxRho, rho*tauRho); 
        mu = min(maxMu, mu*tauMu);
        gamma = max(minGamma, gamma*tauGamma);
    end

%         end
%     end

    if ~isQuiet
        if numStill == 3
            fprintf('\n------------------------------------------------------------------------------\n');
            fprintf('Stoped.\n');
            fprintf('It seems cost does not change anymore.\n')
            cost = cost(1:floor(i/vizIter));
            return;
        end
        if mod(i, vizIter) == 0
            % stop criterion
            cost(i/vizIter) = 0.5*norm2(X - A*B) + gamma*norm21(B(g));
                
            if i/vizIter == 1
                fprintf('\n%6d|%1.1e|%1.1e|%1.1e|%1.2e|%10d|%1.1e|%1.1e|%1.1e|%+1.2e|%.2f', i, gamma, rho, mu, ...
                    norm(X - A*B, 'fro'), numel(find(sum(abs(C(g))>1e-5, 1) ~= 0)), ...
                    norm(B - C, 'fro'), norm(A-D, 'fro'), cost(i/vizIter), cost(i/vizIter), toc);
            else
                fprintf('\n%6d|%1.1e|%1.1e|%1.1e|%1.2e|%10d|%1.1e|%1.1e|%1.1e|%+1.2e|%.2f', i, gamma, rho, mu, ...
                    norm(X - A*B, 'fro'), numel(find(sum(abs(C(g))~=0, 1) ~= 0)), ...
                    norm(B - C, 'fro'), norm(A-D, 'fro'), cost(i/vizIter), (cost(i/vizIter) - cost(i/vizIter - 1))/vizIter, toc);
                if abs(cost(i/vizIter) - cost(i/vizIter-1)) < stopCritieria*cost(i/vizIter)*vizIter
                    numStill = numStill + 1;
                else
                    numStill = 0;
                end
            end
        end
    end
end
if ~isQuiet
    fprintf('\n------------------------------------------------------------------------------\n');
    fprintf('Reach max Iteration.  It may be caused by a bad initialization.\n');
    fprintf('Restart.\n');
end
[A, B,cost, C, D] = gs_Factorization(X, L, sizeG, option);
