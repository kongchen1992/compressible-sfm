% [B, cost] = gs_Coding(X, A, sizeG) performs group sparse 
%   coding using Alternating Direction Method of Multipliers.  X is Data 
%   matrix. A is Codebook matrix. sizeG indicates the size of each group. 
%   
%   Specificially, the group sparse coding problem is group LASSO problem
%   arg min_{B} 0.5*||X - A*B||_{F}^{2} + gamma*||B||_{21}
%
%   [B, cost] = gs_Coding(X, sizeG, facOpt) 
%   specifies optional parameters from facOpt.  Check facOption for detail.
%
%   Example:
%       X = randn(2000, 40);
%       opt = facOption('gamma', 1, 'maxIter', 1e6);
%       [B, cost] = gs_Coding(X, A, [3, 2], opt);
%
% Written by Chen Kong 2014

function [B,cost, C] = gs_Coding(X, A, sizeG, option)
% Get dimensions
p = sizeG(1); q = sizeG(2);
[m, qn] = size(X);
n = qn/q;
L = size(A, 2)/p;

gamma = option.gamma;
maxIter = option.maxIter;
isQuiet = option.isQuiet;
rho = option.initRho;
maxRho = option.maxRho;
tauRho = option.tauRho;
stopCritieria = option.stopCritieria;
stUpdate = option.stUpdate;

% Initialize all
L_B = zeros(L*p, n*q);% Initialize Lagragian to be nothing (seems to work well)
I = speye(L*p); % Set the sparse identity matrix

% B = option.initB;
B = randn(L*p, n*q);
C = randn(L*p, n*q);
g = getGroupIdx(size(B), sizeG);

% Set the norm functions
norm2 = @(x) x(:)'*x(:); 
norm21 = @(x) sum(sqrt(sum(x.^2, 1)));

cost = zeros(maxIter/200, 1);
% costL = zeros(maxIter, 1);

if ~isQuiet
    fprintf('\nUsing Alternating Direction Method of Multipliers.');
    fprintf('\n------------------------------------------------------------------------------')
    fprintf('\nDim of X:\t%d-by-%d,\tDim of block:\t%d-by-%d,', ...
        size(X, 1), size(X, 2), sizeG(1), sizeG(2));
    fprintf('\nDim of A:\t%d-by-%d,\tDim of B:\t%d-by-%d,', ...
        m, L*p, size(B, 1), size(B, 2));
    fprintf('\ngamma = \t%0.4f,\tmaxIter = \t%d.', gamma, maxIter);
    fprintf('\n****************************************************************************\n');
    fprintf('   it     mu    |X-AB|_2   |C|_g0    |B-C|    cost   d(cost) time\n');
    fprintf('------------------------------------------------------------------------------');
end

numStill = 0; tic;

for i = 1:maxIter

    % Solve sub-problem to solve B
    B = (A'*A + rho*I)\(A'*X + rho*C - L_B);
    
    % Solve sub-problem to solve C
    C(g) = fast_sthresh(B(g) + 1/rho*L_B(g), gamma/rho);
    
    % Update the Lagrangian
    L_B = L_B + rho*(B - C);
    
    % get the current cost
    % since it's computational expensive to compute cost per iteration, we
    % compute cost per 100 iteration.
%     cost(i) = 0.5*norm2(X - A*B) + gamma*norm21(B(g)); 
    
%     if i > 1
%         if abs(cost(i) - cost(i-1)) < 1e-5*cost(i)
    if i > stUpdate
        rho = min(maxRho, rho*tauRho);
    end
%         end
%     end
    
    if ~isQuiet
        if numStill == 3
            fprintf('\n------------------------------------------------------------------------------\n');
            fprintf('Stoped.\n');
            fprintf('It seems cost does not change anymore.\n')
            cost = cost(1:floor(i/200));
            break;
        end
        if mod(i, 200) == 0
            % stop criterion
            cost(i/200) = 0.5*norm2(X - A*B) + gamma*norm21(B(g));
                
            if i/200 == 1
                fprintf('\n%6d|%1.1e|%1.2e|%10d|%1.1e|%1.1e|%+1.2e|%.2f', i, rho, ...
                    norm(X - A*B, 'fro'), numel(find(sum(abs(C(g))>1e-5, 1) ~= 0)), ...
                    norm(B - C, 'fro'), cost(i/200), cost(i/200), toc);
            else
                fprintf('\n%6d|%1.1e|%1.2e|%10d|%1.1e|%1.1e|%+1.2e|%.2f', i, rho, ...
                    norm(X - A*B, 'fro'), numel(find(sum(abs(C(g))>1e-5, 1) ~= 0)), ...
                    norm(B - C, 'fro'),  cost(i/200), (cost(i/200) - cost(i/200 - 1))/200, toc);
                if abs(cost(i/200) - cost(i/200-1)) < stopCritieria*cost(i/200)*200
                    numStill = numStill + 1;
                else
                    numStill = 0;
                end
            end
        end
    end
end
if ~isQuiet
    if i == maxIter
        fprintf('\n------------------------------------------------------------------------------\n');
        fprintf('Stoped.\n');
        fprintf('It reaches max Iteration.\n');
    end
%     if numel(find(sum(abs(B(g))>1e-5) ~= 0)) > 0.6*numel(B)/prod(sizeG)
%         fprintf('The codebook is not that group sparse. Seems something wrong.\n');
%     end
end