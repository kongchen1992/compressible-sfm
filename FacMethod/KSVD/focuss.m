%FOCUSS Focal Unversdetermined System Solver.
%   []
%   [...] = focuss(..., Name, Value,...) is used to specify some parameters:
%
function [xk, qk] = focuss(A, b, option)
if nargin < 3
    option = focussOption();
end

l = option.l;
epsilon = option.epsilon;
MaxIter = option.MaxIter;
InitX = option.InitX;
Tikhonov = option.Tikhonov;
lambda = option.lambda;
TSVD = option.TSVD;
t = option.t;
Trunc = option.Trunc;
BlurC = option.BlurC;
HardThres = option.HardThres;
p = option.p;
if isempty(option.Wak)
    Areshape = reshape(A, size(A, 1)*p, size(A, 2)/p);
    wak = max(abs(Areshape), [], 1);
    wak = ones(size(wak)) ./ wak;
    Wak = diag(kron(wak, ones(1, p)));
else
    Wak = option.Wak;
    wak = diag(Wak)';
    wak = wak(1:p:end);
end

xk = cell(MaxIter+1, 1);
qk = cell(MaxIter, 1);

% temporary, need to be removed later
Fqcell = zeros(MaxIter, size(A, 2)/p);

% Initialization
if isempty(InitX)
    % Compute A\b only when InitX is not provided.
    % Computation saving.
    InitX = A\b;
end
% Blur the sparse initial point.
FX = norm(InitX, 'fro');
if FX < 1e-5
    xk{1} = InitX + BlurC*randn(size(InitX));
else
    xk{1} = InitX + BlurC*FX*randn(size(InitX));
end

for k = 1:MaxIter
    % Compute Wpk
    if p == 1 && size(xk{k}, 2) == 1
        Wpk = diag(xk{k}.^l);
    else
        xkReshape = reshape((xk{k}.^l)', size(xk{k}, 2)*p, size(xk{k}, 1)/p);
        wpk = sqrt(sum(xkReshape.^2, 1)); % maybe exist better design
        Wpk = diag(kron(wpk, ones(1, p)));
    end

    % Compute qk
    Aw = A*Wak*Wpk;
    if Tikhonov
        % Tikhonov Regularization
        if isempty(lambda)

        end
        qk{k} = Aw'*((Aw*Aw' + lambda*eye(size(Aw, 1)))\b);
    elseif TSVD
        % Truncated SVD
        [U, S, V] = svd(Aw, 'econ');
        invS = diag(1/diag(S(1:t, 1:t)));
        qk{k} = V(:, 1:t)*invS*U(:, 1:t)'*b;
    elseif Trunc
        qk{k} = zeros(size(Aw, 2), size(b, 2));
        qk{k}(diag(Wpk)~=0, :) = (Aw(:, diag(Wpk)~=0))\b;
    else
        % no regularization
        qk{k} = (Aw)\b;
    end

    % Compute Frobenius norm of each block in q
    if p == 1 && size(xk{k}, 2) == 1
        Fq = qk{k};
    else
        qkReshape = reshape((qk{k})', size(qk{k}, 2)*p, size(qk{k}, 1)/p);
        Fq = sqrt(sum(qkReshape.^2, 1)).*wak;
        Fqcell(k, :) = Fq;
    end

    % hard thresholding
    if HardThres && p == 1
        qk{k}(abs(Fq) < epsilon) = 0;
    elseif HardThres
        qk{k}(kron(Fq, ones(p, 1)) < epsilon, :) = 0;
    end

    % Compute x_k
    xk{k+1} = Wak*Wpk*qk{k};
    
    % Check if converged
    indZeros = (abs(Fq) < 1e-6);
    indOnes = (abs(Fq - 1) < 1e-6);
    if sum(indZeros + indOnes) == numel(Fq)
        break;
    end
end

xk = xk(2:k+1);
qk = qk(1:k);
