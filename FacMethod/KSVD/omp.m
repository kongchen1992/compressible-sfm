%OMP Orthogonal Matching Pursuit
%   [S_HAT, ERR] = OPM(PHI, V, M) finds a M-P-sparse estimateion S_HAT for
%   signal V via OMP method.  PHI is the known dictionary. ERR is the error
%   after each iteration.  Notation: M-P-sparse means there are m active
%   blocks whose 1st dimension is P.
%
%   Notation: Make sure Phi is normalized.

function [s_hat, err] = omp(Phi, v, m, sizeG, U)
p = sizeG(1); q = sizeG(2);
K = size(Phi,2)/p;
R = v; 
Lambda = zeros(size(Phi, 2)/p, 1);
LAMBDA = zeros(size(Phi, 2), 1);
err = zeros(m, 1);
for i = 1:m
    
    % find the biggest lambda
    if p == 1 && size(v, 2) == 1
        dotProds = abs((Phi'*R));
    else
        URf = reshape(R'*U, p*q, K);
        dotProds = sum(URf.^2, 1);
    end
    dotProds(Lambda==1) = 0;
    [~,I] = max(dotProds);
    Lambda(I) = 1;
    LAMBDA(p*(I-1)+1:p*I) = 1;
    Phit = Phi(:, LAMBDA~=0);
    
    xt = Phit\v;
    
    R = v - Phit*xt;
    
    err(i) = norm(R, 'fro');
end
s_hat = zeros(size(Phi, 2), size(v, 2));
s_hat(LAMBDA~=0, :) = xt;
