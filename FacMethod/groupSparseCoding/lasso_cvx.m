% Function to perform LASSO regression using CVX
%
% arg min_{B} 0.5*||X - A*B||_{2}^{2} + gamma*||B||_{1}
%
% Usage:- [B,cost] = lasso_admm(X, A, gamma) 
%
% where:- <in>
%         b = bias vector
%         lambda = weighting on the l1 penalty
%         <out>
%         x = solution  
%
% Written by Simon Lucey 2012

function B = lasso_cvx(X, A, gamma)

% Get the size of B
c = size(X,2); 
r = size(A,2); 

% Now vectorize everything
I = speye(c); 
J = kron(I,A); 
x = X(:); 

cvx_begin  
    variable b(r*c,1);
    minimize(0.5*sum_square(x - J*b) + gamma*norm(b,1)); 
cvx_end
B = reshape(b,[r,c]); 

