%EXTACTSP Extract Sparse Pattern
%   OMEGA = EXTRACTSP(X, sizeG) creats a matrix OMEGA that extracts active 
%   part by X*OMEGA i.e. XOMEGA are active parts of X. X is a row (group)
%   sparse vector. sizeG indicates the size of each block in X.
%
%   [OMEGA, PATTERN] = EXTRACTSP(X, sizeG) also creates a matrix PATTERN 
%   indicating the sparse pattern of X.

function [Omega, pat] = extractSP(X, sizeG)

assert(size(X, 1) == sizeG(1));

reshapX = reshape(X, sizeG(1)*sizeG(2), numel(X)/(sizeG(1)*sizeG(2)));
pat = max(abs(reshapX), [], 1) >= 1e-8;

indActive = find(pat);
Omega = zeros(size(X, 2)/sizeG(2), numel(indActive));
Omega(sub2ind(size(Omega), indActive, 1:numel(indActive))) = 1;

pat = kron(pat, ones(1, sizeG(2)));
Omega = kron(Omega, eye(sizeG(2)));