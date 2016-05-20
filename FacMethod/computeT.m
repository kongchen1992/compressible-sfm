%COMPUTET computes the transformation between Bhat and B.
function T = computeT(B, Bhat)

K = size(B, 2)/3;
Uj = cell(K, 1);
for j = 1:K
    [Uj{j}, ~, ~] = svd(B(:, 3*(j-1)+1:3*j), 'econ');
end

Bj = cell(K, 1);
for j = 1:K
    Bj{j} = B(:, 3*j-2:3*j);
end

T = zeros(size(B, 2));
valid = 1:numel(Uj);
for j = 1:K
    U = cat(2, Uj{valid});
    t = omp(cat(2, Bj{valid}), Bhat(:, 3*(j-1)+1:3*j), 1, [3,3], U);
    I = ceil(find(t(:, 1)~=0, 1)/3);
    T(3*valid(I)-2:3*valid(I), 3*(j-1)+1:3*j) = t(3*I-2:3*I, :);
    valid = valid([1:I-1, I+1:end]);
end
