function [errR, errS, Srot] = evalResults(Rihat, Shat, Ri, S)
F = numel(Rihat);
useInds = randsample(F, 2);
subR_hat = Rihat(useInds);

Rot = cat(1, subR_hat{:}) \ cat(1, Ri{useInds});
[U, ~, V] = svd(Rot);
Rot = U*V';
if det(Rot) < 0
    Rot = -Rot;
end
errplus = zeros(F, 1);
errminus = zeros(F, 1);
for j = 1:F
    errplus(j) = norm(Rihat{j}*Rot + Ri{j}, 'fro');
    errminus(j) = norm(Rihat{j}*Rot - Ri{j}, 'fro');
end

if max(min([errplus, errminus], [], 2)) > 1e-1
    subR_hat{1} = -subR_hat{1};
    Rot = cat(1, subR_hat{:}) \ cat(1, Ri{useInds});
    [U, ~, V] = svd(Rot);
    Rot = U*V';
    if det(Rot) < 0
        Rot = -Rot;
    end
    errplus = zeros(F, 1);
    errminus = zeros(F, 1);
    for j = 1:F
        errplus(j) = norm(Rihat{j}*Rot + Ri{j}, 'fro');
        errminus(j) = norm(Rihat{j}*Rot - Ri{j}, 'fro');
    end
end
fprintf('\nEvaluating results without any refine ...\n')
errR = min([errplus, errminus], [], 2);
fprintf('\tThe error in camera matrix is %.2e\n', mean(errR));

P = size(S, 2);
for i = 1:F
    errplus(i) = sum(sqrt(sum((S(3*i-2:3*i, :)+Rot'*Shat(3*i-2:3*i, :)).^2)))/P;
    errminus(i) = sum(sqrt(sum((S(3*i-2:3*i, :)-Rot'*Shat(3*i-2:3*i, :)).^2)))/P;
end
errS = min([errplus, errminus], [], 2);
fprintf('\tThe error in shape matrix is %.2e\n', mean(errS));

Srot = zeros(size(Shat));
for i = 1:F
    if errplus(i) < errminus(i)
        Srot(3*i-2:3*i, :) = -Rot'*Shat(3*i-2:3*i, :);
    else
        Srot(3*i-2:3*i, :) = Rot'*Shat(3*i-2:3*i, :);
    end
end
