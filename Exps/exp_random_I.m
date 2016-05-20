%This script is the 3rd Experiments of synthetic data.
function exp_random_I(K, Ls, F, P, I)
%F = 500; P = 100; L = 20; K = 2;
for l = 1:numel(Ls)
    L = Ls(l);
    Bsharp = zeros(L, 3*P);
    B = randn(3*L, P);
    B = B/diag(sqrt(sum(B.^2, 1)));
    Bsharp(:, 1:3:end) = B(1:3:end, :);
    Bsharp(:, 2:3:end) = B(2:3:end, :);
    Bsharp(:, 3:3:end) = B(3:3:end, :); 

    Csharp = zeros(F, L);
    for i = 1:F
        Csharp(i, randsample(L, K)) = randn(K, 1);
    end
    Ri = camera_setup(F, 'continuous');
    R = blkdiag(Ri{:});

    Ssharp = Csharp*Bsharp;
    C = kron(Csharp, eye(3));
    Pi = R * C;
    S = C * B;
    W = R * S;

    [~, nproc] = unix('nproc');
    NumCore = str2num(nproc)/2;

    admmOpt = admmOption('MaxIter', 3000);
    ksvdOpt = ksvdOption('IsOMP', true, ...
        'lambdaTau', 0.9, 'MaxIter', 100, 's', K, 'NumCore', NumCore, ...
        'IsPrune', false, 'IsCut', false, 'ThresReplace', 0.8*F*K/L);
    facOpt = facOption('ksvdOption', ksvdOpt, 'admmOption', admmOpt);
    [Pihat, Bhat, results] = facMethod(W, L, facOpt);
    if results.cost(end) > 1e-3
        exp_random_I(K, L, F, P, I);
        continue;
    end

    T = computeT(B', Bhat');
    figure;
    imagesc(T);
    errB = norm(Bhat - T'*B, 'fro');
    errPi = norm(Pihat*T' - Pi, 'fro');
    fprintf('The representation error of B is %.4e\n', errB);
    fprintf('The representation error of Pi is %.4e\n', errPi);

    [Rhat, Shat, Rihat, debugDetail] = solveProblem(Pihat, Bhat, [], 0);
    [errR, errS, Srot] = evalResults(Rihat, Shat, Ri, S);
    Dir = '~/GitHub/Chen-CVPR-2016/data/RANDOM';
    file = fullfile(Dir, sprintf('L%03dF%04dP%03d_%02d.mat', L, F, P, I));
    save(file);
    close all;
end
