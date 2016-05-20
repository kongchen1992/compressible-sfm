function exp_noise_I(ratios, I)
Dir = '~/GitHub/Chen-CVPR-2016/data/Noise';
if ~exist(Dir, 'dir')
    mkdir(Dir);
end
file = fullfile(Dir, sprintf('noise_%02d.mat', I));
K = 2;
L = 5;
F = 100;
P = 30;

if exist(file, 'file')
    load(file);
    W = GroundTruth.W;
    B = GroundTruth.B;
    Pi = GroundTruth.Pi;
    S = GroundTruth.S;
    Ri = GroundTruth.Ri;
else
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
    noise = randn(size(W));
end

numR = numel(ratios);
for i = 1:numR
    ratio = ratios(i);
    fprintf('RATIO of NOISE is %.4f\n', ratio);
    W_noise = W + (10^ratio)*norm(W, 'fro')/norm(noise, 'fro')*noise;

    [~, nproc] = unix('nproc');
    NumCore = str2num(nproc)/2;

    admmOpt = admmOption('MaxIter', 3000);
    ksvdOpt = ksvdOption('IsOMP', true, ...
        'lambdaTau', 0.9, 'MaxIter', 100, 's', K, 'NumCore', NumCore, ...
        'IsPrune', false, 'IsCut', false, 'ThresReplace', 0.8*F*K/L);
    facOpt = facOption('ksvdOption', ksvdOpt, 'admmOption', admmOpt);
    [Pihat, Bhat, ~] = facMethod(W_noise, L, facOpt);

    T = computeT(B', Bhat');
    figure;
    imagesc(T);
    errB = norm(Bhat - T'*B, 'fro');
    errPi = norm(Pihat*T' - Pi, 'fro');
    fprintf('The representation error of B is %.4e\n', errB);
    fprintf('The representation error of Pi is %.4e\n', errPi);

    [~, Shat, Rihat, ~] = solveProblem(Pihat, Bhat, [], 0);
    [errR, errS, Srot] = evalResults(Rihat, Shat, Ri, S);
    Results(i).Pihat = Pihat;
    Results(i).Bhat = Bhat;
    Results(i).Shat = Shat;
    Results(i).Rihat = Rihat;
    Results(i).errR = errR;
    Results(i).errS = errS;
    Results(i).Srot = Srot;
    close all;
end

GroundTruth.Pi = Pi;
GroundTruth.B = B;
GroundTruth.W = W;
GroundTruth.S = S;
GroundTruth.Ri = Ri;
save(file, 'Results', 'GroundTruth', 'noise', '-v7.3');
