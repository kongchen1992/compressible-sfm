function exp_jaws_structures()
Dir = '~/GitHub/Chen-CVPR-2016/data/jaws';
file = fullfile(Dir, 'jaws.mat');
data = load(file);
S = data.S;
F = size(S, 1)/3;
P = size(S, 2);
K = 2;
L = 20;

Ssharp = zeros(F, 3*P);
Ssharp(:, 1:3:end) = S(1:3:end, :);
Ssharp(:, 2:3:end) = S(2:3:end, :);
Ssharp(:, 3:3:end) = S(3:3:end, :);

Sprior = diag(sqrt(sum(Ssharp.^2, 1)));
SsharpP = Ssharp/Sprior;

admmOpt = admmOption('MaxIter', 5000, 'initRho', 1e-1, 'initMu', 1e2, ...
    'maxRho', 1e1, 'maxMu', 1e6, 'tauRho', 1.001, 'tauMu', 1.01, 'gamma', 10);
[~, ~, ~, Xinit, Dinit] = gs_Factorization(Ssharp', L, [1, 1], admmOpt);

[~, nproc] = unix('nproc');
focussOpt = focussOption('p', 1, 'Trunc', true, 'Tikhonov', true, ...
    'lambda', 6e-3, 'BlurC', 1);
ksvdOpt = ksvdOption('IsFOCUSS', true, 'focussOption', focussOpt, ...
    'lambdaTau', 0.9, 'MaxIter', 30, 'sMin', 1, 'sMax', K, 'IsPrune', false, ...
    'LocalTrap', 1e-1, 'ThresReplace', 0.8*F*2/L, 'NumCore', str2num(nproc)/2, ...
    'IsCut', false);

ksvdOpt.Dinit = Dinit;
ksvdOpt.Xinit = Xinit;

[Dksvd, Xhat, details] = ksvd(SsharpP', L, [1, 1], ksvdOpt);

Dhat = Sprior'*Dksvd;
Ssharphat = (Dhat*Xhat)';

Shat = zeros(size(S));
Shat(1:3:end, :) = Ssharphat(:, 1:3:end);
Shat(2:3:end, :) = Ssharphat(:, 2:3:end);
Shat(3:3:end, :) = Ssharphat(:, 3:3:end);

[errS, Srot] = compareStructs(S, Shat);

file = fullfile(Dir, 'structures.mat');
save(file);
