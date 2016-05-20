function [Pihat, Bhat, S_clean, Ri] = debug_jaws_factorization()
[~, nproc] = unix('nproc');
NumCore = str2num(nproc)/2;
Dir = '~/GitHub/Chen-CVPR-2016/data/jaws';
file = fullfile(Dir, 'structures.mat');
data = load(file);

F = data.F;
P = data.P;
K = data.K;
L = data.L;
Bsharp = data.Dhat';
Csharp = data.Xhat';

noiseR = 0;
Bsharp = Bsharp + noiseR*std(Bsharp(:))*randn(size(Bsharp));
Bprior = diag(sqrt(sum(Bsharp.^2, 2)));
Bsharp = inv(Bprior)*Bsharp;
Csharp = Csharp*Bprior;
Ssharp_clean = Csharp*Bsharp;

fprintf('The coherence of dictionary is %.4f\n', coherence(Bsharp'));

S_clean = zeros(3*F, P);
S_clean(1:3:end, :) = Ssharp_clean(:, 1:3:end);
S_clean(2:3:end, :) = Ssharp_clean(:, 2:3:end);
S_clean(3:3:end, :) = Ssharp_clean(:, 3:3:end);

Ri = camera_setup(F, 'continuous');
R = blkdiag(Ri{:});
W_clean = R*S_clean;

B = zeros(size(Bsharp, 1)*3, size(Bsharp, 2)/3);
B(1:3:end, :) = Bsharp(:, 1:3:end);
B(2:3:end, :) = Bsharp(:, 2:3:end);
B(3:3:end, :) = Bsharp(:, 3:3:end);

Pi = R*kron(Csharp, eye(3));

UA = sum(Csharp~=0, 1);
ratioR = min(UA)/(K*F/L);
sMin = 1;
PriorConst = mean(sqrt(sum(W_clean.^2, 1))./sqrt(sum(S_clean.^2, 1)));
Sprior = diag(sqrt(sum(S_clean.^2, 1)));

admmOpt = admmOption('MaxIter', 5000, 'gamma', 0.06, ...
    'initRho', 1e-1, 'initMu', 1e1, 'maxRho', 1e1, 'maxMu', 1e5, ...
    'tauRho', 1.001, 'tauMu', 1.01);

focussOpt = focussOption('p', 3, 'Trunc', true, 'Tikhonov', true, ...
    'lambda', 6e-6, 'BlurC', 1);
ksvdOpt = ksvdOption('IsFOCUSS', true, 'focussOption', focussOpt, ...
    'lambdaTau', 0.9, 'MaxIter',500, 'sMin', sMin, 'sMax', K, ...
    'IsPrune', false, 'LocalTrap', 1e-1, 'ThresReplace', ratioR*F*K/L, ...
    'NumCore', NumCore, 'IsCut', false, 'LocalIter', 20, ...
    'ThresPrune', 0.5);
facOpt = facOption('ksvdOption', ksvdOpt, 'admmOption', admmOpt, ...
    'PriorConst', PriorConst);
[Pihat, Bhat, results] = facMethod(W_clean, L, facOpt);

T = computeT(B', Bhat');
figure;
imagesc(T);
errB = norm(Bhat - T'*B, 'fro');
errPi = norm(Pihat*T' - Pi, 'fro');
fprintf('The representation error of W is %.4e\n', ...
    norm(W_clean - Pihat*Bhat, 'fro'));
fprintf('The representation error of B is %.4e\n', errB);
fprintf('The representation error of Pi is %.4e\n', errPi);

file = fullfile(Dir, 'factorization.mat');
save(file);
