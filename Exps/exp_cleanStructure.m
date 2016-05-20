function [IsOK, S_clean, Ri] = exp_cleanStructure(configI)
[~, nproc] = unix('nproc');
NumCore = str2num(nproc)/2;
Dir = '~/GitHub/Chen-CVPR-2016/data/Structures';
file = fullfile(Dir, sprintf('config%d_Structures.mat', configI));
data = load(file);

P = data.P;
K = data.K;
L = data.L;
F = data.F;
Ssharp = data.Ssharp;
Bsharp = data.Bsharp;
noiseR = 0.1;
Bsharp = Bsharp + noiseR*std(Bsharp(:))*randn(size(Bsharp));

Bprior = inv(diag(sqrt(sum(Bsharp.^2, 2))));
Bsharp = Bprior*Bsharp;

Csharp = zeros(size(data.Csharp));
for i = 1:F
    Csharp(i, randsample(L, K)) = randn(K, 1);
end

Ssharp_clean = Csharp*Bsharp;
err = Ssharp - Ssharp_clean;
fprintf('Norm of Err is %.2f\n', norm(err, 'fro'));
fprintf('Coherence of dictionary is %.2f\n', coherence(Bsharp));
fprintf('Condition of dictionary is %.2f\n', cond(Bsharp));

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

ratioR = 0.8;
sMin = 1;
PriorConst = 1.3;
Sprior = diag(sqrt(sum(S_clean.^2, 1)));

admmOpt = admmOption('MaxIter', 5000, 'gamma', 0.04, ...
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
    'PriorConst', PriorConst, 'Sprior', Sprior);
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

if errPi < 1e-2
    IsOK = true;
else
    IsOK = false;
end

% Use OMP to re-generate new Csharp
for i = 1:F
    xOMP = omp(Bsharp', Ssharp(i, :)', K, [1, 1], []);
    Csharp(i, :) = xOMP';
end

Ssharp_clean = Csharp*Bsharp;
err = Ssharp - Ssharp_clean;
fprintf('Norm of Err is %.2f\n', norm(err, 'fro'));

S_clean = zeros(3*F, P);
S_clean(1:3:end, :) = Ssharp_clean(:, 1:3:end);
S_clean(2:3:end, :) = Ssharp_clean(:, 2:3:end);
S_clean(3:3:end, :) = Ssharp_clean(:, 3:3:end);

Ri = camera_setup(F, 'continuous');
R = blkdiag(Ri{:});
W_clean = R*S_clean;

file = fullfile(Dir, sprintf('config%d_cleanStructures.mat', configI));
save(file, 'Bsharp', 'L', 'K', 'F', 'P', 'Ssharp_clean', 'Ssharp', ...
    'W_clean', 'S_clean', 'Ri', 'R', 'Csharp');
