% Title
function exp_mocap(configI, configRank)
configFile = sprintf('Exps/config_mocap_%s.txt', configRank);
fileID = fopen(configFile);
Config = textscan(fileID, '%d-%d-%d', 'CommentStyle', '#', 'Delimiter', '\n');
Data{1} = load_CMU_MoCap(Config{1}(3*configI-2), Config{2}(3*configI-2), ...
    Config{3}(3*configI-2));
Data{2} = load_CMU_MoCap(Config{1}(3*configI-1), Config{2}(3*configI-1), ...
    Config{3}(3*configI-1));
Data{3} = load_CMU_MoCap(Config{1}(3*configI-0), Config{2}(3*configI-0), ...
    Config{3}(3*configI-0));

switch configRank
case 'full'
    K = 33;
    PriorConst = 3.5;
case 'high'
    K = 20;
    PriorConst = 4;
case 'low'
    K = 10;
    PriorConst = 5;
end

Bsharp = [Data{1}.gt.Ssharp(1:K, :), Data{2}.gt.Ssharp(1:K, :), ...
    Data{3}.gt.Ssharp(1:K, :)];
M = size(Bsharp, 2)/3;
s = 2;
F = 500;

Bclean = zeros(3*K, M);
Bclean(1:3:end, :) = Bsharp(:, 1:3:end);
Bclean(2:3:end, :) = Bsharp(:, 2:3:end);
Bclean(3:3:end, :) = Bsharp(:, 3:3:end);

% add some noise into dictionary for better condition
B = Bclean + 0.1*std(Bclean(:))*randn(size(Bclean));
%Wprior = diag(sqrt(sum(B.^2, 1)));

Csharp = zeros(F, K);
for i = 1:F
    sp = randperm(K, s);
    Csharp(i, sp(1)) = rand(1, 1);
    Csharp(i, sp(2)) = 1 - Csharp(i, sp(1));
end
Ri = camera_setup(F, 'continuous');
R = blkdiag(Ri{:});
C = kron(Csharp, eye(3));
Pi = R * C;
S = C * B;
W = R * S;

% Get number of processors on current machine
[~, nproc] = unix('nproc');

admmOpt = admmOption('MaxIter', 5000);
focussOpt = focussOption('p', 3, 'Trunc', true, 'Tikhonov', true, ...
    'lambda', 6e-3, 'BlurC', 1);
ksvdOpt = ksvdOption('IsFOCUSS', true, 'focussOption', focussOpt, ...
    'lambdaTau', 0.9, 'MaxIter',500, 'sMin', 1, 'sMax', s, 'IsPrune', false, ...
    'LocalTrap', 1e-1, 'ThresReplace', 0.8*F*2/K, 'NumCore', str2num(nproc)/2, ...
    'IsCut', false);
facOpt = facOption('ksvdOption', ksvdOpt, 'admmOption', admmOpt, ...
    'PriorConst', PriorConst);
[Pihat, Bhat, results] = facMethod(W, K, facOpt);

%% Solve Camera Matrix and 3D Structure
%[Rhat, Shat, Rihat] = solveProblem(Pihat, Bhat, [], 0); % NRSfM
%
%% Evaluate Results
%errR = compareRotations(cat(1, Ri{:}), cat(1, Rihat{:}));
%[errS, Srot] = compareStructs(S, Shat);
%fprintf('\tThe error in camera matrix is %.2e\n', mean(errR));
%fprintf('\tThe error in shape matrix is %.2e\n', mean(errS));
ResultPath;
dir = fullfile(MOCAPDIR, sprintf('%sRank', configRank));
if ~exist(dir, 'dir')
    mkdir(dir);
end
save(fullfile(dir, sprintf('mocap_%02d.mat', configI)));
