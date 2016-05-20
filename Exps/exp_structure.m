% This script is for experiments on compressibility of a certain video
function exp_structure(configI)
configFile = 'Exps/config_real.txt';
fileID = fopen(configFile);
Config = textscan(fileID, '%d-%d-%d', 'CommentStyle', '#', 'Delimiter', '\n');
Data{1} = load_CMU_MoCap(Config{1}(3*configI-2), Config{2}(3*configI-2), 1);
Data{2} = load_CMU_MoCap(Config{1}(3*configI-1), Config{2}(3*configI-1), 1);
Data{3} = load_CMU_MoCap(Config{1}(3*configI-0), Config{2}(3*configI-0), 1);

F = Config{3}(3*configI-2);
K = 2;
L = 20;

indx = cell(3, 1);
if F < Data{1}.F && F < Data{2}.F && F < Data{3}.F
    IsRepeat = false;
    disp('Not Repeat');
else
    IsRepeat = true;
    disp('Repeat');
end
for i = 1:3
    indx{i} = randsample(Data{i}.F, F, IsRepeat);
end
Ssharp = [Data{1}.gt.Ssharp(indx{1}, :), Data{2}.gt.Ssharp(indx{2}, :), ...
    Data{3}.gt.Ssharp(indx{3}, :)];
P = size(Ssharp, 2)/3;
S = zeros(3*F, P);
S(1:3:end, :) = Ssharp(:, 1:3:end);
S(2:3:end, :) = Ssharp(:, 2:3:end);
S(3:3:end, :) = Ssharp(:, 3:3:end);

Sprior = diag(sqrt(sum(Ssharp.^2, 1)));
SsharpP = Ssharp/Sprior;

admmOpt = admmOption('MaxIter', 5000);
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

Dir = '~/GitHub/Chen-CVPR-2016/data/Structures';
if ~exist(Dir, 'dir');
    mkdir(Dir);
end
file = fullfile(Dir, sprintf('config%d.mat', configI));
save(file, 'Xhat', 'Dhat', 'K', 'L', 'P', 'Ssharp');
