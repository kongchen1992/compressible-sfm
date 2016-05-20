function exp_real(configI)
Dir = '~/GitHub/Chen-CVPR-2016/data/Structures';

structure_file = fullfile(Dir, sprintf('config%d.mat', configI));
if ~exist(structure_file, 'file')
    disp('Decomposing Structures...');

    configFile = 'Exps/config_real.txt';
    fileID = fopen(configFile);
    Config = textscan(fileID, '%d-%d-%d', 'CommentStyle', '#', 'Delimiter', '\n');
    Data{1} = load_CMU_MoCap(Config{1}(3*configI-2), Config{2}(3*configI-2), 1);
    Data{2} = load_CMU_MoCap(Config{1}(3*configI-1), Config{2}(3*configI-1), 1);
    Data{3} = load_CMU_MoCap(Config{1}(3*configI-0), Config{2}(3*configI-0), 1);

    F = Config{3}(3*configI-2);
    K = 2;
    L = 100;

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
    indDhat = randsample(F, L);
    Dhat = Ssharp(indDhat, :)';
    normedDhat = Dhat*inv(diag(sqrt(sum(Dhat.^2, 1))));
    Inner = abs(normedDhat'*normedDhat);
    Inner = triu(Inner, 1);
    [~, IInner] = sort(Inner(:), 'descend');
    [I, J] = ind2sub(size(Inner), IInner(1:30));
    indxL = unique([I;J]);
    L = 25;
    Dhat = Dhat(:, indxL(1:L));
    Xhat = zeros(L, F);

    % Use OMP to re-generate new Csharp
    for i = 1:F
        Xhat(:, i) = omp(Dhat, Ssharp(i, :)', K, [1, 1], []);
    end

    save(structure_file, 'Xhat', 'Dhat', 'K', 'L', 'P', 'Ssharp');
else
    disp('Structures Composition Existed');
end

selectFrames_file = fullfile(Dir, sprintf('config%d_Structures.mat', configI));
if ~exist(selectFrames_file, 'file')
    disp('Selecting Frames...');
    exp_selectFrames(configI);
else
    disp('Frames have been selected already');
end

cleanStructure_file = fullfile(Dir, sprintf('config%d_cleanStructures.mat', configI));
if ~exist(cleanStructure_file, 'file')
    disp('Generating clean Structures...');
    data = load(selectFrames_file);
    Bsharp = data.Bsharp;
    Bprior = diag(sqrt(sum(Bsharp.^2, 2)));
    Bsharp = inv(Bprior)*Bsharp;
    Csharp = data.Csharp;
    Csharp = Csharp*Bprior;
    L = data.L;
    K = data.K;
    F = data.F;
    P = data.P;
    Ssharp = data.Ssharp;

    Ssharp_clean = Csharp*Bsharp;
    err = Ssharp - Ssharp_clean;
    fprintf('Norm of Err is %.2f\n', norm(err, 'fro'));
    fprintf('Coherence of dictionary is %.4f\n', coherence(Bsharp'));
    fprintf('Condition of dictionary is %.4f\n', cond(Bsharp));

    S_clean = zeros(3*F, P);
    S_clean(1:3:end, :) = Ssharp_clean(:, 1:3:end);
    S_clean(2:3:end, :) = Ssharp_clean(:, 2:3:end);
    S_clean(3:3:end, :) = Ssharp_clean(:, 3:3:end);

    Ri = camera_setup(F, 'continuous');
    R = blkdiag(Ri{:});
    W_clean = R*S_clean;
    save(cleanStructure_file, 'Bsharp', 'L', 'K', 'F', 'P', 'Ssharp_clean', ...
    'Ssharp', 'W_clean', 'S_clean', 'Ri', 'R', 'Csharp');
else
    disp('Clean Structures Existed');
    data = load(cleanStructure_file);
    S = data.S_clean;
    Ri = data.Ri;
end

fac_file = fullfile(Dir, sprintf('config%d_factorization.mat', configI));
if ~exist(fac_file, 'file')
    disp('Factorizing clean W...');
    [Pihat, Bhat] = exp_factorize(configI);
else
    disp('Factorization Existed');
    data = load(fac_file);
    Pihat = data.Pihat;
    Bhat = data.Bhat;
end

disp('Recovering Cameras and Structures from factorization...');
[~, Shat, Rihat, ~] = solveProblem(Pihat, Bhat, [], 0);

disp('Evaluating Results...');
[errR, errS, Srot] = evalResults(Rihat, Shat, Ri, S);

results_file = fullfile(Dir, sprintf('config%d_results.mat', configI));
save(results_file, 'Shat', 'Rihat', 'errR', 'errS', 'Srot');
