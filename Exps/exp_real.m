function exp_real(configI)
Dir = '~/GitHub/Chen-CVPR-2016/data/Structures';

structure_file = fullfile(Dir, sprintf('config%d.mat', configI));
if ~exist(structure_file, 'file')
    disp('Decomposing Structures...');
    exp_structure(configI);
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
    IsIncoherent = false;
    while ~IsIncoherent
        [IsIncoherent, S, Ri] = exp_cleanStructure(configI);
    end
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
