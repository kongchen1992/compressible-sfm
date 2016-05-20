Dir = '~/GitHub/Chen-CVPR-2016/data/jaws';
structure_file = fullfile(Dir, 'structures.mat');
if ~exist(structure_file, 'file')
    disp('Generating Structures...');
    exp_jaws_structures();
else
    disp('Structures Exist');
end

fac_file = fullfile(Dir, 'factorization.mat');
if ~exist(fac_file)
    disp('Factorizing W...');
    [Pihat, Bhat, S, Ri] = debug_jaws_factorization()
else
    disp('Factorization Exist');
    data = load(fac_file);
    S = data.S_clean;
    Ri = data.Ri;
    Pihat = data.Pihat;
    Bhat = data.Bhat;
end

disp('Recovering Cameras and Structures from factorization...');
[~, Shat, Rihat, ~] = solveProblem(Pihat, Bhat, [], 0);

disp('Evaluating Results...');
[errR, errS, Srot] = evalResults(Rihat, Shat, Ri, S);

results_file = fullfile(Dir, 'results.mat');
save(results_file, 'Shat', 'Rihat', 'errR', 'errS', 'Srot');
