function viz_results(configI, I)
Dir = '~/GitHub/Chen-CVPR-2016/data/Structures';
results_file = fullfile(Dir, sprintf('config%d_results.mat', configI));
data = load(results_file);
dai_file = fullfile(Dir, sprintf('config%d_dai.mat', configI));
data_dai = load(dai_file);
cleanStructure_file = fullfile(Dir, sprintf('config%d_cleanStructures.mat', configI));
data_structure = load(cleanStructure_file);
S = data_structure.S_clean;
Srot = data.Srot;
Sdai = data_dai.Shat_BMM;
[~, Srot_dai] = compareStructs(S, Sdai);

fprintf('The error of proposed method is %.4e\n', mean(data.errS));
fprintf('The error of Dai et.al. is %.4e\n', data_dai.Shape_Err_BMM);

viz_S(Srot(:, 31*I-30:31*I), 'GT', S(:, 31*I-30:31*I), 'R0', eye(3), 'BL', Srot_dai(:, 31*I-30:31*I), 'Dataset', 'MoCap', ...
    'PauseTime', 0, 'MarkerSize', 20);
