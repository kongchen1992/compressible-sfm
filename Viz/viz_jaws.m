Dir = '~/GitHub/Chen-CVPR-2016/data/jaws';
fac_file = fullfile(Dir, 'factorization.mat');
fac_data = load(fac_file);
S = fac_data.S_clean;
Ri = fac_data.Ri;

results_file = fullfile(Dir, 'results.mat');
res_data = load(results_file);
Srot = res_data.Srot;

dai_file = fullfile(Dir, 'dai.mat');
dai_data = load(dai_file);
Sdai = dai_data.Shat_BMM;
%[~, Srot_dai] = compareStructs(S, Sdai);

fprintf('The error of proposed method is %.4e\n', mean(res_data.errS));
fprintf('The error of Dai et.al. is %.4e\n', dai_data.Shape_Err_BMM);

for i = 1:size(Srot, 1)/3
    viz_S(Srot(3*i-2:3*i, :), 'GT', S(3*i-2:3*i, :), 'R0', eye(3), ...
        'BL', Sdai(3*i-2:3*i, :), 'PauseTime', 0, 'MarkerSize', 20);
end
