function exp_jaws_dai
Dir = '~/GitHub/Chen-CVPR-2016/data/jaws';
fac_file = fullfile(Dir, 'factorization.mat');
data = load(fac_file);

Rs = cat(1, data.Ri{:});
rotStruct = 0;
S_GT = data.S_clean;
W = data.W_clean;
Ks = 2:data.L;

[~, NumProc] = unix('nproc');
NumCore = min(str2num(NumProc)/2, numel(Ks));
poolobj = parpool('local', NumCore);

numK = numel(Ks);
Shat_BMM_K = cell(numK, 1);
R_Recover_K = cell(numK, 1);
Shape_Err_BMM_K = cell(numK, 1);
Rotation_Err_K = cell(numK, 1);

parfor i = 1:numK
    K = Ks(i);
    [Shat_BMM_K{i}, ~, ~, ~, R_Recover_K{i}, Shape_Err_BMM_K{i}, ...
        ~, ~, Rotation_Err_K{i}] = ...
        NRSFM_BMM(W, K, rotStruct, S_GT, Rs);
end

Shape_Err = cellfun(@mean, Shape_Err_BMM_K);
[~, ind] = min(Shape_Err);
bestK_S = Ks(ind);
Shape_Err_BMM = Shape_Err_BMM_K{ind};
Shat_BMM = Shat_BMM_K{ind};

R_Err = cellfun(@mean, Rotation_Err_K);
[~, ind] = min(R_Err);
bestK_R = Ks(ind);
Rotation_Err = Rotation_Err_K{ind};
R_Recover = R_Recover_K{ind};

fprintf('\tThe error in camera matrix is %.2e\n', mean(Rotation_Err));
fprintf('\tThe BMM error in shape matrix is %.2e\n', mean(Shape_Err_BMM));
file = fullfile(Dir, sprintf('dai.mat'));
save(file, 'Shat_BMM', 'R_Recover', 'bestK_S', 'bestK_R', ...
    'Shape_Err_BMM', 'Rotation_Err');
