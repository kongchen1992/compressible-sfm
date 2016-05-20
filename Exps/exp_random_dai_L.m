function exp_random_dai_L(K, L, F, P, I)
Dir = '~/GitHub/Chen-CVPR-2016/data/RANDOM';
file = fullfile(Dir, sprintf('L%03dF%04dP%03d_%02d.mat', L, F, P, I));
data=load(file);

tic;

Rs = cat(1, data.Ri{:});
rotStruct = 0;
S_GT = data.S;
Ls = 2:min([size(data.S, 2)/3, L, 20]);
Rotation_Err_L = cell(numel(Ls), 1);
Shape_Err_BMM_L = cell(numel(Ls), 1);

[~, NumProc] = unix('nproc');
NumCore = min(str2num(NumProc)/2, numel(Ls));
poolobj = parpool('local', NumCore);

parfor l = 1:numel(Ls)
    L_dai = Ls(l);
    [~, ~, ~, ~, ~, Shape_Err_BMM_L{l}, ~, ~, Rotation_Err_L{l}] = ...
        NRSFM_BMM(data.W, L_dai, rotStruct, S_GT, Rs);
end

delete(gcp);

Shape_Err = cellfun(@mean, Shape_Err_BMM_L);
[~, ind] = min(Shape_Err);
bestL_S = Ls(ind);
Shape_Err_BMM = Shape_Err_BMM_L{ind};

R_Err = cellfun(@mean, Rotation_Err_L);
[~, ind] = min(R_Err);
bestL_R = Ls(ind);
Rotation_Err = Rotation_Err_L{ind};

time = toc;

fprintf('\tThe error in camera matrix is %.2e\n', mean(Rotation_Err));
fprintf('\tThe BMM error in shape matrix is %.2e\n', mean(Shape_Err_BMM));
fprintf('\tThe best L_dai_R is %d for L = %d\n', bestL_R, L);
fprintf('\tThe best L_dai_S is %d for L = %d\n', bestL_S, L);
file = fullfile(Dir, sprintf('dai_L%03dF%04dP%03d_%02d.mat', L, F, P, I));
save(file);
