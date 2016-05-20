function exp_noise_dai_I(ratios, I)
Dir = '~/GitHub/Chen-CVPR-2016/data/Noise';
file = fullfile(Dir, sprintf('noise_%02d.mat', I));
load(file);
W = GroundTruth.W;

Rs = cat(1, GroundTruth.Ri{:});
rotStruct = 0;
S_GT = GroundTruth.S;
Ls = 3:5;
Rotation_Err_L = cell(numel(Ls), 1);
Shape_Err_BMM_L = cell(numel(Ls), 1);

for i = 1:numel(ratios)
    [~, NumProc] = unix('nproc');
    NumCore = min(str2num(NumProc)/2, numel(Ls));
    poolobj = parpool('local', NumCore);

    ratio = ratios(i);
    fprintf('RATIO of NOISE is %.4f\n', ratio);
    W_noise = W + (10^ratio)*norm(W, 'fro')/norm(noise, 'fro')*noise;

    parfor l = 1:numel(Ls)
        L_dai = Ls(l);
        [~, ~, ~, ~, ~, Shape_Err_BMM_L{l}, ~, ~, Rotation_Err_L{l}] = ...
            NRSFM_BMM(W_noise, L_dai, rotStruct, S_GT, Rs);
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

    fprintf('\tThe error in camera matrix is %.2e\n', mean(Rotation_Err));
    fprintf('\tThe BMM error in shape matrix is %.2e\n', mean(Shape_Err_BMM));
    fprintf('\tThe best L_dai_R is %d for L = %d\n', bestL_R, 5);
    fprintf('\tThe best L_dai_S is %d for L = %d\n', bestL_S, 5);

    Results_dai(i).Rotation_Err = Rotation_Err;
    Results_dai(i).Shape_Err_BMM = Shape_Err_BMM;
end

file = fullfile(Dir, sprintf('noise_dai_%02d.mat', I));
save(file, 'Results_dai');
