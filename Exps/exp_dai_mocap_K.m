% Title
function exp_dai_mocap_K(configI, configRank, K)

ResultPath;
dir = fullfile(MOCAPDIR, sprintf('%sRank', configRank));
Data = load(fullfile(dir, sprintf('mocap_%02d.mat', configI)));

Rs = cat(1, Data.Ri{:});
rotStruct = 0;
S_GT = Data.S;

[Shat_BMM, Shat_PI, Shat_Improved, Rsh, R_Recover, Shape_Err_BMM, ...
    Shape_Err_PI, Shape_Err_Smooth, Rotation_Err] = ...
    NRSFM_BMM(Data.W, K, rotStruct, S_GT, Rs);

fprintf('\tThe error in camera matrix is %.2e\n', mean(Rotation_Err));
fprintf('\tThe BMM error in shape matrix is %.2e\n', mean(Shape_Err_BMM));
fprintf('\tThe PI error in shape matrix is %.2e\n', mean(Shape_Err_PI));
fprintf('\tThe Smooth error in shape matrix is %.2e\n', mean(Shape_Err_Smooth));
save(fullfile(dir, sprintf('dai_%02d_%02d.mat', configI, K)));
