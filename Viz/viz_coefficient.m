function viz_coefficient(K, L)
%Dir = '~/GitHub/Chen-CVPR-2016/data/Compressibility';
Dir = '~/GitHub/Chen-CVPR-2016/data/SampleCompressibility';
file = fullfile(Dir, sprintf('K%dL%02d.mat', K, L));
load(file);
figure;
imagesc(abs(Xhat));
alpha 0.6
%colormap 'hot';
xlabel('Frames', 'FontWeight', 'bold', ...
    'FontName', 'New Roman Times');
ylabel('Bases', 'FontWeight', 'bold', ...
    'FontName', 'New Roman Times');
