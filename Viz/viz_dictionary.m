function viz_dictionary(K, L)
Dir = '~/GitHub/Chen-CVPR-2016/data/Compressibility';
file = fullfile(Dir, sprintf('K%dL%02d.mat', K, L));
load(file);
Bsharp = Dhat';
B = zeros(3*L, P);
B(1:3:end, :) = Bsharp(:, 1:3:end);
B(2:3:end, :) = Bsharp(:, 2:3:end);
B(3:3:end, :) = Bsharp(:, 3:3:end);

viz_S(B, 'Dataset', 'MoCap', 'PauseTime', 0);
