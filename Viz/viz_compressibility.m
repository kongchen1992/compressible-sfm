% This script is used for generating figures for paper
Dir = '~/GitHub/Chen-CVPR-2016/data/Compressibility';
%Ks = 2:7;
%Ls = 10:5:45;
Ks = [3, 4];
Ls = [30, 40];
errS = zeros(numel(Ks), numel(Ls));
cohereceS = zeros(numel(Ks), numel(Ls));
for i = 1:numel(Ks)
    K = Ks(i);
    for j = 1:numel(Ls)
        L = Ls(j);
        file = fullfile(Dir, sprintf('K%dL%02d_3S.mat', K, L));
        data = load(file);
        %errS(i, j) = mean(data.errS);
        errS(i, j) = norm(data.Ssharp - data.Ssharphat, 'fro')...
            /norm(data.Ssharp, 'fro');
        coherenceS(i, j) = coherence(data.Dksvd);
    end
end

[X, Y] = meshgrid(Ls, Ks);
[Xq, Yq] = meshgrid(min(Ls):0.1:max(Ls), min(Ks):0.1:max(Ks));
Vq = interp2(X, Y, errS, Xq, Yq, 'cubic');

figure;
surf(Xq, Yq, Vq, 'EdgeColor', 'none');
xlabel('L', 'FontSize', 20, 'FontWeight', 'bold', ...
    'FontName', 'New Roman Times');
ylabel('K', 'FontSize', 20, 'FontWeight', 'bold', ...
    'FontName', 'New Roman Times');
zlabel('Represent Error', 'FontSize', 20, 'FontWeight', 'bold', ...
    'FontName', 'New Roman Times');
set(gca, 'fontSize', 20);

Vq = interp2(X, Y, coherenceS, Xq, Yq, 'cubic');

figure;
%surf(X, Y, errS);
surf(Xq, Yq, Vq, 'EdgeColor', 'none');
%title('Coherence vs. K and L');
xlabel('L', 'FontSize', 20, 'FontWeight', 'bold', ...
    'FontName', 'New Roman Times');
ylabel('K', 'FontSize', 20, 'FontWeight', 'bold', ...
    'FontName', 'New Roman Times');
zlabel('Coherence', 'FontSize', 20, 'FontWeight', 'bold', ...
    'FontName', 'New Roman Times');
%set(gca, 'ztick', []);
