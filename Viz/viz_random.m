function [errR, errR_dai, errS, errS_dai] = viz_random(Is)
Dir = '~/GitHub/Chen-CVPR-2016/data/RANDOM';
rand_file = fullfile(Dir, 'random.mat');
K = 2;
Ls = 3:12;
F = 100;
P = 30;
if exist(rand_file, 'file')
    load(rand_file);
else
    errR = zeros(numel(Ls), numel(Is));
    errR_dai = zeros(numel(Ls), numel(Is));
    errS = zeros(numel(Ls), numel(Is));
    errS_dai = zeros(numel(Ls), numel(Is));
    for l = 1:numel(Ls)
        L = Ls(l);
        for i = 1:numel(Is);
            I = Is(i);
            file = fullfile(Dir, sprintf('L%03dF%04dP%03d_%02d.mat', L, F, P, I));
            data = load(file);
            errR(l, i) = mean(data.errR);
            errS(l, i) = mean(data.errS);

            file = fullfile(Dir, sprintf('dai_L%03dF%04dP%03d_%02d.mat', L, F, P, I));
            data =load(file);
            errR_dai(l, i) = mean(data.Rotation_Err);
            errS_dai(l, i) = mean(data.Shape_Err_BMM);
        end
    end
    save(rand_file, 'errR', 'errS', 'errR_dai', 'errS_dai');
end

errR = mean(errR, 2);
errR_dai = mean(errR_dai, 2);
errS = mean(errS, 2);
errS_dai = mean(errS_dai, 2);

minX = min(Ls);
maxX = max(Ls);
minY = -0.1;
maxY = max([errR; errR_dai])+0.1;
figure;
hold on;
plot(Ls, [errR, errR_dai], '+', 'MarkerSize', 10, 'LineWidth', 2);
f = fit(Ls', errR, 'poly3');
f_dai = fit(Ls', errR_dai, 'poly3');
h = plot(f, 'b:');
set(h, 'LineWidth', 2);
h = plot(f_dai, 'r:');
set(h, 'LineWidth', 2);
axis([minX, maxX, minY, maxY]);
xlabel('L', 'FontSize', 20, 'FontWeight', 'Bold', ...
    'FontName', 'Times New Roman')
ylabel('Error of Recovered Cameras', 'FontSize', 20, 'FontWeight', 'Bold', ...
    'FontName', 'Times New Roman');
legend({'The proposed Method', 'Dai. et al'}, 'Location', 'northwest', ...
    'FontSize',20, 'FontName', 'Times New Roman');
grid on;

minX = min(Ls);
maxX = max(Ls);
minY = -0.1;
maxY = max([errS; errS_dai])+0.1;
figure;
hold on;
plot(Ls, [errS, errS_dai], '+', 'MarkerSize', 10, 'LineWidth', 2);
f = fit(Ls', errS, 'poly3');
f_dai = fit(Ls', errS_dai, 'poly3');
h = plot(f, 'b:');
set(h, 'LineWidth', 2);
h = plot(f_dai, 'r:');
set(h, 'LineWidth', 2);
axis([minX, maxX, minY, maxY]);
xlabel('L', 'FontSize', 20, 'FontWeight', 'Bold', ...
    'FontName', 'Times New Roman')
ylabel('Error of Recovered Structures', 'FontSize', 20, 'FontWeight', 'Bold', ...
    'FontName', 'Times New Roman');
legend({'The proposed Method', 'Dai. et al'}, 'Location', 'northwest', ...
    'FontSize',20, 'FontName', 'Times New Roman');
grid on;
