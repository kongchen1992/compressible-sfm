function viz_noise
Ps = [1, 3:5, 7:8, 10, 12:23];
%Ps = 1:23;
Is = 1:5;
Dir = '~/GitHub/Chen-CVPR-2016/data/Noise';
for i = 1:numel(Is)
    I = Is(i);
    file = fullfile(Dir, sprintf('noise_%02d.mat', I));
    dai_file = fullfile(Dir, sprintf('noise_dai_%02d.mat', I));

    if exist(file, 'file')
        load(file);
        ErrR{i} = arrayfun(@(x) mean(x.errR), Results);
        ErrS{i} = arrayfun(@(x) mean(x.errS), Results);
    end

    if exist(dai_file, 'file')
        load(dai_file);
        ErrR_dai{i} = arrayfun(@(x) x.Rotation_Err, Results_dai);
        ErrS_dai{i} = arrayfun(@(x) x.Shape_Err_BMM, Results_dai);
    end
end
errR = mean(cat(1, ErrR{:}), 1)';
errR_dai = mean(cat(1, ErrR_dai{:}), 1)';
errS = mean(cat(1, ErrS{:}), 1)';
errS_dai = mean(cat(1, ErrS_dai{:}), 1)';

%ratios = -5:0.1:-0.5;
ratios = -5:0.2:-0.5;

ratios = ratios(Ps);
errR = log(errR(Ps))/log(10);
errR_dai = log(errR_dai(Ps))/log(10);
errS = log(errS(Ps))/log(10);
errS_dai = log(errS_dai(Ps))/log(10);
%ratios = ratios(Ps);
%errR = errR(Ps);
%errR_dai = errR_dai(Ps);
%errS = errS(Ps);
%errS_dai = errS_dai(Ps);

figure;
hold on;
minX = min(ratios);
maxX = max(ratios);
minY = min([errR; errR_dai])-1;
maxY = max([errR; errR_dai])+1;
%minY = min([errR; errR_dai])-0.1;
%maxY = max([errR; errR_dai])+0.1;
plot(ratios, [errR, errR_dai], 'x', 'MarkerSize', 10, 'LineWidth', 2);
axis([minX, maxX, minY, maxY]);
%f = fit(ratios', errR, 'poly4');
%f_dai = fit(ratios', errR_dai, 'poly4');
%h = plot(f, 'b:');
%set(h, 'LineWidth', 2);
%h = plot(f_dai, 'r:');
%set(h, 'LineWidth', 2);
[param,stat]=sigm_fit(ratios,errR_dai,[],[],true, 'r:', ...
    'MarkerSize', 10, 'LineWidth', 2);
[param,stat]=sigm_fit(ratios,errR,[],[],true, 'b:', ...
    'MarkerSize', 10, 'LineWidth', 2);
xlabel('Noise Ratio', 'FontSize', 20, 'FontWeight', 'Bold', ...
    'FontName', 'Times New Roman');
ylabel('Error of Recovered Cameras', 'FontSize', 20, 'FontWeight', 'Bold', ...
    'FontName', 'Times New Roman');
legend({'The proposed Method', 'Dai. et al'}, 'Location', 'southeast', ...
    'FontSize',20, 'FontName', 'Times New Roman');
grid on;

figure;
hold on;
minX = min(ratios);
maxX = max(ratios);
minY = min([errS; errS_dai])-1;
maxY = max([errS; errS_dai])+1;
%minY = min([errR; errR_dai])-0.1;
%maxY = max([errR; errR_dai])+0.1;
plot(ratios, [errS, errS_dai], 'x', 'MarkerSize', 10, 'LineWidth', 2);
axis([minX, maxX, minY, maxY]);
%f = fit(ratios', errR, 'poly4');
%f_dai = fit(ratios', errR_dai, 'poly4');
%h = plot(f, 'b:');
%set(h, 'LineWidth', 2);
%h = plot(f_dai, 'r:');
%set(h, 'LineWidth', 2);
[param,stat]=sigm_fit(ratios,errS_dai,[],[],true, 'r:', ...
    'MarkerSize', 10, 'LineWidth', 2);
[param,stat]=sigm_fit(ratios,errS,[],[],true, 'b:', ...
    'MarkerSize', 10, 'LineWidth', 2);
xlabel('Noise Ratio', 'FontSize', 20, 'FontWeight', 'Bold', ...
    'FontName', 'Times New Roman')
ylabel('Error of Recovered Structures', 'FontSize', 20, 'FontWeight', 'Bold', ...
    'FontName', 'Times New Roman');
legend({'The proposed Method', 'Dai. et al'}, 'Location', 'southeast', ...
    'FontSize',20, 'FontName', 'Times New Roman');
grid on;
