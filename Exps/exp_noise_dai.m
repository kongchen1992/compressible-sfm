function exp_noise_dai(Is)
ratios = -5:0.2:-0.5;
for i = 1:numel(Is)
    I = Is(i);
    exp_noise_dai_I(ratios, I);
end
