function exp_noise(Is)
%Is = 2:5;
ratios = -5:0.2:-0.5;
for i = 1:numel(Is)
    I = Is(i);
    exp_noise_I(ratios, I);
end
