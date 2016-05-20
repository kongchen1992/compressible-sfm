%This script is the 3rd Experiments of synthetic data.
function exp_random(Is)
K = 2;
Ls = 3:12;
F = 100;
P = 30;
for i = 1:numel(Is)
    I = Is(i);
    exp_random_I(K, Ls, F, P, I);
end
