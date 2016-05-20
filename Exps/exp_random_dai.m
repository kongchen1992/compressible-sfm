function exp_random_dai(K, Ls, F, P, Is)
for i = 1:numel(Is)
    I = Is(i);
    for l = 1:numel(Ls)
        L = Ls(l);
        exp_random_dai_L(K, L, F, P, I);
    end
end
