function exp_selectFrames(configI)
Dir = '~/GitHub/Chen-CVPR-2016/data/Structures';
file = fullfile(Dir, sprintf('config%d.mat', configI));
data = load(file);
K = data.K;
L = data.L;

Csharp = data.Xhat';
Ssharp = data.Ssharp;
UtiAtom = sum(Csharp~=0, 1);
[UA, UI] = sort(UtiAtom, 'descend');

maxUti = round(1.2*K*600/L);
UA = UA(UA>maxUti);
UI = UI(1:numel(UA));

for i = 1:numel(UI)
    I = UI(i);
    numDel = UA(i) - maxUti;
    for j = i+1:numel(UI);
        if numDel == 0
            break;
        end
        J = UI(j);
        indx_toDel = (sum(Csharp(:, [I, J])~=0, 2) == 2);
        indxtoDel = find(indx_toDel);
        num_toDel = numel(indxtoDel);
        if num_toDel > numDel
            indx = randsample(indxtoDel, numDel);
            indx_toDel = zeros(size(indx_toDel));
            indx_toDel(indx) = 1;
        end
        num_toDel = numel(find(indx_toDel));
        Ssharp = Ssharp(~indx_toDel, :);
        Csharp = Csharp(~indx_toDel, :);
        numDel = numDel - num_toDel;
        UA(i) = UA(i) - num_toDel;
        UA(j) = UA(j) - num_toDel;
    end
end

Bsharp = data.Dhat';
F = size(Ssharp, 1);

UtiAtom = sum(Csharp~=0, 1);
idx = find(UtiAtom < 0.5*K*F/L, 1);
while ~isempty(idx)
    indx_toDel = (Csharp(:, idx) ~= 0);
    Csharp = Csharp(~indx_toDel, [1:idx-1, idx+1:end]);
    Ssharp = Ssharp(~indx_toDel, :);
    Bsharp = Bsharp([1:idx-1, idx+1:end], :);
    UtiAtom = sum(Csharp~=0, 1);
    F = size(Ssharp, 1);
    L = size(Bsharp, 1);
    idx = find(UtiAtom < 0.5*K*F/L, 1);
end

P = data.P;
file = fullfile(Dir, sprintf('config%d_Structures.mat', configI));
save(file, 'Csharp', 'Bsharp', 'Ssharp', 'P', 'K', 'L', 'F');

fprintf('New F = %d\n', F);
UA = sum(Csharp~=0, 1);
fprintf('%d ', UA);
fprintf('\n%.4f\n', min(UA)/(K*F/L));
