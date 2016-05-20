% function loads CMU-MoCap data, processes data and projects it to images.
function Data = load_CMU_MoCap(subject, sequences, sampleF)
if nargin < 3
    sampleF = 1; % Sample frequency.
end
globalPath;
Ss = cell(numel(sequences), 1);
for i = 1:numel(sequences)
    seq = sequences(i);
    file = fullfile(cmuMoCapDir, 'aux', sprintf('%02d', subject), ...
        sprintf('%02d_%02d', subject, seq));
    aux = load(file);
    Ss{i} = cat(1, aux.all_points.pts);
end
rawS = cat(1, Ss{:});
idx = [1:3*sampleF:size(rawS, 1); 2:3*sampleF:size(rawS, 1)+1;...
    3:3*sampleF:size(rawS, 1)+2];
idx = reshape(idx, numel(idx), 1);
rawS = rawS(idx, :);

s = mean(std(rawS, 1, 2));
S = rawS/s;
sm = mean(S,2);
S = S - sm*ones(1,size(S,2));

F = size(S, 1)/3;
p = size(S, 2);

Ri = camera_setup(F, 'vertical');
R = blkdiag(Ri{:});

W = R*S;

Ssharp = zeros(F, 3*p);
Ssharp(:, 1:3:end) = S(1:3:end, :);
Ssharp(:, 2:3:end) = S(2:3:end, :);
Ssharp(:, 3:3:end) = S(3:3:end, :);

Data.gt.S = S;
Data.gt.Ri = Ri;
Data.gt.R = R;
Data.gt.rawS = rawS;
Data.gt.Ssharp = Ssharp;

Data.F = F;
Data.p = p;
Data.W = W;
