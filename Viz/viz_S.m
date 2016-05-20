% Function for visualizing 3D video.
function viz_S(S, varargin)
ip = inputParser;

defaultPauseTime = 0.08;
defaultGT = [];
defaultR0 = [];
defaultDataset = '';
defaultMarkerSize = 10;
defaultIsTitle = true;
defaultBaseLine = [];

addOptional(ip, 'PauseTime', defaultPauseTime, @isnumeric);
addOptional(ip, 'GT', defaultGT, @isnumeric);
addOptional(ip, 'R0', defaultR0, @isnumeric);
addOptional(ip, 'Dataset', defaultDataset, @isstr);
addOptional(ip, 'MarkerSize', defaultMarkerSize, @isnumeric);
addOptional(ip, 'IsTitle', defaultIsTitle, @islogical);
addOptional(ip, 'BL', defaultBaseLine, @isnumeric);

parse(ip, varargin{:});
PauseTime = ip.Results.PauseTime;
GT = ip.Results.GT;
R0 = ip.Results.R0;
Dataset = ip.Results.Dataset;
MarkerSize = ip.Results.MarkerSize;
IsTitle = ip.Results.IsTitle;
BL = ip.Results.BL;

if isempty(R0)
    if isempty(GT)
        R0 = eye(3);
    else
        R0 = findRotation(GT, S);
    end
end

switch Dataset
case 'MoCap'
    ConnectPts = [1,2; 2,3; 3,4; 4,5; 5,6; 1,7; 7,8; 8,9; 9,10; 10,11; 1,12; ...
        12,13; 13,14; 14,15; 15,16; 16,17; 15,18; 18,19; 19,20; 20, 21; 21,22; ...
        22,23; 23,24; 15,25; 25,26; 26,27; 27,28; 28,29; 29,30; 30,31];
otherwise
    ConnectPts = [];
end
% figure;
maxC = 1.05*max(R0'*[reshape(S, 3, numel(S)/3), reshape(GT, 3, numel(GT)/3)], [], 2);
minC = 1.05*min(R0'*[reshape(S, 3, numel(S)/3), reshape(GT, 3, numel(GT)/3)], [], 2);
for i = 1:size(S, 1)/3
    Sv(3*(i-1)+1:3*(i-1)+3, :) = R0'*S(3*(i-1)+1:3*(i-1)+3, :);
    plot3(Sv(3*(i-1)+1, :), Sv(3*(i-1)+2, :), Sv(3*(i-1)+3, :), 'r.', ...
        'markersize', MarkerSize*2);
    hold on;
    if ~isempty(ConnectPts)
        for j = 1:size(ConnectPts, 1)
            k1 = ConnectPts(j, 1);
            k2 = ConnectPts(j, 2);
            line([Sv(3*(i-1)+1, k1); Sv(3*(i-1)+1, k2)], [Sv(3*(i-1)+2, k1); ...
                Sv(3*(i-1)+2, k2)], [Sv(3*(i-1)+3, k1); Sv(3*(i-1)+3, k2)]);
        end
    end
    if ~isempty(BL)
        plot3(BL(3*(i-1)+1, :), BL(3*(i-1)+2, :), BL(3*(i-1)+3, :), 'g.', ...
            'markersize', MarkerSize*2);
        if ~isempty(ConnectPts)
            for j = 1:size(ConnectPts, 1)
                k1 = ConnectPts(j, 1);
                k2 = ConnectPts(j, 2);
                line([BL(3*(i-1)+1, k1); BL(3*(i-1)+1, k2)], [BL(3*(i-1)+2, k1); ...
                    BL(3*(i-1)+2, k2)], [BL(3*(i-1)+3, k1); BL(3*(i-1)+3,k2)],...
                    'Color', 'g');
            end
        end
    end
    axis equal;
    axis off;
    axis([minC(1),maxC(1), minC(2),maxC(2), minC(3),maxC(3)]);
    if IsTitle
        title(sprintf('3D    f = %d', i));
    end
%     line([0;maxC(1)], [0;0], [0;0]);
%     line([0;0], [0;maxC(2)], [0;0]);
%     line([0;0], [0;0], [0;maxC(3)]);
    if ~isempty(GT)
        plot3(GT(3*(i-1)+1, :), GT(3*(i-1)+2, :), GT(3*(i-1)+3, :), 'bo', ...
            'markersize', MarkerSize);
        for j = 1:size(S, 2)
            line([GT(3*(i-1)+1, j); Sv(3*(i-1)+1, j)], [GT(3*(i-1)+2, j); ...
                Sv(3*(i-1)+2, j)], [GT(3*(i-1)+3, j); Sv(3*(i-1)+3, j)]);
        end
    end
    hold off;
    if PauseTime == 0
        pause
    else
        pause(PauseTime);
    end
end 
