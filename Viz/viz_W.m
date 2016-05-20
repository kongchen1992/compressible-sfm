% Function for visualizing images series.

function viz_W(W, GT)
if nargin < 2
    GT = [];
end
% figure;
if isempty(GT)
    mm = 1.1*minmax(reshape(W, 2, numel(W)/2));
else
    mm = 1.1*minmax([reshape(W, 2, numel(W)/2), reshape(GT, 2, numel(GT)/2)]);
end
for i = 1:size(W, 1)/2
    plot(W(2*(i-1)+1, :), W(2*(i-1)+2, :), 'r.', 'markersize', 8);
    if ~isempty(GT)
        hold on;
        plot(GT(2*(i-1)+1, :), GT(2*(i-1)+2, :), 'bo', 'markersize', 8);
        for j = 1:size(W, 2)
            line([GT(2*(i-1)+1, j); W(2*(i-1)+1, j)], [GT(2*(i-1)+2, j); ...
                W(2*(i-1)+2, j)])
        end
        hold off;
    end
    axis equal
    axis(reshape(mm', 1, 4));
    title(sprintf('2D    f = %d', i));
    pause(0.05)
end