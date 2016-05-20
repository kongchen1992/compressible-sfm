function g = getGroupIdx(sizeX, sizeG)
% g = GetGroupIdx(sizeX, sizeG) seperate Matrix X into several group G.
% 
%   Inpute: -sizeX is dimension of data
%           -sizeG is dimension of grouping of data
%   
%   Outpute:    -g, g(:, i) denotes an index set corresponding to the i-th
%                   group.
%
%   Example: 
%   g = getGroupIdx([4,6], [2,3]);
%
%   for i = 1:size(g, 2)
%       a(g(:, i)) = i;
%   end
%
%   a =
%
%      1     1     1     2     2     2
%      1     1     1     2     2     2
%      3     3     3     4     4     4
%      3     3     3     4     4     4


%% no overlap-grouping
assert(mod(sizeX(1), sizeG(1)) == 0);
assert(mod(sizeX(2), sizeG(2)) == 0);

%% main
N = prod(sizeX)/prod(sizeG);
g = cell(prod(sizeG), N);
k = 0;

for i = 1:sizeX(1)/sizeG(1)
    gi = (i-1)*sizeG(1)+(1:sizeG(1));
    for j = 1:sizeX(2)/sizeG(2)
        gj = (j-1)*sizeG(2)+(1:sizeG(2));
        [Gi, Gj] = meshgrid(gi, gj);
        k = k+1;
        g{k} = sub2ind(sizeX, Gi(:), Gj(:));
    end
end
assert(k == N);

g = cat(2, g{:});