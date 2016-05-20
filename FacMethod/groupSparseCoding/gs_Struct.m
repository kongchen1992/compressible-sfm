% Function for checking group sparse structure of matrix
%   structure = gs_Struct(Matrix, BlockSize)
%   
%   Output:
%       -structure: a matrix recording group sparse structure of input
%                   matrix M.  the value of each elements here denotes how
%                   many active elemetns in corresponding block in M.

function s = gs_Struct(M, groupSize)
g = getGroupIdx(size(M), groupSize);
s = sum(double((M(g)~=0)), 1);
s = reshape(s, [size(M, 2)/groupSize(2), size(M, 1)/groupSize(1)]);
s = s';