function viz_S_sharp(S_sharp, varargin)
% Don't use property 'GT'! doesn't work.

[F, p] = size(S_sharp); p = p/3;

S = zeros(3*F, p);
S(1:3:end, :) = S_sharp(:, 1:3:end);
S(2:3:end, :) = S_sharp(:, 2:3:end);
S(3:3:end, :) = S_sharp(:, 3:3:end);

viz_S(S, varargin{:});

%if ~isempty(GT_sharp)
%    GT = zeros(3*F, p);
%    GT(1:3:end, :) = GT_sharp(:, 1:3:end);
%    GT(2:3:end, :) = GT_sharp(:, 2:3:end);
%    GT(3:3:end, :) = GT_sharp(:, 3:3:end);
%else
%    viz_S(S, sec);
%end
