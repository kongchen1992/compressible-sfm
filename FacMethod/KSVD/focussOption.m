%FOCUSSOPTION specifies the parameters in FOCUSS.
%       Name                Value
%       'Wak'               An additonal weight matrix that is dependent of the
%                           a posteriori constraints, the same notation in
%                           paper.
%                           Default identity matrix.
%       'l'                 A postive integer, the same notation in paper.
%                           Default 1.
%       'epsilon'           A small number as the threshold for hard
%                           thresholding on p.
%                           Default 1e-5.
%       'MaxIter'           A positive number denoting how many iterations
%                           FOCUSS run most.
%       'InitX'             A (block) vector as an initial point. If this field
%                           is empty, FOCUSS will use minimum norm algorithm as
%                           initial point.
%                           Default A\b.
%       'Tikhonov'          True or false, indicating if FOCUSS use Tikhonov
%                           Regularization method.
%                           Default false;
%       'lambda'            A number as Tikhonov regularization constant. If it
%                           were empty, Tikhonov Regularization would find
%                           optimal lambda by searching the 'corner' of an
%                           L-curve.
%                           Default 2;
%       'TSVD'              True or false, indicating if FOCUSS use Truncated
%                           Singular Value decomposition.
%                           Default false;
%       't'                 A number as TSVD constant. If it were empty, TSVD
%                           would find optimal t by searching the 'corner' of an
%                           L-curve. 
%                           Default 20.
%       'Trunc'             True or false, indicating if FOCUSS truncates the
%                           the all-zero columns in Aw to make Aw better
%                           conditioned.
%                           Default false.
%       'BlurC'             A positive number as ratio of blurring.
%                           Default 0.5;
%       'HardThres'         True or false, indicating if FOCUSS use hard
%                           thresholding to save computation time.
%                           Default true.
%       'p'                 A positive integer as the 1-st dimension of block in
%                           coefficient matrix. This parameter is designed for
%                           group sparse coding.
%                           Default 1.

function option = focussOption(varargin)
ip = inputParser;

defaultWak = [];
defaultl = 1;
defaultepsilon = 1e-5;
defaultMaxIter = 10;
defaultInitX = [];
defaultTikhonov = false;
defaultlambda = 2;
defaultTSVD = false;
defaultt = 20;
defaultTrunc = false;
defaultBlurC = 0.5;
defaultHardThres = true;
defaultp = 1;

addOptional(ip, 'Wak', defaultWak, @isnumeric);
addOptional(ip, 'l', defaultl, @isnumeric);
addOptional(ip, 'epsilon', defaultepsilon, @isnumeric);
addOptional(ip, 'MaxIter', defaultMaxIter, @isnumeric);
addOptional(ip, 'InitX', defaultInitX, @isnumeric);
addOptional(ip, 'Tikhonov', defaultTikhonov, @islogical);
addOptional(ip, 'lambda', defaultlambda, @isnumeric);
addOptional(ip, 'TSVD', defaultTSVD, @islogical);
addOptional(ip, 't', defaultt, @isnumeric);
addOptional(ip, 'Trunc', defaultTrunc, @islogical);
addOptional(ip, 'BlurC', defaultBlurC, @isnumeric);
addOptional(ip, 'HardThres', defaultHardThres, @islogical);
addOptional(ip, 'p', defaultp, @isnumeric);

parse(ip, varargin{:});
option = ip.Results;
