%KSVDOPTION specifies some optional parameters in K-SVD.
%       Name                Value
%       'Dinit'             A matrix used as the initilization of D.
%                           Default is a randomly generated matrix.
%       'Xinit'             A matrix used as the initilization of X.
%                           Default is minimum norm solution.
%       's'                 A positive number as the number of active elements
%                           in signal.
%                           Default is 2.
%       'sMax'              A positive number as the upper bound of the number
%                           of active elements in signal.
%                           Default is 2.
%       'sMin'              A positive number as the lower bound of the number
%                           of active elements in signal.
%                           Default is 1.
%       'IsOMP'             True or false, indicating if OMP is used for sparse
%                           coding.
%                           Default is false.
%       'IsFOCUSS'          True or false, indicating if FOCUSS is used for
%                           sparse coding.
%                           Default is false.
%       'focussOption'      A structure containing the values of parameters
%                           necessary for FOCUSS. It is generated by
%                           focussOption.m
%       'MaxIter'           A positive number denoting the maximum iterations.
%                           Default is 1e2.
%       'NumCore'           A positive number denoting the numbers of cores to
%                           use in KSVD.  If NumCore > 1, KSVD will exploit
%                           parallel computing.
%                           Default is 1.
%       'StopCri'           A positive number denoting the stop criterion. If
%                           the cost is less than StopCri, KSVD will stop.
%                           Default is 2e-4.
%       'IsViz'             True or False, indicating if KSVD plots figures or
%                           not.
%                           Default is false.
%       'LocalTrap'         A positive number. If the difference of cost between
%                           each iterations is less than LocalTrap, it is
%                           considered that algorithm fall into local trap.
%                           Default is 1e-1.
%       'IsReplace'         True or False, indicating whether KSVD uses strategy
%                           'replace' to avoid local trap or not.  Check the 
%                           paper below for more information of strategy
%                           'replace'.
%                           Default is true.
%       'IsPrune'           True or False, indicating whether KSVD uses the
%                           strategy 'prune' to avoid local trap or not. Check
%                           the paper below for more information of the strategy
%                           'prune'.
%                           Default is false.
%       'IsCut'		        True or False, indicating whether KSVD uses the strategy
%                           'cut' to cut the most used atoms in dictionary.
%                           Default is true.
%       'ThresReplace'      A positive number denoting which atom is not 
%                           being used enough.
%                           Default is 0.2*N*2/K (0.2 times average numbers
%                           that atoms are used).
%       'ThresPrune'        A positive number denoting the threshold of prune.
%                           Default is 0.5.
%       'ThresCut'          A positive number denoting shich atome is being used 
%                           too much.
%                           Default is 2*N*2/K (double average numbers that
%                           atoms are used).
%       'LocalIter'         A positive number denoting how many iterations
%                           need to be done after using strategy getting
%                           rid of local trap before the next.
%       'AtomSet'           A vector indicating which atoms are fixed.
%       'HoldOMP'           True or false, a switch of 2 mode of OMP. The basic
%                           idea is that OMP might choose different atoms in
%                           each iterations. But if the new candidate atoms
%                           cannot model signal as good as the old one, holdOMP
%                           make OMP holds old atom to get rid of increasing
%                           objective function.
%       'lambdaTau'         A positive number in (0,1), indicating the speed of
%                           shrinkage of lambda.]
%       'MinLambda'         A small positive number, as the minimum value of
%                           lambda.
%                           Default is 1e-15.
%       'DMask'             A matrix serving as the mask of dictionary,
%                           indicating which element is supposed to be active.
%                           Default is fullly active.

function option = ksvdOption(varargin)
% parse varargin
ip = inputParser;

defaultDinit = [];
defaultXinit = [];
defaults = 2;
defaultsMin = 1;
defaultsMax = defaults;
defaultIsOMP = false;
defaultIsFOCUSS = false;
defaultfocussOption = focussOption();
defaultMaxIter = 1e2;
defaultNumCore = 1;
defaultStopCri = 2e-4;
defaultIsViz = true;
defaultLocalTrap = 1e-3;
defaultIsReplace = true;
defaultIsPrune = false;
defaultIsCut = true;
defaultThresReplace = [];
defaultThresPrune = 0.5;
defaultThresCut = [];
defaultLocalIter = 10;
defaultAtomSet = [];
defaultHoldOMP = false;
defaultlambdaTau = 0.9;
defaultMinLambda = 1e-15;
defaultDMask = [];

addOptional(ip, 'Dinit', defaultDinit, @isnumeric);
addOptional(ip, 'Xinit', defaultXinit, @isnumeric);
addOptional(ip, 's', defaults, @isnumeric);
addOptional(ip, 'sMin', defaultsMin, @isnumeric);
addOptional(ip, 'sMax', defaultsMax, @isnumeric);
addOptional(ip, 'IsOMP', defaultIsOMP, @islogical);
addOptional(ip, 'IsFOCUSS', defaultIsFOCUSS, @islogical);
addOptional(ip, 'focussOption', defaultfocussOption, @isstruct);
addOptional(ip, 'MaxIter', defaultMaxIter, @isnumeric);
addOptional(ip, 'NumCore', defaultNumCore, @isnumeric);
addOptional(ip, 'StopCri', defaultStopCri, @isnumeric);
addOptional(ip, 'IsViz', defaultIsViz, @islogical);
addOptional(ip, 'LocalTrap', defaultLocalTrap, @isnumeric);
addOptional(ip, 'IsReplace', defaultIsReplace, @islogical);
addOptional(ip, 'IsPrune', defaultIsPrune, @islogical);
addOptional(ip, 'IsCut', defaultIsCut, @islogical);
addOptional(ip, 'ThresReplace', defaultThresReplace, @isnumeric);
addOptional(ip, 'ThresPrune', defaultThresPrune, @isnumeric);
addOptional(ip, 'ThresCut', defaultThresCut, @isnumeric);
addOptional(ip, 'LocalIter', defaultLocalIter, @isnumeric);
addOptional(ip, 'AtomSet', defaultAtomSet, @isnumeric);
addOptional(ip, 'HoldOMP', defaultHoldOMP, @islogical);
addOptional(ip, 'lambdaTau', defaultlambdaTau, @isnumeric);
addOptional(ip, 'MinLambda', defaultMinLambda, @isnumeric);
addOptional(ip, 'DMask', defaultDMask, @isnumeric);

parse(ip, varargin{:});
option = ip.Results;

if ~(option.IsOMP || option.IsFOCUSS)
    option.IsOMP = true;
end