% Create factorization options
%       Name                Value
%       'gamma'             A positive number balancing cost function. 
%                           Default is 1.
%       'maxIter'           A positive number specifying the max number of 
%                           iterations.  The algorithm will stop when 
%                           reaching maxIter, even if it doesn't converge.
%                           Default is 3e5.
%       'isQuiet'           True or False, determining whether the function
%                           prints information on screen.
%                           Default is false;
%       'initRho'           A positive number, specifying the initial value of 
%                           Rho. Check paper for more information of 'Rho'.
%                           Default is 1e-1.
%       'initMu'            A positive number, specifying the initial value of Mu.
%                           Check paper for more information of 'Mu'.  
%                           Default is 1e-2.
%       'maxRho'            A positive number, specifying the max value of Rho.
%                           Default is 1e1;
%       'maxMu'             A positive number, specifying the max value of Mu.
%                           Default is 1e6;
%       'tauRho'            A positive number in [1.01, 1.5], specifying the speed
%                           of Rho's increase.  Check paper for more
%                           information.  
%                           Default is 1.001.
%       'tauMu'             A positive number in [1.01, 1.5], specifying the speed
%                           of Mu's increase.  Check paper for more
%                           information.  
%                           Default is 1.01.
%       'stUpdate'          A positive number specifying the number of
%                           iteration that mu and rho begin to be updated.
%                           The intuition is when mu and rho is small, it's
%                           faster to move to a good initial point. and
%                           then from this point, let's try to converge.
%                           Default 0.
%       'stopCritieria'     A positive number specifying the stop
%                           critieria. If the d(cost) = cost(i) - cost(i-1)
%                           is less than stopCretieria*cost(i), and this
%                           situation holds for 100 iterations.  We say it
%                           coverges.  
%                           Default is 1e-7.
%       'vizIter'           A positive number specifying how many iterations
%                           ADMM print information on screen once.
%                           Default 100.
%       'tauGamma'          A positive number, specifying the speed of gamma's
%                           decreasing.
%                           Default 1.
%       'minGamm'           A positive number, specifying the limit of gamma.
%                           Default 1e-8.
%                           
function opt = admmOption(varargin)
ip = inputParser;

defaultGamma = 1;
defaultMaxIter = 3e5;
defaultIsQuiet = false;
defaultInitRho = 1e-1;
defaultInitMu = 1e2;
defaultMaxRho = 1e1;
defaultMaxMu = 1e6;
defaultTauRho = 1.001;
defaultTauMu = 1.01;
defaultinitA = [];
defaultinitB = [];
defaultStartUpdate = 0;
defaultstopCritieria = 1e-7;
defaultvizIter = 100;
defaultTauGamma = 1;
defaultMinGamma = 1e-8;

addOptional(ip, 'gamma', defaultGamma, @isnumeric);
addOptional(ip, 'maxIter', defaultMaxIter, @isnumeric);
addOptional(ip, 'isQuiet', defaultIsQuiet, @islogical);
addOptional(ip, 'initRho', defaultInitRho, @isnumeric);
addOptional(ip, 'initMu', defaultInitMu, @isnumeric);
addOptional(ip, 'maxRho', defaultMaxRho, @isnumeric);
addOptional(ip, 'maxMu', defaultMaxMu, @isnumeric);
addOptional(ip, 'tauRho', defaultTauRho, @isnumeric);
addOptional(ip, 'tauMu', defaultTauMu, @isnumeric);
addOptional(ip, 'stUpdate', defaultStartUpdate, @isnumeric);
addOptional(ip, 'initB', defaultinitB, @ismatrix);
addOptional(ip, 'initA', defaultinitA, @ismatrix);
addOptional(ip, 'stopCritieria', defaultstopCritieria, @isnumeric);
addOptional(ip, 'vizIter', defaultvizIter, @isnumeric);
addOptional(ip, 'tauGamma', defaultTauGamma, @isnumeric);
addOptional(ip, 'minGamma', defaultMinGamma, @isnumeric);

parse(ip, varargin{:});

opt.gamma = ip.Results.gamma;
opt.maxIter = ip.Results.maxIter;
opt.isQuiet = ip.Results.isQuiet;
opt.initMu = ip.Results.initMu;
opt.initRho = ip.Results.initRho;
opt.maxRho = ip.Results.maxRho;
opt.maxMu = ip.Results.maxMu;
opt.tauRho = ip.Results.tauRho;
opt.tauMu = ip.Results.tauMu;
opt.stopCritieria = ip.Results.stopCritieria;
opt.initA = ip.Results.initA;
opt.initB = ip.Results.initB;
opt.stUpdate = ip.Results.stUpdate;
opt.vizIter = ip.Results.vizIter;
opt.tauGamma = ip.Results.tauGamma;
opt.minGamma = ip.Results.minGamma;
