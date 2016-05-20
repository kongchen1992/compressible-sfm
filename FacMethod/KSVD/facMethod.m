%FACMETHOD create facMethod options.
%       Name                Value
%       ''
%       'myFun_KSVD'        A function handle for sparse coding.
%                           Default is @omp.
%       'MaxIter_KSVD'      A positive number denoting the maximum iterations in
%                           KSVD.
%                           Default is 1e3.
%       'NumCore_KSVD'      A positive number denoting the numbers of cores to
%                           use in KSVD.  If NumCore > 1, KSVD will exploit
%                           parallel computing. Usually no need to touch.
%                           Default is 1.
%       'StopCri_KSVD'      A positive number denoting the stop criterion in
%                           KSVD. If the cost is less than StopCri, KSVD will stop.
%                           Default is 2e-4.
%       'IsViz_KSVD'        True or False, indicating if KSVD plots figures or
%                           not.
%                           Default is true.
%       'LocalTrap_KSVD'    A positive number. If the difference of cost between
%                           each iterations in KSVD is less than LocalTrap, it is
%                           considered that KSVD fall into local trap.
%                           Default is 1e-3.
%       'IsReplace_KSVD'    True or False, indicating whether KSVD uses strategy
%                           'replace' to avoid local trap or not.  
%                           Default is true.
%       'IsPrune_KSVD'      True or False, indicating whether KSVD uses the
%                           strategy 'prune' to avoid local trap or not. 
%                           Default is true.
%       'ThresReplace_KSVD' A positive number denoting which atom is not 
%                           being used enough.
%                           Default is 0.2*N*2/K (0.2 times average numbers
%                           that atoms are used).
%       'ThresPrune_KSVD'   A positive number denoting the threshold of prune.
%                           Default is 0.6.
%       'Recursive'         True or False, indicating whether do ADMM KSVD
%                           again and again to get results.
%       'ThresScore'        A positive number between [0, 1], above which
%                           the algorithm stops.
%       'MaxRecr'           A positive integer indicating how many times
%                           KSVD run recursively for finding better
%                           dictionary.
%       'ThresGood'         A positive number give a threshold for deciding
%                           if the block in T is satisfied theorem 1 in
%                           blockKSVD.pdf document.
%       'HoldOMP'           True or false, a switch of 2 mode of OMP. The basic
%                           idea is that OMP might choose different atoms in
%                           each iterations. But if the new candidate atoms
%                           cannot model signal as good as the old one, holdOMP
%                           make OMP holds old atom to get rid of increasing
%                           objective function.
