%FACMETHOD   Factorize the mesurement matrix by block sparse coding
%   [Pihat, Bhat] = FACMETHOD(W, K) factorize the measurement matrix W, a
%   2F-by-P matrix, into a 2F-by-3K matrix Pihat and 3K-by-P matrix Bhat.t
%
%   [Pihat, Bhat, Score] = FACMETHOD(..., 'Recursive', true) returns a score
%   of the factorization showing how good the dictionary is.
%
%   [Pihat, Bhat] = FACMETHOD(W, K, option) factorized with
%   specified parameters.

function [Pihat, Bhat, Score] = facMethod(W, K, option)
if nargin < 3
    option = facOption();
end

admmOpt = option.admmOption;
ksvdOpt = option.ksvdOption;
Recursive = option.Recursive;
ThresScore = option.ThresScore;
MaxRecr = option.MaxRecr;
ThresGood = option.ThresGood;
PriorConst = option.PriorConst;
Sprior = option.Sprior;

if ~isempty(Sprior)
    W = W/Sprior;
elseif ~isempty(PriorConst)
    Wprior = diag(sqrt(sum(W.^2, 1))/PriorConst);
    W = W/Wprior;
end

%% Main
if ~Recursive
    % GSD
    [~, ~, ~, Xinit, Dinit] = gs_Factorization(W', K, [3, 2], admmOpt);

    % KSVD
    ksvdOpt.Dinit = Dinit;
    ksvdOpt.Xinit = Xinit;
    [Dhat, Xhat, details] = ksvd(W', K, [3, 2], ksvdOpt);
    details.Xinit = Xinit;
    details.Dinit = Dinit;
    if isempty(Sprior) && ~isempty(PriorConst)
        details.Wprior = Wprior;
    end
    Score = details;
    Pihat = Xhat';
    Bhat = Dhat';
    if ~isempty(Sprior)
        Bhat = Bhat * Sprior;
    elseif ~isempty(PriorConst)
        Bhat = Bhat * Wprior;
    end
else

    % GSD
    Dinit = cell(2, 1);
    for iGSD = 1:2
        [~, ~, ~, ~, Dinit{iGSD}] = gs_Factorization(W', K, [3, 2], admmOpt);
    end

    % KSVD
    fprintf('\n\nRecursive Factorization Method: \n')
    fprintf('---------------------------------------------------\n')
    Dhat = cell(2, 1);
    Xhat = cell(2, 1);
    results = cell(2, 1);
    fixAtomSet = cell(2, 1);
    Score = 0;
    for r = 1:MaxRecr
        for iKSVD = 1:2
            fprintf('BEGIN_KSVD%d\n', iKSVD);
            ksvdOpt.Dinit = Dinit{iKSVD};
            [Dhat{iKSVD}, Xhat{iKSVD}, results{iKSVD}] = ksvd(W', K, [3, 2], ksvdOpt);
            fprintf('END_KSVD%d\n\n', iKSVD);
        end

        % Find `good' atoms
        T = computeT(Dhat{1}, Dhat{2});
        Tstruct = gs_Struct(T, [3, 3]);
        for iT = 1:size(Tstruct, 1)
            jT = find(Tstruct(iT, :), 1);
            blockT = T(iT*3-2:iT*3, jT*3-2:jT*3);
            if norm(abs(blockT) - eye(3)) < ThresGood
                if ismember(iT, fixAtomSet{1})
                    assert(ismember(jT, fixAtomSet{2}));
                else
                    fixAtomSet{1} = [fixAtomSet{1}; iT];
                    fixAtomSet{2} = [fixAtomSet{2}; jT];
                end
            end
        end

        % update Initial point Dinit
        if isempty(fixAtomSet{1})
            % need to re-run GSD to get a new initial point
            fprintf('Restart ADMM for new initial points.\n')
            for iGSD = 1:2
                [~, ~, ~, ~, Dinit{iGSD}] = gs_Factorization(W', [3, 2], admmOpt);
            end
            continue;
        else
            % randomly generate not fixed atoms
            % NOTE: here might use some other strategies.
            for iDinit = 1:2
                Dinit{iDinit} = Dhat{iDinit};
                freeAtom = setdiff(1:K, fixAtomSet{iDinit});
                Dinit{iDinit}(:, [3*freeAtom-2, 3*freeAtom-1, 3*freeAtom]) = ...
                    randn(size(Dinit{iDinit}, 1), 3*numel(freeAtom));
            end
        end

        for iViz = 1:2
            fprintf('--Fixed Atoms in D{%d}:%d', iViz, fixAtomSet{iViz}(1));
            fprintf(',%d', fixAtomSet{iViz}(2:end));
            fprintf('.\n');
        end

        Score = numel(fixAtomSet{1})/K;
        if Score >= ThresScore
            fprintf('--SCORE = %.2f >= ThresScore = %.2f\n', Score, ThresScore);
            fprintf('Recursive Factorization Method END.\n')
            Pihat = Xhat{1}';
            Bhat = Dhat{1}';
            if ~isempty(PriorConst)
                Bhat = Bhat * Wprior;
            end
            return;
        else
            fprintf('--SCORE = %.2f < ThresScore = %.2f\n', Score, ThresScore);
            fprintf('---------------------------------------------------\n')
        end
    end
    fprintf('Reach Maximum Iterations.\n')
    fprintf('Recursive Factorization Method END.\n')
    Pihat = Xhat{1}';
    Bhat = Dhat{1}';
    if ~isempty(PriorConst)
        Bhat = Bhat * Wprior;
    end
end
