%KSVD learns an analysis dictionary using KSVD method.
%   [D, X] = KSVD(Y, K, sizeG) learns a dictionary D and block-sparse
%   codes X from data Y, where K determines the number of atoms in D, sizeG
%   denotes the size of block in X.
%
%   [D, X, cost] = KSVD(...) returns the cost values along the iterations.
%
%   [D, X, cost, Xpath, Dpath] = KSVD(...) returns the Xs and Ds in each
%   iterations, called solution path.
%
%   KSVDOPTION is used to specify some parameters in KSVD.
%   Example:
%       opt = ksvdOption('MaxIter', 100);
%       [D, X, cost] = ksvd(..., opt);
%
%
%   Algorithm from paper 'K-SVD: An algorithm for Designing Overcomplete 
%   Dictionaries for Sparse Representation'.  Chen generalizes it to group 
%   sparse situation.
%
function [bestD, bestX, result] = ksvd(Y, K, sizeG, option)
if nargin < 4
    option = ksvdOption();
end

M = size(Y, 1); p = sizeG(1); q = sizeG(2);
N = size(Y, 2)/q;

if isempty(option.Dinit)
    if isempty(option.DMask)
        D = randn(M, p*K);
    else
        D = randn(M, p*K).*option.DMask;
    end
    D = D/diag(sqrt(sum(D.^2, 1)));
else
    D = option.Dinit;
    D = D/diag(sqrt(sum(D.^2, 1)));
end
X = option.Xinit;
s = option.s;
sMin = option.sMin;
sMax = option.sMax;
IsOMP = option.IsOMP;
IsFOCUSS = option.IsFOCUSS;
focussOpt = option.focussOption;
MaxIter = option.MaxIter;
NumCore = option.NumCore;
StopCri = option.StopCri;
IsViz = option.IsViz;
LocalTrap = option.LocalTrap;
IsReplace = option.IsReplace;
IsPrune = option.IsPrune;
IsCut = option.IsCut;
if isempty(option.ThresReplace)
    ThresReplace = 0.2*N*2/K;
else
    ThresReplace = option.ThresReplace;
end
ThresPrune = option.ThresPrune;
if isempty(option.ThresCut)
    ThresCut = 2*N*2/K;
else
    ThresCut = option.ThresCut;
end
LocalIter = option.LocalIter;
AtomSet = option.AtomSet;
HoldOMP = option.HoldOMP;
lambdaTau = option.lambdaTau;
MinLambda = option.MinLambda;
DMask = option.DMask;

readyReplace = true;
readyPrune = true;
readyCut = true;

IndexLocal = 1;

if NumCore > 1
    poolobj = parpool('local', NumCore);
    XnewCell = cell(N, 1);
    Ycell = cell(N, 1);
    focussOptCell = cell(N, 1);
    for j = 1:N
        Ycell{j} = Y(:, sizeG(2)*(j-1)+1:sizeG(2)*j);
    end
end

cost = zeros(MaxIter, 1);
Xpath = cell(MaxIter, 1);
Dpath = cell(MaxIter, 1);

bestCost = inf;
bestD = [];
bestX = [];

%temperary variables, will be removed
Lambda = ones(MaxIter, N);
focussIters = zeros(MaxIter, N);

if IsViz
    figure; hold on;
    axis([1, inf, 0, inf])
    xlabel('iterations')
    ylabel('cost')
end

for i = 1:MaxIter


    % Sparse Coding Stage
    if IsOMP
        if isempty(X)
            X = zeros(p*K, q*N);
        end
        if p ~= 1 || q ~= 1
            Uj = cell(K, 1);
            for j = 1:K
                [Uj{j}, ~, ~] = svd(D(:, p*(j-1)+1:p*j), 'econ');
            end
            U = cat(2, Uj{:});
        else
            U = [];
        end
        if NumCore > 1
            % prepare data for parallel
            Ycell = cell(N, 1);
            XnewCell = cell(N, 1);
            for j = 1:N
                Ycell{j} = Y(:, sizeG(2)*(j-1)+1:sizeG(2)*j);
            end
            % OMP
            parfor j = 1:N
                XnewCell{j} = omp(D, Ycell{j}, s, sizeG, U);
            end
            Xnew = cat(2, XnewCell{:});
        else
            for j = 1:N
                candidate = omp(D, Y(:, sizeG(2)*(j-1)+1:sizeG(2)*j), s, sizeG, U);
                if HoldOMP
                    errNew = norm(Y(:, sizeG(2)*(j-1)+1:sizeG(2)*j) - D*candidate, 'fro');
                    errOld = norm(Y(:, sizeG(2)*(j-1)+1:sizeG(2)*j) - ...
                        D*X(:, sizeG(2)*(j-1)+1:sizeG(2)*j), 'fro');
                    if errNew < errOld
                        Xnew(:, sizeG(2)*(j-1)+1:sizeG(2)*j) = candidate;
                    else
                        blockStruc = gs_Struct(X(:, sizeG(2)*(j-1)+1:sizeG(2)*j), [3,2]);
                        oldAtoms = find(blockStruc~=0);
                        candidate = zeros(size(candidate));
                        indOld = zeros(size(D, 2), 1);
                        for iOldAtom = 1:s
                            indOld((oldAtoms(iOldAtom)-1)*sizeG(1)+1:...
                                oldAtoms(iOldAtom) *sizeG(1)) = 1;
                        end
                        candidate(indOld==1, :) = D(:, indOld==1)\Y(:,sizeG(2)*(j-1)+1:sizeG(2)*j);
                        Xnew(:, sizeG(2)*(j-1)+1:sizeG(2)*j) = candidate;
                    end
                else
                    Xnew(:, sizeG(2)*(j-1)+1:sizeG(2)*j) = candidate;
                end
            end
        end

    elseif IsFOCUSS
        % Compute Wak beforehand.
        Dreshape = reshape(D, size(D, 1)*p, size(D, 2)/p);
        wak = max(abs(Dreshape), [], 1);
        wak = ones(size(wak)) ./ wak;
        Wak = diag(kron(wak, ones(1, p)));
        focussOpt.Wak = Wak;
        if focussOpt.lambda >= MinLambda
            focussOpt.lambda = focussOpt.lambda*lambdaTau;
        end
        if NumCore > 1
            if ~isempty(X)
                for j = 1:N
                    focussOptCell{j} = focussOpt;
                    focussOptCell{j}.InitX = X(:, sizeG(2)*(j-1)+1:sizeG(2)*j);
                end
            else
                for j = 1:N
                    focussOptCell{j} = focussOpt;
                end
            end
            parfor j = 1:N
                [xFOCUSS, ~] = focuss(D, Ycell{j}, focussOptCell{j});
                if sizeG(1) == 1 && sizeG(2) == 1
                    gsStruct = (xFOCUSS{end} ~= 0 );
                else
                    gsStruct = gs_Struct(xFOCUSS{end}, sizeG);
                end
                nGS = numel(find(gsStruct~=0));
                while nGS > sMax || nGS < sMin
                    if nGS > sMax
                        focussOptCell{j}.lambda = focussOptCell{j}.lambda*1.5;
                    elseif nGS < sMin
                        if focussOptCell{j}.Tikhonov
                            if focussOptCell{j}.lambda < 1e-12
                                focussOptCell{j}.Tikhonov = false;
                            else
                                focussOptCell{j}.lambda = focussOptCell{j}.lambda*0.5;
                            end
                        else
                            break;
                        end
                    end
                    [xFOCUSS, ~] = focuss(D, Ycell{j}, focussOptCell{j});
                    if sizeG(1) == 1 && sizeG(2) == 1
                        gsStruct = (xFOCUSS{end} ~= 0 );
                    else
                        gsStruct = gs_Struct(xFOCUSS{end}, sizeG);
                    end
                    nGS = numel(find(gsStruct~=0));
                end
                XnewCell{j} = xFOCUSS{end};
            end
            Xnew = cat(2, XnewCell{:});
            for j = 1:N
                Lambda(i,j) = focussOptCell{j}.lambda;
            end

        else
            for j = 1:N
                if ~isempty(X)
                    focussOpt.InitX = X(:, sizeG(2)*(j-1)+1:sizeG(2)*j);
                end
                [xFOCUSS, ~] = focuss(D, Y(:, sizeG(2)*(j-1)+1:sizeG(2)*j), focussOpt);
                if sizeG(1) == 1 && sizeG(2) == 1
                    gsStruct = (xFOCUSS{end} ~= 0 );
                else
                    gsStruct = gs_Struct(xFOCUSS{end}, sizeG);
                end
                focussOptInner = focussOpt;
                nGS = numel(find(gsStruct~=0));
                focussIters(i, j) = focussIters(i, j) + 1;
                while nGS > sMax || nGS < sMin
                    if nGS > sMax
                        focussOptInner.lambda = focussOptInner.lambda*1.5;
                    elseif nGS < sMin
                        if focussOptInner.Tikhonov
                            if focussOptInner.lambda < 1e-12
                                focussOptInner.Tikhonov = false;
                            else
                                focussOptInner.lambda = focussOptInner.lambda*0.5;
                            end
                        else
                            break;
                        end
                    end
                    [xFOCUSS, ~] = focuss(D, Y(:, sizeG(2)*(j-1)+1:sizeG(2)*j), focussOptInner);
                    focussIters(i, j) = focussIters(i, j) + 1;
                    if sizeG(1) == 1 && sizeG(2) == 1
                        gsStruct = (xFOCUSS{end} ~= 0 );
                    else
                        gsStruct = gs_Struct(xFOCUSS{end}, sizeG);
                    end
                    nGS = numel(find(gsStruct~=0));
                end
                Xnew(:, sizeG(2)*(j-1)+1:sizeG(2)*j) = xFOCUSS{end};
                Lambda(i, j) = focussOptInner.lambda;
            end
        end
    end

    if i == 120
        a = 1;
    end
    % Codebook Update Stage
    for j = 1:K
        if ismember(j, AtomSet)
            continue;
        end
        [Omega, pat] = extractSP(Xnew(p*(j-1)+1:p*j, :), sizeG);
        if isempty(Omega)
            continue;
        end
        Ek = Y - D(:, [1:p*(j-1), p*j+1:end])*Xnew([1:p*(j-1), p*j+1:end], :);
        EkR = Ek*Omega;
        if isempty(DMask)
            [U, S, V] = svd(EkR);
            D(:, p*(j-1)+1:p*j) = U(:, 1:p);
            Xnew(p*(j-1)+1:p*j, pat~=0) = S(1:p, :)*V';
        else
            [U, S, V] = svd(EkR(DMask(:, p*j)~=0, :));
            D(DMask(:, p*j)~=0, p*(j-1)+1:p*j) = U(:, 1:p);
            Xnew(p*(j-1)+1:p*j, pat~=0) = S(1:p, :)*V';
        end
    end

    cost(i) = norm(Y - D*Xnew, 'fro');

    if IsViz
        plot(i, cost(i), '.');
        drawnow;
    else
        fprintf('%d\t%.5f\n', i, cost(i));
    end

    if cost(i) < bestCost
        bestD = D;
        bestX = Xnew;
        bestCost = cost(i);
        Xpath{i} = bestX;
        Dpath{i} = bestD;
    end


    if cost(i) < StopCri
        fprintf('Problem solved.\n');
        if NumCore > 1
            delete(poolobj);
        end
        result.cost = cost(1:i);
        result.Xpath = Xpath(1:i);
        result.Dpath = Dpath(1:i);
        result.Lambda = Lambda(1:i, :);
        result.focussIters = focussIters(1:i, :);
        return;
    elseif i >1 && ( (sizeG(1)==1 && sizeG(2)==1) || (sizeG(1)==3 && sizeG(2)==2))

        if abs( cost(i)-cost(i-1)) < LocalTrap
            if IsReplace || IsPrune || IsCut
                % use technique to get rid of local trap
                Nums_usedAtom = sum(Xnew~=0, 2);
                if sizeG(1) == 3
                    Nums_usedAtom = reshape(Nums_usedAtom, 3, numel(Nums_usedAtom)/3);
                    Nums_usedAtom = max(Nums_usedAtom, [], 1);
                end
                Errors = sum((Y - D*Xnew).^2);
                if sizeG(1) == 3
                    Errors = reshape(Errors, 2, numel(Errors)/2);
                    Errors = sum(Errors, 1);
                end
                IDs_invoidAtom = [];
                [~, IDs_badY] = sort(Errors, 'descend');
            end

            if IsReplace && readyReplace && i - IndexLocal > LocalIter
                % replace the atom in Dictionary that are not used enough
                %   with the signal that are not represented well.
                IDs_invoidAtom = find(Nums_usedAtom < ThresReplace);
                % make sure that atoms in AtomSet are not removed.
                IDs_invoidAtom = setdiff(IDs_invoidAtom, AtomSet);
                if ~isempty(IDs_invoidAtom)
                    fprintf('[Replace]\t');
                    fprintf('Replaced atoms: ');
                    fprintf('%d ', IDs_invoidAtom);
                    fprintf('   at  iteration %d\n', i);
                    if p == 3 && q == 2
                        %badY = zeros(size(Y, 1), 3*numel(IDs_invoidAtom));
                        D(:, IDs_invoidAtom*3-2) = ...
                            Y(:, IDs_badY(1:numel(IDs_invoidAtom))*2-1);
                        D(:, IDs_invoidAtom*3-1) = ...
                            Y(:, IDs_badY(1:numel(IDs_invoidAtom))*2);
                        D(:, IDs_invoidAtom*3) = ...
                            randn(size(Y, 1), numel(IDs_invoidAtom));
                        if ~isempty(DMask)
                            D = D.*DMask;
                        end
                        %normalize 
                        D(:, IDs_invoidAtom*3-2) = D(:, IDs_invoidAtom*3-2)/...
                            diag(sqrt(sum(D(:, IDs_invoidAtom*3-2).^2, 1)));
                        D(:, IDs_invoidAtom*3-1) = D(:, IDs_invoidAtom*3-1)/...
                            diag(sqrt(sum(D(:, IDs_invoidAtom*3-1).^2, 1)));
                        D(:, IDs_invoidAtom*3) = D(:, IDs_invoidAtom*3)/...
                            diag(sqrt(sum(D(:, IDs_invoidAtom*3).^2, 1)));
                    elseif p ==1 && q == 1
                        D(:, IDs_invoidAtom) = ...
                            Y(:, IDs_badY(1:numel(IDs_invoidAtom)));
                        % normalize
                        D(:, IDs_invoidAtom) = D(:, IDs_invoidAtom)/ ...
                            diag(sqrt(sum(D(:, IDs_invoidAtom).^2, 1)));
                    end
                end
                if IsPrune
                    readyReplace = false;
                    readyPrune = true;
                    readyCut = false;
                elseif IsCut
                    readyReplace = false;
                    readyPrune = false;
                    readyCut = true;
                else
                    readyReplace = true;
                    readyPrune = false;
                    readyCut = false;
                end
                IndexLocal = i;
            end

            if IsPrune && readyPrune && i - IndexLocal > LocalIter
                % prune the dictionary from having too-close atoms
                if sizeG(1) == 1 && sizeG(2) == 1
                    CorrAtoms = triu(D'*D, 1);
                    [X_tcAtom, ~] = find(abs(CorrAtoms) > ThresPrune);
                    X_tcAtom = unique(X_tcAtom);
                    % Make sure atoms in AtomSet are not pruned.
                    X_tcAtom = setdiff(X_tcAtom, AtomSet);
                    % should be ascend sorted automatically. In case not, I use sort here.
                    X_tcAtom = sort(X_tcAtom, 1, 'ascend'); 
                    if ~isempty(X_tcAtom)
                        fprintf('[Prune]\t\t')
                        fprintf('Pruned atoms:')
                        fprintf('%d ', X_tcAtom);
                        fprintf('   at  iteration %d\n', i);
                        fprintf('\n');
                        D(:, X_tcAtom) = Y(:, IDs_badY(numel(IDs_invoidAtom)+1:...
                            numel(IDs_invoidAtom)+numel(X_tcAtom)));
                        % Normalize
                        D(:, X_tcAtom) = D(:, X_tcAtom)./ ...
                            repmat(sqrt(sum(D(:,X_tcAtom).^2, 1)),size(D,1),1);
                    end
                elseif sizeG(1) == 3 && sizeG(2) == 2
                    CorrAtoms = D'*D;
                    g = getGroupIdx(size(CorrAtoms), [3, 3]);
                    CorrAtoms = reshape(sum(CorrAtoms(g([1, 5, 9], :)), 1),K,K);
                    CorrAtoms = triu(CorrAtoms, 1)/3;
                    [X_tcAtom, ~] = find(abs(CorrAtoms) > ThresPrune);
                    X_tcAtom = unique(X_tcAtom);
                    % Make sure atoms in AtomSet are not pruned.
                    X_tcAtom = setdiff(X_tcAtom, AtomSet);
                    % should be ascend sorted automatically. In case not, I use sort here.
                    X_tcAtom = sort(X_tcAtom, 1, 'ascend'); 
                    if ~isempty(X_tcAtom)
                        fprintf('[Prune]\t\t')
                        fprintf('Pruned atoms:') 
                        fprintf('%d ', X_tcAtom);
                        fprintf('   at  iteration %d\n', i);
                        fprintf('\n');
                        D(:, X_tcAtom*3-2) = Y(:, IDs_badY(numel(IDs_invoidAtom)...
                            +1:numel(IDs_invoidAtom)+numel(X_tcAtom))*2-1);
                        D(:, X_tcAtom*3-1) = Y(:, IDs_badY(numel(IDs_invoidAtom)...
                            +1:numel(IDs_invoidAtom)+numel(X_tcAtom))*2);
                        D(:, X_tcAtom*3) = randn(size(Y, 1), numel(X_tcAtom));
                        % Normalize
                        D(:, X_tcAtom*3-2) = D(:, X_tcAtom*3-2)./ ...
                            repmat(sqrt(sum(D(:,X_tcAtom*3-2).^2, 1)),size(D,1),1);
                        D(:, X_tcAtom*3-1) = D(:, X_tcAtom*3-1)./ ...
                            repmat(sqrt(sum(D(:,X_tcAtom*3-1).^2, 1)),size(D,1),1);
                        D(:, X_tcAtom*3) = D(:, X_tcAtom*3)./ ...
                            repmat(sqrt(sum(D(:,X_tcAtom*3).^2, 1)),size(D,1),1);
                    end
                end
                if IsCut
                    readyReplace = false;
                    readyPrune = false;
                    readyCut = true;
                elseif IsReplace
                    readyReplace = true;
                    readyPrune = false;
                    readyCut = false;
                else
                    readyReplace = false;
                    readyPrune = true;
                    readyCut = false;
                end
                IndexLocal = i;
            end

            if IsCut && readyCut && i - IndexLocal > LocalIter
                % Cut the atom being used too much
                IDs_invoidAtom = find(Nums_usedAtom > ThresCut);
                % make sure that atoms in AtomSet are not removed.
                IDs_invoidAtom = setdiff(IDs_invoidAtom, AtomSet);
                if ~isempty(IDs_invoidAtom)
                    fprintf('[Cut]\t');
                    fprintf('Cut atoms: ');
                    fprintf('%d ', IDs_invoidAtom);
                    fprintf('   at  iteration %d\n', i);
                    if p == 3 && q == 2
                        D(:, IDs_invoidAtom*3-2) = ...
                            Y(:, IDs_badY(1:numel(IDs_invoidAtom))*2-1);
                        D(:, IDs_invoidAtom*3-1) = ...
                            Y(:, IDs_badY(1:numel(IDs_invoidAtom))*2);
                        D(:, IDs_invoidAtom*3) = ...
                            randn(size(Y, 1), numel(IDs_invoidAtom));
                        %normalize 
                        D(:, X_tcAtom*3-2) = D(:, X_tcAtom*3-2)./ ...
                            repmat(sqrt(sum(D(:,X_tcAtom*3-2).^2, 1)),size(D,1),1);
                        D(:, X_tcAtom*3-1) = D(:, X_tcAtom*3-1)./ ...
                            repmat(sqrt(sum(D(:,X_tcAtom*3-1).^2, 1)),size(D,1),1);
                        D(:, X_tcAtom*3) = D(:, X_tcAtom*3)./ ...
                            repmat(sqrt(sum(D(:,X_tcAtom*3).^2, 1)),size(D,1),1);
                    elseif p ==1 && q == 1
                        D(:, IDs_invoidAtom) = ...
                            Y(:, IDs_badY(1:numel(IDs_invoidAtom)));
                        % normalize
                        D(:, IDs_invoidAtom) = D(:, IDs_invoidAtom)/ ...
                            norm(D(:, IDs_invoidAtom));
                    end
                end
                if IsReplace
                    readyReplace = true;
                    readyPrune = false;
                    readyCut = false;
                elseif IsPrune
                    readyReplace = false;
                    readyPrune = true;
                    readyCut = false;
                else
                    readyReplace = false;
                    readyPrune = false;
                    readyCut = true;
                end
                IndexLocal = i;
            end
        end
    end
    X = Xnew;
end


fprintf('Reach Maximum Iterations.\n');
if NumCore > 1
    delete(poolobj);
end

result.cost = cost;
result.Xpath = Xpath;
result.Dpath = Dpath;
result.Lambda = Lambda;
result.focussIters = focussIters;
