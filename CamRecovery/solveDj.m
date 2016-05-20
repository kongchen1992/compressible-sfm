function Dj = solveDj(Aij, algorithm)%, Qj_gt)
if nargin < 2
    algorithm = 'eig';
end

switch algorithm
    case 'eig'
        Aj = cat(1, Aij{:});
        [eigV, eigD] = eig(Aj'*Aj);
        q = eigV(:, 1);
        Qj  =reshape(q, 3, 3);
        % correct symmetric property
        Qj = 0.5*(Qj+Qj');

        [U, S, ~] = svd(Qj);
        Dj = U*sqrt(S);

    case 'cvx'
        Aj = cat(1, Aij{:});
        cvx_begin sdp
            variable Q(3, 3) symmetric
            minimize norm(Aj*vec(Q));
            Q >= 0
            sum(vec(Q)) == 1;
        cvx_end
        [U, S, ~] = svd(Q);
        Dj = U*sqrt(S);

    case 'RANSAC'
        iter = 1e4; bestfit = 0; Dj = [];
        for i = 1:iter
            idx = randperm(numel(Aij), 5);
            M = solveDj(Aij(idx));
            if isempty(M)
                continue;
            end
            n = numel(find(abs(cat(1, Aij{:})*reshape(M*M', 9, 1))<1e-6));
            if n > bestfit
                bestfit = n;
                Dj = M;
            end
        end
    
    otherwise
        error('Wrong algorithm name.')
end
