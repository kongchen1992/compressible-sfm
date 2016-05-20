function r = evalD(D_ij, delta)
delDDdel = delta*(D_ij)*D_ij'*delta';

diagdDDd = diag(delDDdel);
diagPdDDd = diag(delDDdel, 1);

r = 1 - diagdDDd(2:2:end)./diagdDDd(1:2:end);
r = [r; 2*diagPdDDd(1:2:end)./diagdDDd(1:2:end)];