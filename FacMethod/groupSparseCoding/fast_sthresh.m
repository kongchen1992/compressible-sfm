% Set the fast soft thresholding function
function r = fast_sthresh(x,th)
normx = sqrt(sum(x.^2, 1));
r = repmat(max(normx - th,0)./normx, size(x, 1), 1).*x;