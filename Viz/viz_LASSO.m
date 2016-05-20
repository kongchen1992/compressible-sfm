function viz_LASSO(y, A, Xs, x0, option)
if nargin < 5
    viz_LASSO(y, A, Xs, x0, 'ErrorY - L1_norm');
    viz_LASSO(y, A, Xs, x0, 'ErrorY - L0_norm');
    viz_LASSO(y, A, Xs, x0, 'X - L1_norm');
    viz_LASSO(y, A, Xs, x0, 'X - L0_norm');
    viz_LASSO(y, A, Xs, x0, 'ErrorX - L1_norm');
    viz_LASSO(y, A, Xs, x0, 'ErrorX - L0_norm');
    return;
end

figure;
switch option
    case 'ErrorY - L1_norm'
        l1norm = sum(abs(Xs));
        err = sqrt(sum((repmat(y, 1, size(Xs, 2)) - A*Xs).^2, 1));
        plot(l1norm, err, '*-');
        title('Error vs l1 norm of code');
        xlabel('||x||_1');
        ylabel('||y - Ax||_2');
        
    case 'ErrorY - L0_norm'
        l0norm = sum(double(Xs ~= 0), 1);
        err = sqrt(sum((repmat(y, 1, size(Xs, 2)) - A*Xs).^2, 1));
        plot(l0norm, err, '*-');
        title('Error vs l0 norm of code');
        xlabel('||x||_0');
        ylabel('||y - Ax||_2');
        
    case 'X - L1_norm'
        l1norm = sum(abs(Xs));
        plot(l1norm, Xs', '*-');
        title('the values of code vs l1 norm of code');
        xlabel('||x||_1');
        ylabel('The values of x');
        
    case 'X - L0_norm'
        l0norm = sum(double(Xs ~= 0), 1);
        plot(l0norm, Xs', '*-');
        title('the values of code vs l0 norm of code');
        xlabel('||x||_0');
        ylabel('The values of x');
        
    case 'ErrorX - L1_norm'
        l1norm = sum(abs(Xs));
        err = sqrt(sum((repmat(x0, 1, size(Xs, 2)) - Xs).^2, 1));
        plot(l1norm, err, '*-');
        title('Error vs l1 norm of code');
        xlabel('||x||_1');
        ylabel('||x-x0||_2');
    case 'ErrorX - L0_norm'
        l0norm = sum(double(Xs ~= 0), 1);
        err = sqrt(sum((repmat(x0, 1, size(Xs, 2)) - Xs).^2, 1));
        plot(l0norm, err, '*-');
        title('Error vs l0 norm of code');
        xlabel('||x||_0');
        ylabel('||x - x0||_2');
end