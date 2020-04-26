function [xrt, n_iter] = func_MDnewton(R, J, xi, tol, maxIter, toggle)

corr = -J(xi)\R(xi);
xi_new = xi;
n = 0; % number of iterations
if (toggle == 1)
    fprintf('%9s %12s   %1s\n','Iteration','Max Error','Current root approximation');
end

while (max(abs(corr)) > tol)
    n = n + 1;
    Rxi = R(xi);
    Jxi = J(xi);
    corr = -Jxi\Rxi;
    xi_new = xi + corr;
    if (toggle == 1)
        g = sprintf('%9.3f', xi);
        fprintf('%6i    %12.2e %s\n', n, max(abs(corr)), g);
    end
    xi = xi_new;
    if (n >= maxIter)
        break;
    end
end
% found solution (or exceeded maxIter)
xrt = xi_new;
n_iter = n;
