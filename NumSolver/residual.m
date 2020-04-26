function res = residual(x)

% x is passed in as a column vector of initial guesses
res = zeros(2,1);
res(1,1) = x(1)^2 + x(2)^2 - 4;
res(2,1) = x(1)*x(2) + x(2) - 1;