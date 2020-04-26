function [y] = func_surr(x,beta)

B0 = beta(1);
B1 = beta(2);
B2 = beta(3);
B3 = beta(4);
B4 = beta(5);

y = B0 + B1*x + B2*x.^2 + B3*x.^3 + B4*x.^4;