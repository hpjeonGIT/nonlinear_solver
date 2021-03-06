function [J, Jt] = df(x1,x2);
J = zeros(2);
J(1,1) = 3*(x1-x2^3+1)^2;
J(1,2) = 3*(x1-x2^3+1)^2*(-3.*x2^2)  - 3*x2^2;
J(2,1) = 2.;
J(2,2) = 3.;
Jt = zeros(2);
Jt(1,1) = 1.;
Jt(1,2) = -3*x2^2 - 1.;
Jt(2,1) = 2./3.;
Jt(2,2) = 1.;
end
