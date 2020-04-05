% sample Newton raphson
% x1 = 0.15
% x2 = 0.35
% F1 = x1    + 0.1*x1*x2 + 0.2*x2**2 - 0.17975
% F2 = x1**2 + 1.5*x1*x2 - 0.7*x2    + 0.41125
x0   = [0, 0]';
xnew = [0, 0]';
res = 1;
i = 0;
while res > 1.e-5
  dx= - df(xnew(1),xnew(2))\f(xnew(1),xnew(2))';
  xnew = xnew + dx;
  res = norm(dx);
  i = i + 1;
end
printf("i=%d\n",i);
xnew  

% using gmres  
x0   = [0, 0]';
xnew = [0, 0]';
res = 1; i=0;
printf("Using GMRES\n")
while res > 1.e-5
 dG = df(xnew(1),xnew(2));
 G = f(xnew(1),xnew(2))';
 [ds, flag, reires,iter,resvec] = gmres(-dG, G, [], rtol=1.e-6, maxit=100);
 xnew = xnew+ ds;
 res = norm(ds);
 i = i + 1;
end
printf("i=%d\n",i);
xnew
