% sample Newton raphson
x0   = [0, 0]';
xnew = [0, 0]';
res = 1;
i = 0;
while res > 1.e-5
  [y, T] = f(xnew(1),xnew(2));
  [J, Jt] = df(xnew(1),xnew(2));
  dx= - J\y;
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
printf("Using GMRES with Netwon\n")
while res > 1.e-5
  [y, T] = f(xnew(1),xnew(2));
  [J, Jt] = df(xnew(1),xnew(2));
  [ds, flag, reires,iter,resvec] = gmres(-J, y, [], rtol=1.e-6, maxit=100);
  xnew = xnew+ ds;
  res = norm(ds);
  i = i + 1;
end
printf("i=%d\n",i);
xnew

% using gmres  
x0   = [0, 0]';
xnew = [0, 0]';
res = 1; i=0;
printf("Using GMRES with ASPIN\n")
while res > 1.e-5
  [y, T] = f(xnew(1),xnew(2));
  [J, Jt] = df(xnew(1),xnew(2));
  [ds, flag, reires,iter,resvec] = gmres(-Jt, T, [], rtol=1.e-6, maxit=100);
  xnew = xnew+ ds;
  res = norm(ds);
  i = i + 1;
end
printf("i=%d\n",i);
xnew
