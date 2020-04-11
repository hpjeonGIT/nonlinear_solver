% sample Newton raphson

% inexact newton
x0   = [0, 0]';
xnew = [0, 0]';
res = 1; i=0;
printf("Using GMRES with Newton\n")
i=0;
eta = 0.1;
norm1 = 1.0; norm2 = 2.0;ds= [0,0]'; res= 1.0;
while i<20 && res > 1.e-4 && norm2 > eta*norm1 %&& i<5 && res > 1.e-5;
 xnew += ds;
 [y, T] = f(xnew(1),xnew(2));
 [J, Jt] = df(xnew(1),xnew(2));
 [ds, flag, reires,iter,resvec] = gmres(-J, y, [], [], maxit=1,[], [], x0=ds);
 %[ds, flag, reires,iter,resvec] = gmres(-Jt, T, [], [], maxit=1,[], [], x0=ds);
 norm2 = norm(J*ds + y)
 norm1 = norm(T) 
 res = norm(ds);
 i = i + 1;
end
fflush(stdout()); 
printf("i=%d\n",i)
xnew

% inexact newton
x0   = [0, 0]';
xnew = [0, 0]';
res = 1; i=0;
printf("Using GMRES with ASPIN\n")
i=0;
eta = 1.e-5;
norm1 = 1.0; norm2 = 2.0;ds= [0,0]'; res= 1.0;
while i<20 && res > 1.e-4 %norm2 > eta*norm1 %&& i<5 && res > 1.e-5;
 xnew += ds;
 [y, T] = f(xnew(1),xnew(2));
 [J, Jt] = df(xnew(1),xnew(2));
 %[ds, flag, reires,iter,resvec] = gmres(-J, y, [], [], maxit=1,[], [], x0=ds);
 [ds, flag, reires,iter,resvec] = gmres(-Jt, T, [], [], maxit=1,[], [], x0=ds);
 norm2 = norm(Jt*ds + T)
 norm1 = norm(T)  
 res = norm(ds);
 i = i + 1;
end
fflush(stdout()); 
printf("i=%d\n",i)
xnew