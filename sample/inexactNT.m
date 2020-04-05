% sample Newton raphson
% x1 = 0.15
% x2 = 0.35
% F1 = x1    + 0.1*x1*x2 + 0.2*x2**2 - 0.17975
% F2 = x1**2 + 1.5*x1*x2 - 0.7*x2    + 0.41125

% inexact newton
x0   = [0, 0]';
xnew = [0, 0]';
res = 1; i=0;
printf("Using GMRES\n")
i=0;
eta = 0.9;
norm1 = 1.0; norm2 = 2.0;ds= [0,0]'; res= 1.0;
while i<5 && res > 1.e-5;
 xnew += ds
 dG = df(xnew(1),xnew(2));
 G = f(xnew(1),xnew(2))';
 [ds, flag, reires,iter,resvec] = gmres(dG, -G, [], [], maxit=1,[], [], x0=ds);
 ds
 norm2 = norm(dG*ds + G)
 norm1 = norm(G)  
 fflush(stdout()); 
 res = norm(ds)
 i = i + 1;
end
printf("i=%d\n",i)
xnew