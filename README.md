## Sample nonlinear solvers
- MATLAB/Octave sample codes for flow in porous media
- Python code might be provided

## Newton Raphson
- Ref: http://www.mj-oystein.no/index_htm_files/ResSimNotes.pdf
- Solve f(x) = 0
- Fig 31. of the reference
- x1 = x0 - f(x0)/f'(x0)
- f'(x0)*(x1-x0) + f(x0) => f'(x_k) dx_k + f(x_k) = 0 : Eqn(44)
- In matrix form:
A = 
df1/dx1 df1/dx2 ... df1/dxN
df2/dx1 df2/dx2 ... df2/dxN
...
dfN/dx1 dfN/dx2 ... dfN/dxN
x = 
dx1, dx2, ... dxN
b = 
-f1, -f2, ... -fsN
Solve Ax=b
- Then update: x = x + alpha *dx. Alpha would be 1 for non-damping cases

## Inexact Newton
- Ref: https://numhpc.org/hiflow3/fileadmin/tutorials/2.0/tut_inexact_newton_method.pdf
- Ref: http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.126.727&rep=rep1&type=pdf
- REF: https://archive.siam.org/books/textbooks/fr16_book.pdf
- MATLAB example code: https://archive.siam.org/books/fa01/
- NR solves f'*dx = -f, f+f'dx=0. Then x=x+dx
  - |f+f'dx| would be the residual
  - |f +f'dx| / |f| could be the ratio of convergence
- Inexact Newton: solve | f + f' dx | <= eta |f| while 0<eta<1
  - In NR, we solve |f + f'dx| = 0 or eta = 0
  - In inexact newton, we can reduce the loops of linear solver as the residual of |f + f'dx| needs to be smaller than eta |f| only
  - x = x + dx, while dx from | f + f' dx | <= eta |f|

## ASPIN
- http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.125.4868&rep=rep1&type=pdf
- Ref: https://arxiv.org/pdf/1607.04254.pdf
- ASPIN is a modified inexact Newton, using nonlinear preconditioner
- Solving |f + f'dx| <= eta |f| => solve |F + F'dx| <=eta |F|
- F = a function of f and x = G(f(x))
- G would be close to inv(f)
- G(f) must be easily computable
- Ref: https://repository.kaust.edu.sa/bitstream/handle/10754/583375/Final%20thesis.pdf?sequence=2&isAllowed=y
  - How Eq.3.7 is made?
- Solve f_i(x - res_i) for local nonlinear function,i=1..N
- F = \sum res_i
- dx is found through \sum J_local^-1 * J_global *dx = res
- Solve IN as | res - \sum J_local^-1 * J_global *dx| <= eta |res|

- Let's split the entire domain into sum of subdomains
  - Locality is gauranteed as sparse system. Just nonlinear
  - Let's introduce F = \sum T_i
  - In subdomain, f_s_i(u - T_i(u)) = 0
- From Gander, on the origins of linear and non-linear preconditioning
Instead of solving:
f1(x1,x2,..., xN) = 0
f2(x1,x,2,...,xN) = 0
...
fN(x1,x,2,...,xN) = 0
ASPIN solves:
x1 = G1(x2,...,xN)
x2 = G2(x1,...,xN)
...
xN = GN(x1,...,xN-1)
- Advantage:
  - Localize Jacobian inverse, reducing computing
- https://link.springer.com/content/pdf/10.1007/s11242-015-0587-5.pdf
  - skogestad called Cai's method as nonlinear domain decomposition preconditioning
  - first level ASPIN - approximated Jacobian
    - C: Compression from the global to local
    - R: reconstruction
    - for Ax=b, Preconditioner P = \sum Ri * inv(Ci A Ri) * Ci
    - Approximated Jacobian = \sum Ri inv(Ci J(u)Ri) Ci J(u)
- https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/SNES/SNESASPIN.html
  - According to this Petsc manual, preconditioner for linear solver is not applied as the nonlinear solver has the preconditioning.
- http://www.ddm.org/DD13/CaiKY.pdf
- https://pde.sciencesconf.org/data/program/cai_Lyon13.pdf
  - For F1(x1,x2) = 0, F2(x1,x2) = 0
  - Let T1 and T2 be F1(x1-T1,x2)=0, F2(x1,x2-T2)=0
  - T1(x1,x2) + T2(x1,x2) = 0 is solved by ASPIN
  - T1,2 often have better conditioning than F1,2 (why?)
    - See sample below (F1,F2, T1, T2). F1 & F2 are 2nd order file T1 and T1 are approximated as first order
- https://kups.ub.uni-koeln.de/9845/1/CDS_TR-2019-17.pdf
  - For F(u) = 0
  - Ri F(u- PiTi(u)) = 0
  - Fa = sum Pi Ti(u) = 0
  - Newton Raphson for Fa
    - u += - inv(DFa(u)) * Fa(u)
- In ASPIN, if the nonlinear function becomes linear, then the nonlinear preconditioner corresponds to block Jacobi preconditioner

## Sample example
- F1(x1,x2) = x1 + 0.1 * x1*x2 + 0.2*x2*x2 - 0.17975 = 0
- F2(x1,x2) = x1*x1 + 1.5*x1*x2 - 0.7 *x2 + 0.41125 = 0
- Let F1(x1-T1,x2) = (x1-T1) + 0.1 * (x1-T1)*x2 + 0.2*x2*x2 - 0.17975 = 0
  - T1 = (x1 + 0.1*x1*x2 + 0.2*x2*x2 -0.17975)/(1 + 0.1*x2)
- Let F2(x1,x2-T2) = x1*x1 + 1.5*x1*(x2-T2) - 0.7 * (x2-T2) + 0.41125 =0
  - T2 = (x1*x1 + 1.5*x1*x2 - 0.7*x2 + 0.41125)/(1.5*x1 -0.7)
- Ref: https://icerm.brown.edu/materials/Slides/tw-15-5

## Sample2 
- Ref: Nonlinear_Schwarz_Preconditioning_]_David_Keyes,_King_Abdullah_University_of_Science_&_Technology.pdf
- ASPIN reduces NR steps as a half
- Inexact NT is 2x faster than exact NR
  - The key of inexact NT is two fold loops 
  - First loop is for the convergence. Residual check might be done
  - Second loop is for the the condition of inequality |f  + df dx | <= eta |f|
  - Iterating until the inequality, find the dx. When dx is found, check if the convergence is met. If not, iterate the first loop
  - The key of inexact NT is that the total cost of the dual loop might be cheaper than the exact NT.

## In reservoir simulator
T= t0, t1, t2, ... t_f
At t_i, 
1) Configure nonlinear system
2) differentiate
3) loop over linear solver
4) when converged, exit linear loop
5) update nonlinear system
6) If not converged, time step cut


## HO expansion
- y = f + f' dx + f''dx*dx/2 = 0
- dx = -f'/f'' +- \sqrt( f'*f' - 2*f''*f)/f''
