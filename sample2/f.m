function [y, T] = f(x1,x2);
  f1 = (x1 - x2^3 + 1)^3 - x2^3;
  f2 = 2*x1 + 3*x2 - 5;
  y = [f1,f2]'; % regular function
  T1 = x1 - x2^3 + 1 - x2;
  T2 = 2.*x1/3. + x2 - 5./3.;
  T = [T1, T2]'; % deduced function from ASPIN
end
%https://icerm.brown.edu/materials/Slides/tw-15-5/Nonlinear_Schwarz_Preconditioning_]_David_Keyes,_King_Abdullah_University_of_Science_&_Technology.pdf
%slide 30