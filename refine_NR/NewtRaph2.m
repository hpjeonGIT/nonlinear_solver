function S=NewtRaph2(Grid,S,Fluid,V,q,T);
N = Grid.Nx*Grid.Ny*Grid.Nz; % number of unknowns
A = GenA(Grid,V,q);          % system matrix
conv=0; IT=0; S00=S;
global prop;
while conv==0;
  dt = T/2^IT;                        % timestep
  dtx = dt./(Grid.V(:)*Grid.por (:)); % timestep / pore volume
  fi = max(q,0).*dtx;                 % injection  
  B=spdiags(dtx,0,N,N)*A;
  I=0;
  while I<2^IT; % loop over sub-timesteps
  S0=S; dsn=1; it=0; I=I+1;
    it_gmres=0;
    while dsn>1e-5 && it<10;           % Newton-Raphson iteration
      [Mw,Mo,dMw,dMo]=RelPerm(S,Fluid); % mobilities and derivatives
      df=dMw./(Mw + Mo)-Mw./(Mw+Mo).^2.*(dMw+dMo);  % df w/ds
      dG=speye(N)-B*spdiags(df,0,N,N); %G’(S)

      fw = Mw./(Mw+Mo);    % fractional flow
      G = S-S0-(B*fw+fi);  % G(s)
      %ds = -dG\G;          % increment ds. Default linear solver
      % in matlab/octave, solve (A*x = b) with 'x = A \ b', rather than 'x = inv (A) * b'.
      % basically it corresponds to mldivide()
      % Therefore, dG*ds + G = 0
      % A => -dG, x => ds, b =>G
      % using gmres with ilu preconditioner
      % inexact newton
      [M1,M2] = ilu(-dG);
      eta = 0.9; %0.5 874sec, 0.9=>.783sec
      [ds, flag, reires,iter,resvec] = gmres(-dG, G, [], rtol=eta, maxit=100,M1,M2);
      it_gmres += sum(iter);                                              
      S = S+ds;            % update S
      dsn = norm(ds);      % norm of increment
      it = it+1;           % number of N-R iterations
    end 
    printf("I= %d NR iters = %d  GMRES iters=%d\n ",I, it, it_gmres);
    fflush(stdout()); 
    prop.nl_it = [prop.nl_it, it];
    prop.lin_it = [prop.lin_it, it_gmres];
    if dsn>1e-3; I=2^IT; S=S00; end
  end
  if dsn<1e-3; 
    conv=1;
  else 
    IT=IT+1; 
  end
  
end