tic
Grid.Nx=60; Grid.hx=20*.3048;   % Dimension in x−direction
Grid.Ny=220; Grid.hy=10*.3048;  % Dimension in y−direction
Grid.Nz=1; Grid.hz=2*.3048;     % Dimension in z−direction
N=Grid.Nx*Grid.Ny*Grid.Nz;      % Number of grid celles
Grid.V=Grid.hx*Grid.hy*Grid.hz; % Volume of each cells
Fluid.vw=3e-4; Fluid.vo=3e-3;   % Viscosities
Fluid.swc=0.2; Fluid.sor=0.2;   % Irreducible saturations
St = 5;                         % Maximum saturation time step
Pt = 100;                       % Pressure time step
ND = 200;                      % Number of days in simulation

Q=zeros(Grid.Nx,Grid.Ny,1);     % Source term for injection
IR=795*(Grid.Nx*Grid.Ny/ (60*220*85)); % and production. Total
Q(1,1,:)=IR; Q(Grid.Nx,Grid.Ny,:)=-IR; Q=Q(:); % rate scaled to one layer

load Udata; Grid.K=KU(:,1:Grid.Nx,1:Grid.Ny,1); % Permeability in layer 1
Por=pU(1:Grid.Nx,1:Grid.Ny,1); % Preprocessed porosity in layer 1
%
% Grid refine: 2x2 or 4x4 or 8x8
Nscale = 1;
KU_ = zeros(3,Grid.Nx*Nscale,Grid.Ny*Nscale,1);
pU_ = zeros(Grid.Nx*Nscale,Grid.Ny*Nscale,1);
for j=1:Grid.Ny*Nscale
  n = round(j/Nscale); if n < 1; n =1;end
  for i=1:Grid.Nx*Nscale
    m = round(i/Nscale); if m < 1; m=1;end
    KU_(:,i,j,1) = KU(:,m,n,1)*(1.0+rand());
    pU_(i,j,1) = pU(m,n,1)*(1.0+rand());
  end
end
KU = KU_; Grid.K = KU;
pU = pU_; Por = pU;
Grid.Nx *= Nscale;Grid.Ny *= Nscale;  % Grid.Nz is not changed
Grid.hx /= Nscale;Grid.hy /= Nscale; Grid.hz /=Nscale; %  hz is changed
N *= Nscale^2;
Grid.V = Grid.hx*Grid.hy*Grid.hz;
Q=zeros(Grid.Nx,Grid.Ny,1);     % Source term for injection
IR=795*(Grid.Nx*Grid.Ny/ (60*220*85*Nscale^3)); % and production. Total
Q(1,1,:)=IR; Q(Grid.Nx,Grid.Ny,:)=-IR; Q=Q(:); % rate scaled to one layer
printf("%d %d \n ", Grid.Nx, Grid.Ny);
Grid.por=max(Por(:),1e-3); 
S=Fluid.swc*ones(N,1); % Initial saturation
Pc =[0; 1]; Tt=0;      % For production curves
global prop;
prop.nl_it = [];
for tp=1:ND/Pt;
  [P,V]=Pres(Grid,S,Fluid,Q);
  for ts=1:Pt/St;
    S=NewtRaph(Grid,S,Fluid,V,Q,St);
    subplot('position' ,[0.05 .1 .4 .8]);
    pcolor(reshape(S,Grid.Nx,Grid.Ny,Grid.Nz)');
    shading flat; caxis([Fluid.swc 1-Fluid.sor ]); 

    [Mw,Mo]=RelPerm(S(N),Fluid); Mt=Mw+Mo;
    Tt=[Tt,(tp-1)*Pt+ts*St];
    Pc=[Pc,[Mw/Mt; Mo/Mt]];
    subplot('position' ,[0.55 .1 .4 .8]);
    plot(Tt,Pc (1,:), Tt,Pc (2,:));
    axis ([0, ND,-0.05,1.05]);
    legend('Water cut','Oil cut' );
    drawnow;
  end
end
toc