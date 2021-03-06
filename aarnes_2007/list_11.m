Grid.Nx=60; Grid.hx=20*.3048;   % Dimension in x−direction
Grid.Ny=220; Grid.hy=10*.3048;  % Dimension in y−direction
Grid.Nz=1; Grid.hz=2*.3048;     % Dimension in z−direction
N=Grid.Nx*Grid.Ny*Grid.Nz;      % Number of grid celles
Grid.V=Grid.hx*Grid.hy*Grid.hz; % Volume of each cells
Fluid.vw=3e-4; Fluid.vo=3e-3;   % Viscosities
Fluid.swc=0.2; Fluid.sor=0.2;   % Irreducible saturations
St = 5;                         % Maximum saturation time step
Pt = 100;                       % Pressure time step
ND = 2000;                      % Number of days in simulation

Q=zeros(Grid.Nx,Grid.Ny,1);     % Source term for injection
IR=795*(Grid.Nx*Grid.Ny/ (60*220*85)); % and production. Total
Q(1,1,:)=IR; Q(Grid.Nx,Grid.Ny,:)=-IR; Q=Q(:); % rate scaled to one layer

load Udata; Grid.K=KU(:,1:Grid.Nx,1:Grid.Ny,1); % Permeability in layer 1
Por=pU(1:Grid.Nx,1:Grid.Ny,1); % Preprocessed porosity in layer 1
Grid.por=max(Por(:),1e-3); 
S=Fluid.swc*ones(N,1); % Initial saturation
Pc =[0; 1]; Tt=0;      % For production curves
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