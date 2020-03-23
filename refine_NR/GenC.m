function C=GenC(Grid)
Nx=Grid.Nx; Ny=Grid.Ny; Nz=Grid.Nz;
C=sparse(0,0);
Nxy=Nx*Ny; N=Nxy*Nz;
vx=ones(Nx,1); vy=ones(Nxy,1); vz=ones(N,1);
% Empty sparse matrix
% Number of grid-points
% Diagonals
for i=1:Ny*Nz
Cx=spdiags([vx,-vx],[-1,0]-(i-1)*Nx,N,Nx-1);
C=[C,Cx];
end % vx-block of C
% create bidiagonal block
% append to C
for i=1:Nz
Cy=spdiags([vy,-vy],[-Nx,0]-(i-1)*Nxy,N,Nxy-Nx);
C=[C,Cy];
end % vy-block of C
% create bidiagonal block
% append to C
C = [C, spdiags([vz,-vz],[-Nxy,0],N,N-Nxy)]; % vz-block of C