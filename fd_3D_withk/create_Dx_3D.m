function [Dx,Dy,Dz] = create_Dx_3D(npts,delx)

D = create_D_matrix(npts,3,delx,1);

Dz = kron(D,kron(speye(npts),speye(npts)));
Dx = kron(speye(npts),kron(D,speye(npts)));
Dy = kron(speye(npts),kron(speye(npts),D));
