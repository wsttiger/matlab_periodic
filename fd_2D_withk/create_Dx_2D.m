function [Dx,Dy] = create_Dx_2D(npts,delx)

D = create_D_matrix(npts,7,delx,1);

Dx = kron(D,speye(npts));
Dy = kron(speye(npts),D);

