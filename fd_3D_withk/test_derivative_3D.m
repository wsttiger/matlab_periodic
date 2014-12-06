clear all;

npts = 80;
L = 5.0;
x = linspace(-L/2,L/2,npts);
delx = x(2)-x(1);
npts = npts-1;
x = x(1:npts);

[X,Y,Z] = meshgrid(x,x,x);
V = cos(2.0*pi*X/L).*cos(2.0*pi*Y/L).*cos(2.0*pi*Z/L);

DxV = -(2.0*pi/L)*sin(2.0*pi*X/L).*cos(2.0*pi*Y/L).*cos(2.0*pi*Z/L);
DyV = -(2.0*pi/L)*cos(2.0*pi*X/L).*sin(2.0*pi*Y/L).*cos(2.0*pi*Z/L);
DzV = -(2.0*pi/L)*cos(2.0*pi*X/L).*cos(2.0*pi*Y/L).*sin(2.0*pi*Z/L);

D = create_D_matrix(npts,7,delx,1);

% meshgrid is a funny function
% y has stride 1
% x has stride npts
% z has stride npts*npts
Dz = kron(D,kron(speye(npts),speye(npts)));
Dx = kron(speye(npts),kron(D,speye(npts)));
Dy = kron(speye(npts),kron(speye(npts),D));

DxV2 = reshape(Dx*V(:), [npts npts npts]);
DyV2 = reshape(Dy*V(:), [npts npts npts]);
DzV2 = reshape(Dz*V(:), [npts npts npts]);


