clear all;

% lattice constant
npts = 132;
L = 5.0;
eps = 1.e-6;
x = linspace(-L/2+eps, L/2+eps, npts);
y = linspace(-L/2+eps, L/2+eps, npts);
z = linspace(-L/2+eps, L/2+eps, npts);
delx = x(2) - x(1);
npts = npts-1;
x = x(1:npts); y = y(1:npts); z = z(1:npts);

% setup k-grid
sidesz = (npts-1)/2;
kx = 2*pi*(-sidesz:sidesz)./delx./npts;
ky = 2*pi*(-sidesz:sidesz)./delx./npts;
kz = 2*pi*(-sidesz:sidesz)./delx./npts;
[KX,KY,KZ] = meshgrid(kx,ky,kz);
K2 = KX.^2 + KY.^2 + KZ.^2;

[X,Y,Z] = meshgrid(x,y,z);
rs = sqrt(X.^2 + Y.^2 + Z.^2);

mu = 8.0;
f = (1/4/pi)*exp(-mu.*rs)./rs;
fk = abs(fftshift(fftn(f)));
fk2 = 1./(K2 + mu*mu);