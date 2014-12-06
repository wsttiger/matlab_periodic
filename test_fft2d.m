clear all;

% lattice constant
npts = 132;
L = 5.0;
eps = 1.e-6;
x = linspace(-L/2+eps, L/2+eps, npts);
y = linspace(-L/2+eps, L/2+eps, npts);
delx = x(2) - x(1);
npts = npts-1;
x = x(1:npts); y = y(1:npts);

% setup k-grid
sidesz = (npts-1)/2;
kx = 2*pi*(-sidesz:sidesz)./delx./npts;
ky = 2*pi*(-sidesz:sidesz)./delx./npts;
[KX,KY] = meshgrid(kx,ky);
K2 = KX.^2 + KY.^2;

[X,Y] = meshgrid(x,y);
rs = sqrt(X.^2 + Y.^2);

% f = exp(-mu*(X.^2+Y.^2));
% fk = abs(fftshift(fft2(f)));
% fk = fk/norm(fk);
% fk2 = exp(-K2/4.0/mu);
% fk2 = fk2/norm(fk2);

mu = 3.0;
f = exp(-mu*rs)./rs;
fk = abs(fftshift(fft2(f)));