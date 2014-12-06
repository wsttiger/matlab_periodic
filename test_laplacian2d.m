clear all;

npts = 102;
L = 5.0;
x = linspace(-L/2,L/2,npts);
delx = x(2)-x(1);
npts = npts-1;
x = x(1:npts);

[X,Y] = meshgrid(x,x);
V = cos(2.0*pi*X/L).*cos(2.0*pi*Y/L);
rho = -(8.0*pi*pi/L/L)*cos(2.0*pi*X/L).*cos(2.0*pi*Y/L);

A = create_laplacian2d(npts,7,delx,1);

