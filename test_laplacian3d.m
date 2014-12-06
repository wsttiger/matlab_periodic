clear all;

npts = 102;
L = 5.0;
x = linspace(-L/2,L/2,npts);
delx = x(2)-x(1);
npts = npts-1;
x = x(1:npts);

[X,Y,Z] = meshgrid(x,x,x);
V = cos(2.0*pi*X/L).*cos(2.0*pi*Y/L).*cos(2.0*pi*Z/L);
rho = -(12.0*pi*pi/L/L)*cos(2.0*pi*X/L).*cos(2.0*pi*Y/L).*cos(2.0*pi*Z/L);

A = create_laplacian3d(npts,7,delx,1);
rho2 = A*V(:);

fprintf('Norm of the error is: %15.8e\n', norm(rho2 - rho(:))*delx*delx*delx);