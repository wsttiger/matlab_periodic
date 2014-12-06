clear all;

% lattice constant
npts = 132;
L = 5.0;
eps = 1.e-6;
x = linspace(-L/2+eps, L/2+eps, npts);
delx = x(2) - x(1);
npts = npts-1;
x = x(1:npts);

% setup k-grid
sidesz = (npts-1)/2;
kx = 2*pi*(-sidesz:sidesz)./delx./npts;

mu = 3.0;
f = exp(-mu*x.^2);
fk = abs(fftshift(fft(f)));
fk2 = (1./sqrt(2*mu))*exp(-kx.^2/4.0/mu);