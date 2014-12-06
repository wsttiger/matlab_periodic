function [X,Y,dx,KX,KY,dk] = make_periodic_2d_grid(npts, L)
x = linspace(-L/2, L/2, npts+1);
dx = x(2) - x(1);
x = x(1:npts);
[X,Y] = meshgrid(x,x);
sidesz = (npts-1)/2;
kx = 2*pi*(-sidesz:sidesz)./dx./npts;
dk = kx(2)-kx(1);
[KX,KY] = meshgrid(kx,kx);
end

