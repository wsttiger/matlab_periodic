clear all;

npts = 222;
L = 5.0;
x = linspace(-L/2,L/2,npts);
delx = x(2)-x(1);
npts = npts-1;
x = x(1:npts);

[X,Y] = meshgrid(x,x);
V = cos(4.0*pi*X/L).*cos(2.0*pi*Y/L);

DxV = -(4.0*pi/L)*sin(4.0*pi*X/L).*cos(2.0*pi*Y/L);
DyV = -(2.0*pi/L)*cos(4.0*pi*X/L).*sin(2.0*pi*Y/L);

[Dx,Dy] = create_Dx_2D(npts,delx);

DxV2 = reshape(Dx*V(:), [npts npts]);
DyV2 = reshape(Dy*V(:), [npts npts]);


