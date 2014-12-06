clear all;
close all;

for npts = [50]
% amplitude
A = 12.5;
% size of the box
L = 5.0
% number of plane waves
%npts = 15;

% generate g-vectors
gvecs = [-(npts-1)/2:(npts-1)/2]*2*pi/L;
% generate potential for -A(cos(2*pi*n*x/L) + 1)
Vconst = -A*L*diag(ones(1,npts));
Vcos = -0.5*A*L*(diag(ones(npts-1,1),1) + diag(ones(npts-1,1),-1));
V = Vconst + Vcos;
% generate T
T = 0.5*L*diag(gvecs.^2);
H = T + V;
evH = eig(H)/L;

if (size(evH,1) >= 10)
    evH(1:10)
else
    evH
end
end