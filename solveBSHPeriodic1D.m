clear all;

% lattice constant
npts = 22;
L = 5.0;
eps = 0.0;
x = linspace(-L/2+eps, L/2+eps, npts);
delx = x(2) - x(1);
npts = npts-1;
x = x(1:npts);

% 1-body potential
alpha = 2.5;
expnt = 1;
% V = -alpha*(cos(2*pi*x/ L) + 1);
% V = -alpha*(exp(-expnt*x.^2) + exp(-expnt*(x-L).^2) + exp(-expnt*(x+L).^2));
V = zeros(1,npts);
for j = -5:5
    V = V - alpha*(exp(-expnt*(x-j*L).^2));
end
% create an initial guess
a = 10.0;
psi = exp(-a*x.^2);
psi = psi./norm(psi(:));
e = -3.0;

% kinetic energy
T = -0.5*sparse(create_laplacian1d(npts, 7, delx));
H = T + spdiags(V(:),0,npts,npts);

% setup k-grid
sidesz = (npts-1)/2;
kx = 2*pi*(-sidesz:sidesz)./delx./npts;
k2 = kx.^2;

alpha = 1.0;
for iter = 1:40
    Vpsi = V.*psi;
    mu = sqrt(-2.0*e);
    
    Vpsi = fftshift(fft(Vpsi));
    Vpsi = ifftshift(Vpsi./(k2 +mu.^2));
    new_psi = -2.0*real(ifft(Vpsi));
    
    psi = alpha*new_psi + (1-alpha)*psi;
    psi = psi/norm(psi(:));
    e = psi(:)'*H*psi(:);
    
    %update the energy
    r = psi-new_psi;
    rnorm = norm(r(:));
    fprintf('rnorm:  %15.8f      e:  %15.8f      eig error: %15.8f\n', ...
        rnorm, e, norm(H*psi(:) - e*psi(:))*sqrt(delx));
end
