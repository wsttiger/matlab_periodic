clear all;

% lattice constant
npts = 42;
nstates = 7;
L = 5.0;
eps = 1.e-6;
x = linspace(-L/2+eps, L/2+eps, npts);
y = linspace(-L/2+eps, L/2+eps, npts);
z = linspace(-L/2+eps, L/2+eps, npts);
delx = x(2) - x(1);
npts = npts-1;
x = x(1:npts); y = y(1:npts); z = z(1:npts);

% 1-body potential
alpha = 2.5;
[X,Y,Z] = meshgrid(x,y,z);
V = -alpha*(cos(2*pi*X/ L).*cos(2*pi*Y/ L).*cos(2*pi*Z/L) + 1);

% create an initial guess
a = 10.0;
psi = randn(npts,npts,npts,nstates);
% psi(:,:,:,1) = exp(-3.0*(X.^2+Y.^2+Z.^2));
% psi(:,:,:,2) = X.*exp(-3.0*(X.^2+Y.^2+Z.^2));
% psi(:,:,:,3) = 7.*exp(-3.0*(X.^2+Y.^2+Z.^2));
% psi(:,:,:,4) = Z.*exp(-3.0*(X.^2+Y.^2+Z.^2));
% psi(:,:,:,5) = X.*Y.*exp(-3.0*(X.^2+Y.^2+Z.^2));
% psi(:,:,:,6) = X.*Z.*exp(-3.0*(X.^2+Y.^2+Z.^2));
% psi(:,:,:,7) = Z.*Y.*exp(-3.0*(X.^2+Y.^2+Z.^2));
psi = reshape(psi,[npts*npts*npts,nstates]);
psi = randn(npts*npts*npts,nstates);
for i = 1:nstates
    psi(:,i) = psi(:,i)/norm(psi(:,i));
end
e = -3.0*ones(nstates,1);

% kinetic energy
T = -0.5*sparse(create_laplacian3d(npts, 7, delx, 1));
V = spdiags(V(:),0,npts*npts*npts,npts*npts*npts);
H = T + V;
clear T;

% Hmo = psi'*H*psi
% [vv,dd] = eig(Hmo);
% e = diag(dd)

% setup k-grid
sidesz = (npts-1)/2;
kx = 2*pi*(-sidesz:sidesz)./delx./npts;
ky = 2*pi*(-sidesz:sidesz)./delx./npts;
kz = 2*pi*(-sidesz:sidesz)./delx./npts;
[KX,KY,KZ] = meshgrid(kx,ky,kz);
K2 = KX.^2 + KY.^2 + KZ.^2;
clear KX KY KZ;

alpha = 1.0;
for iter = 1:145
    Vpsi = V*psi;
    
    Vk = zeros(npts,npts,npts);
    rnorm = 0.0;
    new_psi = zeros(npts*npts*npts,nstates);
    for i = 1:nstates
        mu = sqrt(-2.0*e(i));
        Vk(:) = Vpsi(:,i);
        Vk = fftshift(fftn(Vk));
        Vk = ifftshift(Vk./(K2 +mu.^2));
        t1 = -2.0*real(ifftn(Vk));
        new_psi(:,i) = t1(:);
        psi(:,i) = alpha*new_psi(:,i) + (1-alpha)*psi(:,i);
    end
    
    Sr = psi'*new_psi
    Hmo = psi'*H*psi
    [vv,dd] = eig(Hmo);
    e = diag(dd)
    psi = psi*vv;
    
    %update the energy
%     fprintf('rnorm:  %15.8f      e:  %15.8f      eig error: %15.8f\n', ...
%         rnorm, e, norm(H*psi(:) - e*psi(:))*sqrt(delx));
%    fprintf('e:  %15.8f\n\n', e);
end
