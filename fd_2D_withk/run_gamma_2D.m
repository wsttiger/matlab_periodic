clear all;

npts = 21;
nstates = 6;
L = 5.0;
[X,Y,dx,KX,KY,dk] = make_periodic_2d_grid(npts, L);
% need K.^2 for BSH
K2 = KX.^2 + KY.^2;

% 1-body potential
alpha = 2.5;
V = -alpha*(cos(2*pi*X/L).*cos(2*pi*Y/L) + 1);

% create an initial guess
psi = randn(npts*npts,nstates);
e = -3.0*ones(nstates,1);

% kinetic energy
T = -0.5*sparse(create_laplacian2d(npts, 7, dx, 1));
Vm = spdiags(V(:),0,npts*npts,npts*npts);
H = T + Vm;


alpha = 1;
for iter = 1:400
    Vpsi = Vm*psi;
    Vk = zeros(npts,npts);
    rnorm = 0.0;
    new_psi = zeros(npts*npts,nstates);
    for i = 1:nstates
        new_psi(:,i) = apply_BSH_2D_FFT(Vpsi(:,i), K2, e(i));
        psi(:,i) = alpha*new_psi(:,i) + (1-alpha)*psi(:,i);
    end

    Hmo = psi'*H*psi;
    [vv,dd] = eig(Hmo);
    e = diag(dd)
    psi = psi*vv;
end
