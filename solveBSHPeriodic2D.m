clear all;

% lattice constant
npts = 42;
nstates = 1;
L = 5.0;
x = linspace(-L/2, L/2, npts);
delx = x(2) - x(1);
npts = npts-1;
x = x(1:npts);

% 1-body potential
[X,Y] = meshgrid(x,x);
alpha = 2.5;
V = -alpha*(cos(2*pi*X/L).*cos(2*pi*Y/L) + 1);

% expnt = 0.1;
% V = zeros(npts,npts);
% for n1 = -5:5
%     for n2 = -5:5
%         V = V - alpha*exp(-expnt*(X-n1*L).^2).*exp(-expnt*(Y-n2*L).^2);
%     end
% end

% create an initial guess
psi = randn(npts*npts,nstates);
e = -3.0*ones(nstates,1);

% kinetic energy
T = -0.5*sparse(create_laplacian2d(npts, 7, delx, 1));
Vm = spdiags(V(:),0,npts*npts,npts*npts);
H = T + Vm;

% setup k-grid
sidesz = (npts-1)/2;
kx = 2*pi*(-sidesz:sidesz)./delx./npts;
[KX,KY] = meshgrid(kx,kx);
K2 = KX.^2 + KY.^2;

[Cs,Es] = eigs(H,nstates,'sa');

psi = Cs;
e = diag(Es);

alpha = 1;
for iter = 1:5
    Vpsi = Vm*psi;
    Vk = zeros(npts,npts);
    rnorm = 0.0;
    new_psi = zeros(npts*npts,nstates);
    for i = 1:nstates
        mu = sqrt(-2.0*e(i));
        Vk(:) = Vpsi(:,i);
        Vk = fftshift(fft2(Vk));
        Vk = ifftshift(Vk./(K2 +mu.^2));
        t1 = -2.0*real(ifft2(Vk));
        new_psi(:,i) = t1(:);
        psi(:,i) = alpha*new_psi(:,i) + (1-alpha)*psi(:,i);        
    end

    % DEBUG
    my_error = Cs-psi;
    fprintf('My GOD: norm of my_error is %15.8e\n', norm(my_error));
    
    Hmo = psi'*H*psi;
    [vv,dd] = eig(Hmo);
    e_old = e;
    e = diag(dd);
    resid_e = e - e_old;
%     fprintf('norm of residual energies: %15.8e\n', norm(resid_e));
    psi = psi*vv;
end
