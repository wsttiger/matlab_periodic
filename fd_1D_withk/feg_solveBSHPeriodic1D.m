clear all;

% lattice constant
L = 5.0;

% kpoints
nkpts = 40;
kpts = (2*pi/L)*linspace(-0.5,0.5,nkpts);

npts = 22;
eps = 0.0;
x = linspace(-L/2+eps, L/2+eps, npts)';
delx = x(2) - x(1);
npts = npts-1;
x = x(1:npts);

% potential
V = -10.0*ones(npts,1);

Ek = zeros(npts,nkpts);
psi = zeros(npts,npts,nkpts);
for ik = 1:nkpts
    % kinetic energy
    Tk = -0.5*create_laplacian1d(npts,7,delx) - i*kpts(ik)*create_D_matrix(npts, 7, delx, 1) + 0.5*kpts(ik)^2*speye(npts);
    H = Tk + spdiags(V(:),0,npts,npts);
    [vv,dd] = eig(full(H));
    psi(:,:,ik) = vv;
    Ek(:,ik) = diag(dd);
end

nstates = 7;
% psi_k = psi(:,1:nstates,:);
% e = Ek(1:nstates);
psi_k = rand(npts,nstates,nkpts);
e = -3.0*ones(nstates,nkpts);
e = randn(nstates,nkpts);

% setup k-grid
sidesz = (npts-1)/2;
kx = 2*pi*(-sidesz:sidesz)'./delx./npts;
k2 = kx.^2;

niter = 10;
energies = zeros(nstates,niter);
for iter = 1:niter
    for ik = 1:nkpts
        % residual
        r = zeros(npts,nstates);
        % residual of the orbital charge density --> psi.^2
        r2 = zeros(npts,nstates);
        % the accumulation of the r2 over each state
        r2norm = 0.0;
        % corrections to the energies
        edelta = zeros(nstates,1);
        % the new state
        new_psi_k = zeros(npts, nstates);
        for ist = 1:nstates        
            % compute a shift if one is needed
            shift = 0.0;
            if (e(ist,ik) > -0.0001)
                shift = -0.05 + e(ist,ik);
            end

            % apply the potential (with the shift if needed)
            Vpsi = (V - shift*ones(npts,1)).*psi_k(:,ist,ik);
            Tk = -i*kpts(ik)*create_D_matrix(npts, 7, delx, 1) + 0.5*kpts(ik)^2*speye(npts);
            Vpsi = Vpsi - Tk*psi_k(:,ist,ik);

            % apply the convolution BSH operator
            mu = sqrt(2.0*(shift-e(ist,ik)));
            Vk = fftshift(fft(Vpsi));
            Vk = ifftshift(Vk./(k2 +mu.^2*ones(npts,1)));
            new_psi_k(:,ist) = -2.0.*ifft(Vk);
        end

        % Compute the k-dependent Hamiltonian in the MO basis
        Tk = -0.5*create_laplacian1d(npts,7,delx) - i*kpts(ik)*create_D_matrix(npts, 7, delx, 1) + 0.5*kpts(ik)^2*speye(npts);
        H = Tk + spdiags(V(:),0,npts,npts);
        psi_aug = [psi_k(:,:,ik) new_psi_k]
        rank(psi_aug)
        psi_aug = orth(psi_aug);
        Hmo_k = psi_aug'*H*psi_aug;
        Hmo_k = 0.5*(Hmo_k + Hmo_k')
        [vv,dd] = eig(Hmo_k);

        % update the energies and transform the orbitals accordingly
        t1 = diag(dd);
        e(:,ik) = t1(1:nstates);
        energies(:,iter) = t1(1:nstates);
        t2 = psi_aug*vv;
        psi_k(:,:,ik) = t2(:,1:nstates);
    end
    %fprintf('rnorm:  %15.8f\n', r2norm); 
end
