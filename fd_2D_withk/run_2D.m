clear all;

% lattice constant
L = 5.0;
% number of states
nstates = 100;

% kpoints
nkpts = 5;
q = (2*pi/L)*linspace(0.0,1.0,nkpts);
nkpts = nkpts-1;
q = q(1:nkpts);
[QX,QY] = meshgrid(q,q);

npts = 22;
eps = 0.0;
x = linspace(-L/2+eps, L/2+eps, npts)';
delx = x(2) - x(1);
npts = npts-1;
x = x(1:npts);
[X,Y] = meshgrid(x,x);

% potential
alpha = 2.5;
V = -alpha*(cos(2*pi*X/L).*cos(2*pi*Y/L) + 1) - 0.0;
V = V(:);

Ek = zeros(nstates,nkpts,nkpts);
psi = zeros(npts*npts,nstates,nkpts*nkpts);
[Dx,Dy] = create_Dx_2D(npts, delx);
T0 = -0.5*create_laplacian2d(npts,7,delx, 1);
for ik = 1:nkpts
    for jk = 1:nkpts
        fprintf('%d     %d     \n', ik, jk);
        % kinetic energy
        qx = q(ik); qy = q(jk); q2 = qx^2 + qy^2;
        Tk = T0 - i*qx*Dx - i*qy*Dy + 0.5*q2*speye(npts*npts);
        H = Tk + spdiags(V(:),0,npts*npts,npts*npts);
%        [Cs,Es] = eigs(H,nstates,'sr');
%        psi(:,:,ik,jk) = Cs;
%        Ek(:,ik,jk) = diag(Es);
        [Cs,Es] = eig(full(H),'vector');
        psi(:,:,ik,jk) = Cs(:,1:nstates);;
        Ek(:,ik,jk) = Es(1:nstates);
    end
end

% more realistic initial guess (with 4 orbitals)
%orb1s = exp(-a1*R2)
psi_k = zeros(npts*npts,nstates,nkpts,nkpts);
psi_k = psi(:,1:nstates,:,:);
e = Ek + 1.0*randn(nstates,nkpts,nkpts);
psi_k = rand(npts*npts,nstates,nkpts,nkpts);
% e = -3.0*ones(nstates,nkpts,nkpts);
% e = randn(nstates,nkpts,nkpts);

% setup k-grid
sidesz = (npts-1)/2;
kx = 2*pi*(-sidesz:sidesz)'./delx./npts;
[KX,KY] = meshgrid(kx,kx);
K2 = KX.^2 + KY.^2;

niter = 20;
energies = zeros(nstates,niter);
for iter = 1:niter
    % the accumulation of the r2 over each state
    r2norm_k = 0.0;
    for ik = 1:nkpts
        for jk = 1:nkpts
            % corrections to the energies
            edelta = zeros(nstates,1);
            % the new state
            new_psi_k = zeros(npts*npts, nstates);
            % k-dependent Hamiltonian
            qx = q(ik); qy = q(jk); q2 = qx^2 + qy^2;
            Tk = - i*qx*Dx - i*qy*Dy + 0.5*q2*speye(npts*npts);
            H = T0 + Tk + spdiags(V(:),0,npts*npts,npts*npts);
            for ist = 1:nstates        
                % compute a shift if one is needed
                shift = 0.0;
                if (e(ist,ik,jk) > -0.0001)
                    shift = -0.05 + e(ist,ik,jk);
                end
    
                % apply the potential (with the shift if needed)
                Vpsi = (V - shift*ones(npts*npts,1)).*psi_k(:,ist,ik,jk);
                Vpsi = Vpsi - Tk*psi_k(:,ist,ik);
    
                % apply the convolution BSH operator
                mu = sqrt(2.0*(shift-e(ist,ik,jk)));
                new_psi_k(:,ist) = apply_BSH_2D_FFT(Vpsi, K2, mu);
            end
    
            % Compute the k-dependent Hamiltonian in the MO basis
            psi_aug = [psi_k(:,:,ik,jk) new_psi_k];
            psi_aug = orth(psi_aug);
            Hmo_k = psi_aug'*H*psi_aug;
            Hmo_k = 0.5*(Hmo_k + Hmo_k');
            [vv,dd] = eig(Hmo_k);
    
            % update the energies and transform the orbitals accordingly
            t1 = diag(dd);
            e(:,ik,jk) = t1(1:nstates);
            norm(e(:)-Ek(:))/nstates/nkpts/nkpts
            % energies(:,iter) = t1(1:nstates);
            t2 = psi_aug*vv;
            r2norm_k = r2norm_k + norm(psi_k(:,:,ik,jk)'*psi_k(:,:,ik,jk) - t2(:,1:nstates)'*t2(:,1:nstates));
            psi_k(:,:,ik,jk) = t2(:,1:nstates);
        end
    end
    fprintf('rnorm:  %15.8f\n', r2norm_k); 
end

Ek_display = zeros(size(Ek));
for ik = 1:nkpts
    for jk = 1:nkpts
        Ek_display(:,ik,jk) = sort(Ek(:,ik,jk));
    end
end



