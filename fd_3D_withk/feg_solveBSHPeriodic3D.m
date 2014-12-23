clear all;

% lattice constant
L = 10.0;
% number of states
nstates = 50;

% kpoints
nkpts = 2;
q = (2*pi/L)*linspace(0.21,0.21,nkpts);
nkpts = nkpts-1;
q = q(1:nkpts);
[QX,QY,QZ] = meshgrid(q,q,q);

npts = 22;
eps = 0.0;
x = linspace(-L/2+eps, L/2+eps, npts)';
delx = x(2) - x(1);
npts = npts-1;
x = x(1:npts);
[X,Y,Z] = meshgrid(x,x,x);

% potential
V = -10.0*ones(npts*npts*npts,1);
alpha = 2.5;
% V = -alpha*(cos(2*pi*X/L).*cos(2*pi*Y/L).*cos(2*pi*Z/L) + 1) - 3.0;
% V = V(:);

Ek = zeros(nstates,nkpts,nkpts,nkpts);
Ek_exact = zeros(nkpts,nkpts,nkpts);
psi = zeros(npts*npts*npts,nstates,nkpts,nkpts,nkpts);
[Dx,Dy,Dz] = create_Dx_3D(npts, delx);
T0 = -0.5*create_laplacian3d(npts,7,delx, 1);
for ik = 1:nkpts
    for jk = 1:nkpts
        for kk = 1:nkpts
            fprintf('%d     %d     %d\n', ik, jk, kk);
            % kinetic energy
            qx = q(ik); qy = q(jk); qz = q(kk); q2 = qx^2 + qy^2 + qz^2;
            fprintf('%10.5f     %10.5f     %10.5f     %10.5f\n', qx, qy, qz, q2);
            Tk = T0 - i*qx*Dx - i*qy*Dy - i*qz*Dz + 0.5*q2*speye(npts*npts*npts);
            Hk = Tk + spdiags(V(:),0,npts*npts*npts,npts*npts*npts);
            [Cs,Es] = eigs(Hk,nstates,'sr');
            psi(:,:,ik,jk,kk) = Cs;
            Ek(:,ik,jk,kk) = diag(Es);
            Ek_exact(ik,jk,kk) = -10.0 + 0.5*q2;
        end
    end
end

% psi_k = psi(:,1:nstates,:);
% e = Ek(1:nstates);
psi_k = rand(npts*npts*npts,nstates,nkpts,nkpts,nkpts);
e = -3.0*ones(nstates,nkpts,nkpts,nkpts);
e = randn(nstates,nkpts,nkpts,nkpts);

% setup k-grid
sidesz = (npts-1)/2;
kx = 2*pi*(-sidesz:sidesz)'./delx./npts;
[KX,KY,KZ] = meshgrid(kx,kx,kx);
K2 = KX.^2 + KY.^2 + KZ.^2;

niter = 400;
es = zeros(nstates,nkpts,nkpts,nkpts,niter);
for iter = 1:niter
    % the accumulation of the r2 over each state
    r2norm_k = 0.0;
    for ik = 1:nkpts
        for jk = 1:nkpts
            for kk = 1:nkpts
                % corrections to the energies
                edelta = zeros(nstates,1);
                % the new state
                new_psi_k = zeros(npts*npts*npts, nstates);
                % k-dependent Hamiltonian
                qx = q(ik); qy = q(jk); qz = q(kk); q2 = qx^2 + qy^2 + qz^2;
                Tk = - i*qx*Dx - i*qy*Dy - i*qz*Dz + 0.5*q2*speye(npts*npts*npts);
                H = T0 + Tk + spdiags(V(:),0,npts*npts*npts,npts*npts*npts);
                for ist = 1:nstates        
                    % compute a shift if one is needed
                    shift = 0.0;
                    if (e(ist,ik,jk,kk) > -0.0001)
                        shift = -0.05 + e(ist,ik,jk,kk);
                    end
        
                    % apply the potential (with the shift if needed)
                    Vpsi = (V - shift*ones(npts*npts*npts,1)).*psi_k(:,ist,ik,jk,kk);
                    Vpsi = Vpsi - Tk*psi_k(:,ist,ik,jk,kk);
        
                    % apply the convolution BSH operator
                    mu = sqrt(2.0*(shift-e(ist,ik,jk,kk)));
                    new_psi_k(:,ist) = apply_BSH_3D_FFT(Vpsi, K2, mu);
                end
        
                % Compute the k-dependent Hamiltonian in the MO basis
                psi_aug = [psi_k(:,:,ik,jk,kk) new_psi_k];
                psi_aug = orth(psi_aug);
                Hmo_k = psi_aug'*H*psi_aug;
                Hmo_k = 0.5*(Hmo_k + Hmo_k');
                [vv,dd] = eig(Hmo_k);
        
                % update the energies and transform the orbitals accordingly
                t1 = diag(dd);
                es(:,ik,jk,kk,iter) = e(:,ik,jk,kk)-t1(1:nstates);
                norm(es(:,ik,jk,kk,iter))
                e(:,ik,jk,kk) = t1(1:nstates);
%                 es(:,ik,jk,kk,iter) = t1(1:nstates);
                t2 = psi_aug*vv;
                r2norm_k = r2norm_k + norm(psi_k(:,:,ik,jk,kk)'*psi_k(:,:,ik,jk,kk) - t2(:,1:nstates)'*t2(:,1:nstates));
                psi_k(:,:,ik,jk,kk) = t2(:,1:nstates);
            end
        end
    end
    fprintf('rnorm:  %15.8f\n', r2norm_k); 
end

Ek_display = zeros(size(Ek));
for ik = 1:nkpts
    for jk = 1:nkpts
        Ek_display(:,ik,jk,kk) = sort(Ek(:,ik,jk,kk));
    end
end



