clear all;

L = 5.0;
npts = 81;
nstates = 6;

% kpts = [0.0 0.0;
%         0.0 0.5;
%         0.5 0.0;
%         0.5 0.5];
kpts = [0.5 0.5];
kpts = kpts * (2*pi/L);
nkpts = size(kpts,1);
 
[X,Y,dx,KX,KY,dk] = make_periodic_2d_grid(npts, L);
% need K.^2 for BSH
K2 = KX.^2 + KY.^2;

% 1-body potential
alpha = 2.5;
%V = -alpha*(cos(2*pi*X/L).*cos(2*pi*Y/L) + 1);
V = -10.0*ones(npts,npts);

% create an initial guess
psi = randn(npts*npts,nstates,nkpts);
for ik = 1:nkpts
    for ist = 1:nstates
        tmp = reshape(exp(-3.0*(X.^2 + Y.^2)), [npts*npts 1]);
        tmp = tmp/norm(tmp);
        psi(:,ist,ik) = tmp;
    end
end
e = -3.0*ones(nstates,nkpts);

% kinetic energy
T = -0.5*sparse(create_laplacian2d(npts, 7, dx, 1));
Vm = spdiags(V(:),0,npts*npts,npts*npts);
H = T + Vm;

[Ek, psi_k] = diag_H_k_2D(npts,nstates,dx,H,kpts);

psi = reshape(psi_k, [npts*npts nstates nkpts]);
e = Ek;

alpha = 1;
for iter = 1:0
    e_old = e;
    for ik = 1:nkpts
        psik = psi(:,:,ik);
        Tk = create_k_T_matrix(npts, dx, kpts(ik,:));
        phase1 = exp(i*(kpts(ik,1)*X + kpts(ik,2)*Y));
        phase2 = exp(-i*(kpts(ik,1)*X + kpts(ik,2)*Y));
        Vpsik = Vm*psik;
        Vpsik = Vpsik.*phase1(:);
        new_psik = zeros(npts*npts,nstates);
        for ist = 1:nstates
            new_psik(:,ist) = phase2(:).*apply_BSH_2D_FFT(Vpsik(:,ist), K2, e(ist,ik));
            psik(:,ist) = alpha*new_psik(:,ist) + (1-alpha)*psik(:,ist);
        end
%         Hmo_k = psik'*(H + Tk)*psik;
%         Hmo_k = 0.5*(Hmo_k + Hmo_k');
%         [vv,dd] = eig(Hmo_k);
%         e(:,ik) = diag(dd)
%         psi(:,:,ik) = psik*vv;
    end
    
%     fprintf('norm of residual energies: %15.8e\n', norm(e_old - e));;

end
