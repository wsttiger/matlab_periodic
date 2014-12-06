function [Ek, psik] = diag_H_k_2D(npts,nstates,dx,H,kpoints)

nkpts = size(kpoints,1);
psik = zeros(npts,npts,nstates,nkpts);
Ek = zeros(nstates,nkpts);
for ik = 1:nkpts
    Tk = create_k_T_matrix(npts, dx, kpoints(ik,:));
    Hk = H + Tk;
    [Cs,Es] = eigs(Hk,nstates,'sr');
    psik(:,:,:,ik) = reshape(Cs, [npts npts nstates]);
    Ek(:,ik) = diag(Es);
end

end

