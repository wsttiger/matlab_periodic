function Tk = create_k_T_matrix(n, dx, kpoint)

kx = kpoint(1);
ky = kpoint(2);
k2 = kx*kx + ky*ky;
D = create_D_matrix(n,7,dx,1);
Dx = kron(D,speye(n));
Dy = kron(speye(n),D);
Tk = -i*(kx*Dx + ky*Dy) + 0.5*k2*speye(n*n);
end

