function A = create_laplacian3d(n, stsz, dx, periodic)

offsets = [];
coeffs = [];

if (stsz == 3)
    offsets = [0, 1];
    coeffs = [-2.0, 1];
elseif (stsz == 7)
    offsets = [0, 1, 2, 3];
    coeffs = [-490.0/180.0, 270.0/180.0, -27.0/180.0, 2.0/180.0];
end   

T1N = coeffs(1)*speye(n);

for i = 2:(stsz+1)/2
    T1N = T1N + coeffs(i)*spdiags(ones(n,1),offsets(i),n,n);
    T1N = T1N + coeffs(i)*spdiags(ones(n,1),-offsets(i),n,n);
    if (periodic)
        T1N = T1N + spdiags(coeffs(i)*ones(n,1), n-offsets(i), n, n);
        T1N = T1N + spdiags(coeffs(i)*ones(n,1), -n+offsets(i), n, n);
    end
end

A = kron(speye(n),kron(speye(n),T1N)) + kron(T1N,kron(speye(n),speye(n))) + kron(speye(n),kron(T1N,speye(n)));
A = A./(dx)^2;
