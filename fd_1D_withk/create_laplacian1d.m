function A = create_laplacian1d(n, stsz, dx)

offsets = [];
coeffs = [];

if (stsz == 3)
    offsets = [0, 1];
    coeffs = [-2.0, 1];
elseif (stsz == 7)
    offsets = [0, 1, 2, 3];
    coeffs = [-490.0/180.0, 270.0/180.0, -27.0/180.0, 2.0/180.0];
end   

A = spalloc(n, n, 10*n);

% off-diagonal
for i = 2:(stsz+1)/2
    A = A + spdiags(coeffs(i)*ones(n,1), offsets(i), n, n);
end

% this makes it periodic
for i = 2:(stsz+1)/2
    A = A + spdiags(coeffs(i)*ones(n,1), n-offsets(i), n, n);
end

% make symmetric
A = A'+A;

% diagonal
A = A + spdiags(coeffs(1)*ones(n,1), offsets(1), n, n);

A = A./(dx)^2;
