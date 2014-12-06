function Dx = create_D_matrix(n, stsz, dx, periodic)

offsets = [];
coeffs = [];

if (stsz == 3)
    offsets = [0, 1];
    coeffs = [0.0, 1.0/2.0];
elseif (stsz == 7)
    offsets = [0, 1, 2, 3];
    coeffs = [0.0, 45.0/60.0, -9.0/60.0, 1.0/60.0];
end   

Dx = coeffs(1)*speye(n);

for i = 2:(stsz+1)/2
    Dx = Dx + coeffs(i)*spdiags(ones(n,1),offsets(i),n,n);
    Dx = Dx - coeffs(i)*spdiags(ones(n,1),-offsets(i),n,n);
    if (periodic)
        Dx = Dx - coeffs(i)*spdiags(ones(n,1), n-offsets(i), n, n);
        Dx = Dx + coeffs(i)*spdiags(ones(n,1), -n+offsets(i), n, n);
    end
end

Dx = Dx./dx;

end

