function Q = gram_schmidt(A)

Q = zeros(size(A));
Q(:,1) = A(:,1);
Q(:,1) = Q(:,1)/norm(Q(:,1));
for i = 2:size(A,2)
    Q(:,i) = A(:,i);
    for j = 1:i-1
        Q(:,i) = Q(:,i) - (Q(:,j)'*Q(:,i))*Q(:,j);
    end
    Q(:,i) = Q(:,i)/norm(Q(:,i));    
end

end

