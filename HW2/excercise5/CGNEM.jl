function GCNEM(A::SparseMatrixCSC{Complex{Float64},Int64},b,m)
    n = length(b);
    x_cgne = zeros(n);
    residual_gcne=zeros(m);
    b_norm = norm(b);
    r = A'*b;
    p = r;
    for i in 1:m
        q1 = r'*r;
        q2 = A'*(A*p);
        α = q1/(p'*q2);
        x_cgne = x_cgne+α*p;
        r = r-α*q2;
        β = (r'*r)/q1;
        p = r+β*p;
        residual_gcne[i] = norm(A*x_cgne-b)/b_norm
    end
    return x_cgne,residual_gcne;
end
