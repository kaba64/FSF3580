using LinearAlgebra

include("gram_schmidt_double.jl")

function GMRES(A::SparseMatrixCSC{Float64,Int64},b::Array{Float64,1},m::Int64)

    n=length(b);
    Q=zeros(n,m+1);
    H=zeros(m+1,m);
    error_x=zeros(m);
    residual=zeros(m);
    e_1 = zeros(Int64,m+1,1);
    x = zeros(n);
    x_exact = zeros(n);
    x_exact = A\b;
    x_norm = norm(x_exact)
    b_norm = norm(b);
    Q[:,1]=b/b_norm;
    e_1[1]=1;
    for k=1:m
        w=A*Q[:,k]; # Matrix-vector product with last element
        # Orthogonalize w against columns of Q.
        # Implement this function or replace call with code for orthogonalizatio
        h,β,z=hw1_dgs(Q,w,k);
        #Put Gram-Schmidt coefficients into H
        H[1:(k+1),k]=[h;β];
        # normalize
        Q[:,k+1]=z/β;
        z_s = (H[1:(k+1),1:k]\e_1[1:(k+1)])*b_norm;
        x = Q[:,1:k]*z_s;
        residual[k] = norm(A*x-b)/b_norm;
        error_x[k] = norm(x_exact-x)/x_norm;
    end
    return x,error_x,residual
end
