using LinearAlgebra

include("gram_schmidt_double.jl")

function Benchmark_GMRES(A::SparseMatrixCSC{Float64,Int64},b::Array{Float64,1},m::Int64)

    n=length(b);
    Q=zeros(n,m+1);
    H=zeros(m+1,m);
    e_1 = zeros(Int64,m+1,1);
    x = zeros(n);
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
    end
    z_s = (H[1:(m+1),1:m]\e_1[1:(m+1)])*b_norm;
    x = Q[:,1:m]*z_s;
    return x
end
