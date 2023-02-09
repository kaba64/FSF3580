using LinearAlgebra

include("QR_gram_schmidt_double.jl")

function QR_Factorization(A)

    nx,ny=size(A);
    Q=copy(A);
    R=zeros(nx,ny);
    R[1,1] = norm(A[:,1])
    Q[:,1]=Q[:,1]/R[1,1];
    for k=2:ny
        h,β,z=QR_dgs(Q,A[:,k],k);
        R[1:k,k]=[h;β];
        Q[:,k]=z/β;
    end
    return Q,R
end
