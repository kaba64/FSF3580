using LinearAlgebra

include("Householder.jl")

function Hessenberg_red(A)
    n=size(A,1)
    for k=1:n-2
        x=A[k+1:n,k];
        z,u = Householder(x); # Eq (3-4) in lecture notes
        # Compute P*A
        A[k+1:n,k:n]=A[k+1:n,k:n]-2*u*(u'*A[k+1:n,k:n])
        # Compute P*A*P'
        A[1:n,k+1:n]=A[1:n,k+1:n]-2*(A[1:n,k+1:n]*u)*u'
    end
    return H_red=A
end
