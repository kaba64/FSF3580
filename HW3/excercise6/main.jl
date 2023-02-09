#using Pkg
#Pkg.add("MatrixDepot")
using MatrixDepot
using QuadGK
using LinearAlgebra
using Random
using SparseArrays
using LaTeXStrings
include("PN.jl")

function P(τ)
    return exp(Matrix(-τ.*A))*B*exp(Matrix(τ.*A))
end

A=matrixdepot("neumann",20);
A=A-A';
A=A./(2*norm(A,1));
B=sprandn(20^2,20^2, 0.05);
τ= 1;

println("Time spent in  quadgk : ");
@time quadgk(P,0,τ);
println("Time spent in  the new algorithm : ");
@time PN(4,A,B,τ);
(Pold,erro) = quadgk(P,0,τ);
pnew = PN(4,A,B,τ);
println("The relative error : ", norm((Pold-pnew)/Pold));
