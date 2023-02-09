#using Pkg
#Pkg.add("SparseArrays")
using Plots
using LinearAlgebra
using Random
using SparseArrays
using BenchmarkTools
using LaTeXStrings

include("GMRES.jl")
gr();
n=100; m=100;
alpha_value=[1,5,10,100];
Random.seed!(0);
b = rand(n);
#alpha=5;

for alpha in alpha_value
    A = sprand(n,n,0.5);
    A = A + alpha*sparse(I, n, n);
    A=A/norm(A,1);
    x,error_x,residual = GMRES(A,b,m)
    plot([1:m],error_x,yaxis=:log,label="Error : α=$alpha",legend=:bottomleft);
    scatter!([1:m],residual,yaxis=:log,label="Residula : α=$alpha",legend=:bottomleft);
    xlabel!("m")
    ylabel!("Error");
    savefig("excercise_1a-α=$alpha.png")
end
