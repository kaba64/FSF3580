using MAT
using LinearAlgebra
using SparseArrays
using BenchmarkTools
using Plots
using LaTeXStrings
using Arpack
include("GMRES.jl")
include("CGNEM.jl")
m = 10;
B=matread("Bwedge2.mat")["B"];
b=matread("Bwedge2.mat")["b"];
x,residual = GMRES(B,b[:,1],m);
x_cgne,residual_gcne = GCNEM(B,b[:,1],m);
plot([1:m],residual,yaxis=:log,label="Residula GMRES",legend=:bottomleft);
scatter!([1:m],residual_gcne,yaxis=:log,label="Residula GCNE",legend=:bottomleft);
xlabel!("m")
ylabel!("Residual");
savefig("excercise_5a-1.png")
