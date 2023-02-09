using MAT
using LinearAlgebra
using SparseArrays
using BenchmarkTools
using Plots
using LaTeXStrings
include("GMRES_timing.jl")
include("CGNEM_timing.jl")

m = 100;
B=matread("Bwedge2.mat")["B"];
b=matread("Bwedge2.mat")["b"];


x,residual_time_gmres = GMRES_timing(B,b[:,1],m);
x_cgne, residual_time_cgne = CGNEM_timing(B,b[:,1],m);

scatter(residual_time_gmres[1:m,1],residual_time_gmres[1:m,2],yaxis=:log,label="Residula GMRES",legend=:bottomleft);
scatter!(residual_time_cgne[1:m,1],residual_time_cgne[1:m,2],yaxis=:log,label="Residula CGNE",legend=:bottomleft);
xlabel!("time (s)");
ylabel!("Residual");
savefig("excercise_5a-2.png");
