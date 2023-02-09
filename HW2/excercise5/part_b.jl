using MAT
using LinearAlgebra
using SparseArrays
using BenchmarkTools
using Plots
using LaTeXStrings

function  circlue(x1,x2,radius)
    θ = 0:(2*pi/1000):2pi;
    θ = θ[1:end-1];
    center = (x1 + x2*im) * ones(1000);
    z = center + radius * exp.(θ*im);
    return z;
end

B=matread("Bwedge2.mat")["B"];
b=matread("Bwedge2.mat")["b"];
λexB, vexB = eigen(Matrix(B));
z1 = circlue(36.5,0,21.5);
α1 = round(21.5/36.5;digits = 2);
scatter(real(λexB),imag(λexB),xlim=(0,80),markersize = 6,marker=:circle,labels="eigenvalues of B",legend=:topright,color="MediumPurple");
fig1 = plot!(z1; xlabel="Re", ylabel="Im",labels="α=$α1",legend=:topright, aspectratio=:equal)
xlabel!("Re")
ylabel!("Im");
savefig(fig1,"excercise_5b-1.png")
cod_b = cond(Matrix(B))
α2 = (cod_b-1)/(cod_b+1);
display(α1);
display(α2);
