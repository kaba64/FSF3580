#using Pkg
#Pkg.add("SparseArrays")
using Plots
using LinearAlgebra
using Random
using SparseArrays
using BenchmarkTools
using LaTeXStrings
include("GMRES.jl")

function  circlue(x1,x2,radius)
    θ = 0:(2*pi/1000):2pi;
    θ = θ[1:end-1];
    center = (x1 + x2*im) * ones(1000);
    z = center + radius * exp.(θ*im);
    return z;
end

gr();
n=100;
alpha_value=[1,5,10,100];
m_iter = [100,60,30,20];
Random.seed!(0);
b = rand(n);
z = zeros(3,4);
con_rate = zeros(4);
# z1 = c_real, z2 = c_imag, z3 = radius
z[1,1] = 0.00035; z[2,1] =0.0; z[3,1] =0.0015;
z[1,2] = 0.0017; z[2,2] =0.0; z[3,2] =0.0013;
z[1,3] = 0.00285; z[2,3] =0.0; z[3,3] =0.0011;
z[1,4] = 0.008; z[2,4] =0.0; z[3,4] =0.00031;

for k in 1:length(alpha_value)
    alpha = alpha_value[k];
    A = sprand(n,n,0.5);
    A = A + alpha*sparse(I, n, n);
    A=A/norm(A,1);
    λexa = eigvals(Matrix(A));
    scatter(real(λexa),imag(λexa),markersize = 6,marker=:circle,label="Eigenvalues of A : α=$alpha",color="MediumPurple");
    f = plot!(circlue(z[1,k],z[2,k],z[3,k]); xlabel="Re", ylabel="Im",labels=latexstring("\$\\alpha_{$(k)}\$"),legend=:topright, aspectratio=:equal)
    xlabel!("Re")
    ylabel!("Im");
    savefig(f,"excercise_1b-α=$alpha.png")
end
con_rate[1] = z[3,1]/z[1,1]; con_rate[2] = z[3,2]/z[1,2];
con_rate[3] = z[3,3]/z[1,3]; con_rate[4] = z[3,4]/z[1,4];
# **************** convergence rate **********************
for k in 1:length(alpha_value)
    alpha = alpha_value[k];
    m = m_iter[k];
    A = sprand(n,n,0.5);
    A = A + alpha*sparse(I, n, n);
    A=A/norm(A,1);
    x,error_x,residual = GMRES(A,b,m)
    λexa = eigvals(Matrix(A));
    i=argmax(abs.(λexa));
    coe = (abs(λexa[i]-z[1,k])+ z[3,k])/(abs(λexa[i]));
    plo = coe*(con_rate[k].^((1:length(residual)).-1));
    plot([1:m],plo,yaxis=:log,label="Convergence factor : α=$alpha",legend=:topright);
    f = scatter!([1:m],residual,yaxis=:log,label="Residula : α=$alpha",legend=:topright);
    xlabel!("m")
    ylabel!("Error");
    savefig("excercise_1b2-α=$alpha.png")
end
