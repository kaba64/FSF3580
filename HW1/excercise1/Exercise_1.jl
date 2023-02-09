using Pkg
Pkg.add("Plots")
using Plots
using LinearAlgebra
using LaTeXStrings
include("power_method.jl")
include("Rayleigh_qoutient.jl")

A  = [1 2 3;2 2 2;3 2 9];
λex, vex = eigen(A);
lam_max = maximum(abs.(λex));
gr();
n = 3;
nmax = 15;
v = [1; 1 ;1];
v_p = zeros(n);
v_r = zeros(n);
lam_p = Array{Float64}(undef, nmax);
lam_r = Array{Float64}(undef, nmax);
v = v/norm(v);
lam_p[1] = lam_max-v'*A*v;
lam_r[1] = lam_max-v'*A*v;
max_iter = 1;
lam_p, v_p, max_iter = power_method(A,v,nmax,lam_max);
plot([1:max_iter],abs.(lam_p[1:max_iter]),yaxis=:log,label="PM");
max_iter = 1;
lam_r, v_r, max_iter = Rayleigh_method(A,v,n,nmax,lam_max);
#display(lam_r[1:max_iter]);
plot!([1:max_iter],abs.(lam_r[1:max_iter]),yaxis=:log,label="RQM");
# non-symmetric A
A  = [1 2 4;2 2 2;3 2 9];
λex, vex = eigen(A);
lam_max = maximum(abs.(λex));
lam_r[1] = lam_max-v'*A*v;
max_iter = 1;
lam_r, v_r, max_iter = Rayleigh_method(A,v,n,nmax,lam_max);
t = plot!([1:max_iter],abs.(lam_r[1:max_iter]),yaxis=:log,label="RQM for non-symmetric A");
xlabel!("Number of iteration")
ylabel!(L"log(\lambda^{n}-\lambda)");
t;
#savefig(t,"excercise_1.png")
