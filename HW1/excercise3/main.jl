#Pkg.add("BenchmarkTools")
using Plots
using LinearAlgebra
using MatrixDepot, Random, BenchmarkTools
using Random
include("arnoldi.jl")
include("matrix_km.jl")

Random.seed!(0);
dim = 15
A = matrixdepot("wathen",dim,dim);
b=randn(size(A,2));
m_iter = 80;
gr();
fig = plot(reuse=false);
for m in 1:m_iter
    Q,H=arnoldi(A,b,m);
    #println("1:st should be zero = ", norm(Q*H-A*Q[:,1:m]));
    #println("2:nd Should be zero = ", norm(Q'*Q-I));
    h = H[1:m,1:m];
    Km = matrix_km(A,b,m);
    KmT = transpose(Km);
    AKm = inv(KmT*Km)*(KmT*A*Km);
    λexh, vexh = eigen(h);
    λexk, vexk = eigen(AKm);
    global fig1 = scatter!(m*ones(m),real(λexh),markersize = 6,marker=:circle,legend=false,color="MediumPurple");
    global fig2 = scatter!(m*ones(m),real(λexk),markersize = 5,ylim=(-1000, 1000),markershape=:star4,legend=false,color = "red");
end
fig1;
fig2;
xlabel!("m")
ylabel!("Real part of eigenval. approx.");
#savefig(fig1,"excercise_3b.png")
