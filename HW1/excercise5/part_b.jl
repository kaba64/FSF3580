using MAT
using LinearAlgebra
using Plots
using LaTeXStrings
using Random

include("restart_arnoldi.jl")

B=matread("Bwedge.mat")["B"];
Beigvals = eigvals(B);
Beigvecs = eigvecs(B);
n = size(B)
b = randn(n[1])
b = b./norm(b);

k = 10; # The Ritz vectors corresponding to the k largest Ritz values
p = 20; # For computing an Arnoldi factorization of length k+p
number_elements = 50;
V,H,r,restart_count, RitazVal = restart_arnoldi(B,k,p,1e-10,b,number_elements);

plot([1:restart_count],real(RitazVal[1,1:restart_count]),legend=false);
for j = 2:k
    global f1 = plot!([1:restart_count],real(RitazVal[j,1:restart_count]),legend=false);
end
f1;
xlabel!("Restart")
ylabel!("Real part of the Ritaz values");
savefig(f1,"excercise_5b_r_2.png")

plot([1:restart_count],imag(RitazVal[1,1:restart_count]),legend=false);
for j = 2:k
    global f2 = plot!([1:restart_count],imag(RitazVal[j,1:restart_count]),legend=false);
end
f2;
xlabel!("Restart")
ylabel!("Imaginary part of the Ritaz values");
savefig(f2,"excercise_5b_i_2.png")
