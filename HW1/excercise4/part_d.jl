using MAT
using LinearAlgebra
using Plots
using LaTeXStrings

include("arnoldi.jl")

B=matread("Bwedge.mat")["B"];
Beigvals = eigvals(B);
m = size(Beigvals);
n = size(B)
b = randn(m[1])
b = b./norm(b);
σ = -9.8+2*im;
i=argmin(abs.(Beigvals-σ*ones(m[1])));
xclosest=Beigvals[i];
σ1 = -9.8+1.5*im;
B = inv(B-σ1*Matrix(I,n[1],n[2]));
iteration = [10,20,30];
for k in 1:length(iteration)
    m = iteration[k];
    Q,H = arnoldi(B,b,m);
    Ritz_values = zeros(Complex,m);
    RitzEigvals = eigvals(H[1:m,1:m]);
    for i in 1:m
        global RitzEigvals[i] = 1.0/RitzEigvals[i]+σ1;
    end
    inew=argmin(abs.(xclosest*ones(m)-RitzEigvals));
    println("After $m the error is : ");
    display(abs.(RitzEigvals[inew]-xclosest));
end
