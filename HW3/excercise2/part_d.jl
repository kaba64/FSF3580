using LinearAlgebra
using Random
using BenchmarkTools

include("QR.jl")

eps = 0.4;
sigma = 0.0;
eps_value = [0.4,0.1,0.01,1.e-3,1.e-4,1.e-5,1.e-6,1.e-7,1.e-8,1.e-9,1.e-10];
h_values = zeros(3,length(eps_value));
count = 0;
for eps in eps_value
    global count = count + 1;
    global A = [3 2 ;eps 1];
    global h_values[1,count] = eps;
    sigma = 0
    Q, R = QR_Factorization(A-sigma*I);
    H = R*Q + sigma*I;
    global h_values[2,count] = abs(H[2,1]);
    sigma = A[2,2];
    Q, R = QR_Factorization(A-sigma*I);
    H = R*Q + sigma*I;
    global h_values[3,count] = abs(H[2,1]);
end
display(h_values);
