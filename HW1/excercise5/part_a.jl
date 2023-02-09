using MAT
using LinearAlgebra
using Plots
using LaTeXStrings

include("arnoldi.jl")

function sort_eig(Heigvals,Heigvec,k)
    temp_val = Heigvals;
    temp_vec = Heigvec;
    i=argmax(abs.(temp_val));
    eigvanew = temp_val[i];
    eigvenew = temp_vec[:,i];
    deleteat!(temp_val, findall(x->x==temp_val[i],temp_val));
    temp_vec = temp_vec[:, 1:end .!= i];
    for j in 2:k
        i=argmax(abs.(temp_val));
        eigvanew = hcat(eigvanew,temp_val[i]);
        eigvenew = hcat(eigvenew,temp_vec[:,i]);
        deleteat!(temp_val, findall(x->x==temp_val[i],temp_val));
        temp_vec = temp_vec[:, 1:end .!= i];
    end
    return eigvanew,eigvenew;
end

B=matread("Bwedge.mat")["B"];
Beigvals = eigvals(B);
m = 20;
k = 10;
restart = 100;
n = size(Beigvals)
b = randn(n[1])
b = b./norm(b);
Q,H = arnoldi(B,b,m);
Heigvals, Heigvec= eigen(H[1:m,1:m]);
RitazVal = zeros(Complex,k,restart);
RitazVal[1:k,1] = Heigvals[1:k];
eigvanew,eigvenew = sort_eig(Heigvals,Heigvec,k);
column = ones(k);
for j in 2:restart
    b = (Q[:,1:m]*eigvenew)*column;
    b = b./norm(b);
    global Q,H = arnoldi(B,b,m);
    Heigvals, Heigvec= eigen(H[1:m,1:m]);
    global RitazVal[1:k,j] = Heigvals[1:k];
    global eigvanew,eigvenew = sort_eig(Heigvals,Heigvec,k);
end

plot([1:restart],real(RitazVal[1,1:restart]),legend=false);
for j = 2:k
    global f1 = plot!([1:restart],real(RitazVal[j,1:restart]),legend=false);
end
f1;
xlabel!("Restart")
ylabel!("Real part of the Ritaz values");
savefig(f1,"excercise_5a_r_2.png")

plot([1:restart],imag(RitazVal[1,1:restart]),legend=false);
for j = 2:k
    global f1 = plot!([1:restart],imag(RitazVal[j,1:restart]),legend=false);
end
f1;
xlabel!("Restart")
ylabel!("Imaginary part of the Ritaz values");
savefig(f1,"excercise_5a_i_2.png")
