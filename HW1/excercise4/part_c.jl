using MAT
using LinearAlgebra
using Plots
using LaTeXStrings

include("arnoldi.jl")

B=matread("Bwedge.mat")["B"];
B_eigvals = eigvals(B);
m = size(B_eigvals);
b = randn(m[1])
b = b./norm(b);
#********** Ritz values **************
n_ritz = 50;
Ritz_values = zeros(Complex,n_ritz);
error = 1.0;
m = 0;
eps = 1.e-10;
while error>eps
    global m = m+1;
    Q,H = arnoldi(B,b,m);
    RitzEigvals = eigvals(H[1:m,1:m]);
    global error = abs.(RitzEigvals[1]-B_eigvals[1])
    Ritz_values[m] = RitzEigvals[1];
end
new = B_eigvals[1]*ones(m)
pl = abs.(Ritz_values[1:m]-new)
t = plot([1:m],pl,xlim=(1,m+4),ylim=(1e-13,log(pl[1])),yaxis=:log,labels="Convergence of Ritz values for the Arnoldi method");
xlabel!("m")
ylabel!("Eigenvalue error");
savefig(t,string("excercise_4c_convergence.png"));
#*****************************************************
#m_num = [2, 4, 8, 10, 20, 30, 40];

#for m in m_num
#        Q,H = arnoldi(B,b,m);
#        H_eigvals = eigvals(H[1:m,1:m]);
#        scatter(real(B_eigvals),imag(B_eigvals),xlim=(-60,20),markersize = 6,marker=:circle,labels="Eigenvalues of B",legend=:bottomleft,color="MediumPurple");
#        t = scatter!(real(H_eigvals),imag(H_eigvals),labels=string("Eigenvalues of ","H_{",m,",",m,"}"));
#        xlabel!("Re")
#        ylabel!("Im");
#        savefig(t,string("excercise_4c_m =",m,".png"));
#end
