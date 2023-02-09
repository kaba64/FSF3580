#using Pkg
#Pkg.add("MAT")
using MAT
using LinearAlgebra
using Plots
using LaTeXStrings

B=matread("Bwedge.mat")["B"];
B_eigvals = eigvals(B);
m = size(B_eigvals);
#Q,H=arnoldi(A,b,m);
scatter(real(B_eigvals),imag(B_eigvals),xlim=(-50,10),markersize = 6,marker=:circle,legend=false,color="MediumPurple");
scatter!([real(B_eigvals[1]),real(B_eigvals[m[1]-1])],[imag(B_eigvals[1]),imag(B_eigvals[m[1]-1])],markersize = 6,xlim=(-50,10),markershape=:star4,legend=false,color="red");
fig = scatter!([real(B_eigvals[1]),real(B_eigvals[m[1]])],[imag(B_eigvals[1]),imag(B_eigvals[m[1]])],markersize = 6,xlim=(-50,10),markershape=:star4,legend=false,color="red");
xlabel!("Re")
ylabel!("Im");
savefig(fig,"excercise_4a.png")
