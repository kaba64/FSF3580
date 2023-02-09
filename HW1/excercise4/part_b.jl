using Pkg
#Pkg.add("HDF5")
using HDF5
using MAT
using LinearAlgebra
using Plots
using LaTeXStrings

function alpha(B_eigvals,xr,xi,i,radius)
    x1 = real(B_eigvals[i])-xr;
    x2 = imag(B_eigvals[i])-xi;
    α = radius/sqrt(x1*x1+x2*x2);
    return α;
end

function  circlue(x1,x2,radius)
    θ = 0:(2*pi/1000):2pi;
    θ = θ[1:end-1];
    center = (x1 + x2*im) * ones(1000);
    z = center + radius * exp.(θ*im);
    return z;
end

B=matread("Bwedge.mat")["B"];
B_eigvals = eigvals(B);
m = size(B_eigvals);
z1 = circlue(-5,0.4,14.5);
z2 = circlue(-23.5,9,26);
z3 = circlue(-23,-6,25.5);
α1 = alpha(B_eigvals,-5,0.4,1,14.5);
α2 = alpha(B_eigvals,-23.5,9,m[1]-1,26);
α3 = alpha(B_eigvals,-23,-6,m[1]-1,25.5);
println(α1);
println(α2);
println(α3);
scatter(real(B_eigvals),imag(B_eigvals),xlim=(-60,20),markersize = 6,marker=:circle,labels="Eigenval.",legend=:bottomleft,color="MediumPurple");
plot!(z1; xlabel="Re", ylabel="Im",labels=L"\alpha_{1}",legend=:topright, aspectratio=:equal)
plot!(z2; xlabel="Re", ylabel="Im",labels=L"\alpha_{2}",legend=:topright, aspectratio=:equal)
fig = plot!(z3; xlabel="Re", ylabel="Im",labels=L"\alpha_{3}",legend=:topright, aspectratio=:equal)
xlabel!("Re")
ylabel!("Im");
fig;
savefig(fig,"excercise_4b.png")
