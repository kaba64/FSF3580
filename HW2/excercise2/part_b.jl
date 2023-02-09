#using Pkg
#Pkg.add("Optim")
using LinearAlgebra
using Optim
using SparseArrays
using LaTeXStrings
include("CGM.jl")


A = [2 1 0 0 0 0 0 0 ;
1 2 -1 0 0 0 0 0;
0 -1 2 1 0 0 0 0
0 0 1 2 1 0 0 0
0 0 0 1 2 1 0 0
0 0 0 0 1 2 1 0
0 0 0 0 0 1 2 1
0 0 0 0 0 0 1 2];
b = [1 ;1 ;0 ;0 ;0 ;0 ;0 ;0];
# ************ from part (a)  α=0, β=-1, γ=6 *********************
C=[1 0 0 6
1 0 -1 0
0 1 0 0
0 0 1 7
0 0 0 1
0 0 0 0
0 0 0 0
0 0 0 0];
# CG is the minimizer of ||Ax-b||_A^{-1}; x=Cz
residual_inverse_A(z) = (A*C*z-b)'*(A\(A*C*z-b))
r_bfgs = optimize(residual_inverse_A, [1.0; 1; 1; 1], BFGS())
z=r_bfgs.minimizer;
x_optimized = C*z;
m = 4;
x = GCM(A,b,m);
println("Residual of the solution obtained by the optimization and CGM")
display(norm(x-x_optimized));
