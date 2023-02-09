using LinearAlgebra
using Random
using BenchmarkTools
include("schur_parlett.jl")

A = [1 4 5;3 -1 3;-1 4 5]
f=z->sin(z)
F=schur_parlett(A,f)

display(F);
