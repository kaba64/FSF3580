using LinearAlgebra
using Random
using BenchmarkTools
include("naive_hessenberg_red.jl")
include("Hessenberg_red.jl")
include("alpha_example.jl")

iter_m = [10,100,200,300,400]
times = zeros(3,5);

count = 1
for m in iter_m
    global A = alpha_example(1,m);
    t1 = @belapsed naive_hessenberg_red(A);
    t2 = @belapsed Hessenberg_red(A);
    global times[1,count] = m;
    global times[2,count] = t1;
    global times[3,count] = t2;
    global count = count + 1
    #display(norm(H_red-H));
end
display(times);
