using LinearAlgebra
using Random
using BenchmarkTools
using Plots

include("schur_parlett.jl")
include("power.jl")

A = rand(100,100)+im*rand(100,100)
A=A/norm(A);

iter_n = [5,25,50,100,150,175,250,300,350];
timing = zeros(3,length(iter_n));
N = 5;
# Schur-Parlett method implementation
f=z->z^N;
b_s = @benchmarkable schur_parlett(A,f);
tune!(b_s);
g_t_b_s = run(b_s);
display("Schur-Parlett method")
display((mean(g_t_b_s).time)/1.e6);
# Naive implementation
b_n = @benchmarkable power(A,N);
tune!(b_n);
g_t_b_n = run(b_n);
display("Naive implementation")
display((mean(g_t_b_n).time)/1.e6);
# *** for plotting use this part
#timing = [5 25 50 100 150 175 250 300 350;
#11.5606 11.8796 11.6169 11.5603 12.119 12.7141 13.4390 11.7338 12.1099;
#0.45 2.8219 5.8703 11.647 16.8479 20.2743 30.9794 35.5683 43.1119];
#scatter(timing[1,:],timing[2,:],markersize = 6,marker=:circle,label="Schur-Parlett method",color="MediumPurple");
#plot!(timing[1,:],0.12*timing[1,:],markersize = 6,marker=:star,label="theoretical",color="green");
#f = plot!(timing[1,:],timing[3,:],markersize = 6,marker=:circle,label="Naive implementation",color="red");
#xlabel!("N")
#ylabel!("time(s)");
#savefig(f,"excercise_3b.png")
