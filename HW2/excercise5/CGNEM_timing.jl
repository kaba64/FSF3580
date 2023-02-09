using LinearAlgebra
using BenchmarkTools

include("gram_schmidt_double.jl")
# The tic toc functions are taken from
#https://github.com/JuliaLang/julia/blob/903644385b91ed8d95e5e3a5716c089dd1f1b08a/base/util.jl#L9-L13

function tic()
    t0 = time_ns()
    task_local_storage(:TIMERS, (t0, get(task_local_storage(), :TIMERS, ())))
    return t0
end
function toq()
    t1 = time_ns()
    timers = get(task_local_storage(), :TIMERS, ())
    if timers === ()
        error("toc() without tic()")
    end
    t0 = timers[1]::UInt64
    task_local_storage(:TIMERS, timers[2])
    (t1-t0)/1e9
end
function toc()
    t = toq()
    return t
end

function CGNEM_timing(A::SparseMatrixCSC{Complex{Float64},Int64},b,m)
    #nano = 1.e-9;
    #tt =  Int(time_ns());
    n = length(b);
    x_cgne = zeros(n);
    #residual_gcne=zeros(m);
    b_norm = norm(b);
    r = A'*b;
    p = r;
    residual_time_cgne = zeros(m,2);
    t = 0.0;
    for i in 1:m
        tic();
        q1 = r'*r;
        q2 = A'*(A*p);
        α = q1/(p'*q2);
        x_cgne = x_cgne+α*p;
        r = r-α*q2;
        β = (r'*r)/q1;
        p = r+β*p;
        residual_time_cgne[i,2] = norm(A*x_cgne-b)/b_norm;
        dt = toc();
        t +=dt;
        residual_time_cgne[i,1] = t;
    end
    return x_cgne, residual_time_cgne;
end
