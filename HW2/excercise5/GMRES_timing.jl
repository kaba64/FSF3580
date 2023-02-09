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
    #println("elapsed time: ", t, " seconds")
    return t
end

function GMRES_timing(A::SparseMatrixCSC{Complex{Float64},Int64},b::Array{Float64,1},m::Int64)
    #nano = 1.e-9;
    #tt =  Int(time_ns());
    n=length(b);
    Q=im*zeros(n,m+1);
    H=im*zeros(m+1,m);
    e_1 = zeros(Int64,m+1,1);
    x = zeros(n);
    b_norm = norm(b);
    Q[:,1]=b/b_norm;
    e_1[1]=1;
    residual_time_gmres = zeros(m,2);
    t = 0.0;
    for k=1:m
        tic();
        w=A*Q[:,k]; # Matrix-vector product with last element
        # Orthogonalize w against columns of Q.
        # Implement this function or replace call with code for orthogonalizatio
        h,β,z=hw1_dgs(Q,w,k);
        #Put Gram-Schmidt coefficients into H
        H[1:(k+1),k]=[h;β];
        # normalize
        Q[:,k+1]=z/β;
        z_s = (H[1:(k+1),1:k]\e_1[1:(k+1)])*b_norm;
        x = Q[:,1:k]*z_s;
        residual_time_gmres[k,2] = norm(A*x-b)/b_norm;
        dt = toc();
        t +=dt;
        residual_time_gmres[k,1] = t;
    end
    return x, residual_time_gmres;
end
