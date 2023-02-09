function Rayleigh_method(A,v,n,nmax,lam_max)
    mu::Float64 = v'*A*v;
    v_r = v;
    max_iter = 1;
    tol = 1e-15;
    for i in 2:nmax
        B = A-mu*Matrix{Float64}(I,3,3);
        w = B\v;
        v = w/norm(w);
        #mu = v'*A*v;
        if (abs(lam_max - v'*A*v)<tol || norm(v-v_r)<tol)
            return lam_r, v_r, max_iter;
        end
        v_r = v;
        mu = v'*A*v;
        lam_r[i] = lam_max - mu;
        max_iter+= 1;

    end
    return lam_r, v_r, max_iter;
end
