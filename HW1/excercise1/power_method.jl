function power_method(A,v,nmax,lam_max)
    max_iter = 1;
    v_p = v;
    tol = 1e-15;
    for i in 2:nmax
        v = A*v;
        v = v/norm(v);
        if (abs(lam_max - v'*A*v)<tol || norm(v-v_p)<tol)
            return lam_p, v_p, max_iter;
        end
        v_p = v;
        lam_p[i] = lam_max - v'*A*v;
        max_iter+= 1;

    end
    return lam_p, v_p, max_iter;
end
