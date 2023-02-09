function matrix_km(A,b,m)
    n = length(b)
    Km = zeros(n,m);
    Km[:,1] = b/norm(b);
    w = b;
    for i in 2:m
        w = A*w;
        w = w/norm(w);
        Km[:,i] = w;
    end
    return Km;
end
