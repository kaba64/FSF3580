function GCM(A::Matrix,b,m::Int64)
    n = length(b);
    x = zeros(n);
    r = zeros(n);
    p = zeros(n);
    r = b;
    p = r;
    for i in 1:m
        q1 = r'*r;
        q2 = A*p;
        α = q1/(p'*q2);
        x = x+α*p;
        r = r-α*q2;
        β = α = (r'*r)/q1;
        p = r+β*p;
    end
    return x;
end
