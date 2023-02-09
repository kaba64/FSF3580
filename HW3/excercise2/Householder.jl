using LinearAlgebra

function Householder(x)
    n=length(x);

    norm_x = norm(x);
    rho = -sign(x[1]);
    z = x;
    z[1] = z[1]-rho*norm_x;
    u = z/norm(z);
    return z,u;
end
