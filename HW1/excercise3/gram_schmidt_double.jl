using LinearAlgebra

function hw1_dgs(Q,w,k)
    q = Q[:,1:k];
    h = q'*w
    z = w-q*h;
    g = q'*z;
    z = z-q*g;
    h = h+g;
    β = norm(z);
    return h,β,z;
end
