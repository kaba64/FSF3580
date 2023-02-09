using LinearAlgebra

function QR_dgs(Q,w,k)
    q = Q[:,1:k-1];
    h = q'*w
    z = w-q*h;
    g = q'*z;
    z = z-q*g;
    h = h+g;
    β = norm(z);
    return h,β,z;
end
