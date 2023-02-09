function PN(n,A,B,τ)
    Gi = B;
    global temp = τ;
    P = temp*B;
    for i = 1:n
        Gi = Gi*A-A*Gi;
        global temp = (τ/(i+1))*temp;
        P += temp*Gi;
    end
    return P

end
