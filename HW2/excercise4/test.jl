k = 10.5/2.0;
error = 10e-7;
x = log(0.5*error);
y = log((sqrt(k)-1)/(sqrt(k)+1));
n = x/y;
display(n);
n = n+3;
z = ((sqrt(k)-1)/(sqrt(k)+1))^n
display(2*z);
println(n);
