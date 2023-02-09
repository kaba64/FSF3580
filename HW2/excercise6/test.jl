using Plots

n = 10;
z=zeros(n);
x = 0.008;
for i in 1:n
    z[i] = (1+2*n)/(1-x*n)
end
plot([1:n],z);
display(z);
