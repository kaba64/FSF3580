using LinearAlgebra
using LaTeXStrings
n = 8;
A = [2 1 0 0 0 0 0 0 ;
1 2 -1 0 0 0 0 0;
0 -1 2 1 0 0 0 0
0 0 1 2 1 0 0 0
0 0 0 1 2 1 0 0
0 0 0 0 1 2 1 0
0 0 0 0 0 1 2 1
0 0 0 0 0 0 1 2];
b = [1 ;1 ;0 ;0 ;0 ;0 ;0 ;0];
km = zeros(n,4);
km[:,1] = b;
for k in 2:4
    w = A*km[:,k-1]
    km[:,k] = w;
end
display(km);
