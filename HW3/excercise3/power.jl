function power(A,N)
  P = A
  for k in 2:N
    P = P*A;
  end
  return P
end
