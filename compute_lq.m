Q = diag(1,0,0,0);
R = diag(1);
N = diag(0);
A = (A1 - I)/h;
B = B1;
dlqr(A,B,Q,R,N);
