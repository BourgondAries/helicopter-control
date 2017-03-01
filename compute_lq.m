Q = diag([1,1,0,0]);
R = diag(1);
[k,s,e] = dlqr(A1,B1,Q,R);
