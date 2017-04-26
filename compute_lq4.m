Q = diag([1,0,0,0,0,0]);
R = diag([q1, q2]);
[k,s,e] = dlqr(A1,B1,Q,R);
