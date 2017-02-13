% Find a minimum for a specified problem:
% min x over 1/2x^THx+f^T*x st A*x <= b, Aeq*x=beq, lb<=x<=ub
ub = [ones(1,100)*Inf ones(1,100)*30*pi/180]';
lb = [ones(1,100)*-Inf ones(1,100)*-30*pi/180]';
q = 10.0;
H = [eye(100) zeros(100); zeros(100) eye(100)*q]*2;
% H*[ones(100,1); ones(100,1)*2]
x = quadprog(H,[],[],[],[],[],lb,ub);
