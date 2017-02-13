% Find a minimum for a specified problem:
% min x over 1/2x^THx+f^T*x st A*x <= b, Aeq*x=beq, lb<=x<=ub
ub = [ones(1,100)*Inf ones(1,200)*30*pi/180]';
lb = [ones(1,100)*-Inf ones(1,200)*-30*pi/180]';
q = [0.1 1 10];
A = zeros(300);
A(1,1) = 1;
b = zeros(300,1);
b(1) = pi;
for i=q
	H = blkdiag(eye(100), eye(100)*i, zeros(100))*2;
	x = quadprog(H,[],[],[],A,b,lb,ub);
	fprintf('Optimal value %d\n', x(end));
end
