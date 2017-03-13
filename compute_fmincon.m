% nonlcon is a function returning a [c, ceq] which is a vector of values for a given input x.
N = 15;
q2 = 1;
q1 = q2;
G = [kron(eye(N),[1 0 0 0 0 0]) kron(eye(N),[q1 q2])];
fun = @(x) x'*G*x;

alpha1 = K_1*K_pp;
alpha2 = K_1*K_pd;

% Discrete time system model. x = [lambda r p p_dot]'
A1 = [0   1   0       0        0         0;
      0   0   -K_2    0        0         0;
      0   0   0       1        0         0;
      0   0   -alpha1 -alpha2  0         0;
      0   0   0       0        0         1;
      0   0   0       0        -K_3*K_ep -K_3*K_ed]*h + eye(4);
B1 = [0   0   0       alpha1   0         0;
      0   0   0       0        0         K_3*K_ep]'*h;
fmincon(fun,x0,[],[],Aeq,beq,lb,ub,nonlcon);
