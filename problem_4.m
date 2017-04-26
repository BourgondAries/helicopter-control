% TTK4135 - Helicopter lab
% Hints/template for problem 2.
% Updated spring 2017, Andreas L. Fl?ten

%% Initialization and model definition
init07; % NB: Change this to the init file corresponding to your helicopter
delta_t = 0.25; % sampling time
h = delta_t;
q1 = 1;
q2 = 1;

alpha1 = K_1*K_pp;
alpha2 = K_1*K_pd;

% Discrete time system model. x = [lambda r p p_dot]'
A1 = [0   1   0       0        0         0;
      0   0   -K_2    0        0         0;
      0   0   0       1        0         0;
      0   0   -alpha1 -alpha2  0         0;
      0   0   0       0        0         1;
      0   0   0       0        -K_3*K_ep -K_3*K_ed]*h + eye(6);
B1 = [0   0   0       alpha1   0         0;
      0   0   0       0        0         K_3*K_ep]'*h;

% Number of states and inputs
mx = size(A1,2);                        % Number of states (number of columns in A)
mu = size(B1,2);                        % Number of inputs (number of columns in B)

% Time horizon and initialization
N  = 80;                                % Time horizon for states
M  = N;                                 % Time horizon for inputs
z  = zeros(N*mx+M*mu,1);                % Initialize z for the whole horizon
z0 = z;                                 % Initial value for optimization

% Initial values
x1_0 = pi;                              % Lambda
x2_0 = 0;                               % r
x3_0 = 0;                               % p
x4_0 = 0;                               % p_dot
x5_0 = 0;                               % e
x6_0 = 0;                               % e_dot
x0 = [x1_0 x2_0 x3_0 x4_0 x5_0 x6_0]';  % Initial values


% Bounds
ul      = -30*pi/180;                   % Lower bound on control -- u1
uu      = 30*pi/180;                    % Upper bound on control -- u1

xl      = -Inf*ones(mx,1);              % Lower bound on states (no bound)
xu      = Inf*ones(mx,1);               % Upper bound on states (no bound)
xl(3)   = ul;                           % Lower bound on state x3
xu(3)   = uu;                           % Upper bound on state x3

% Generate constraints on measurements and inputs
[vlb,vub]       = genbegr2(N,M,xl,xu,ul,uu); % hint: genbegr2
vlb(N*mx+M*mu)  = 0;                    % We want the last input to be zero
vub(N*mx+M*mu)  = 0;                    % We want the last input to be zero
vlb(N*mx+M*mu-1)  = 0;                  % We want the last input to be zero
vub(N*mx+M*mu-1)  = 0;                  % We want the last input to be zero

% Generate the matrix Q and the vector c (objecitve function weights in the QP problem)
Q1 = zeros(mx,mx);
Q1(1,1) = 1;                            % Weight on state x1
Q1(2,2) = 0;                            % Weight on state x2
Q1(3,3) = 0;                            % Weight on state x3
Q1(4,4) = 0;                            % Weight on state x4
Q1(5,5) = 0;                            % Weight on state x5
Q1(6,6) = 0;                            % Weight on state x6
P1 = diag([q1, q2]);                    % Weight on input
Q = genq2(Q1,P1,N,M,mu);                % Generate Q
c = zeros(N*mx+M*mu,1);                 % Generate c

fun = @(x) x'*Q*x;

%% Generate system matrixes for linear model
Aeq = gena2(A1,B1,N,mx,mu);           % Generate A, hint: gena2
beq = zeros(mx*N,1);        	  % Generate b
beq(1:mx) = A1*x0; % Initial value

%% Solve QP problem with linear model
options = optimoptions(@fmincon,'MaxIter',27000,'MaxFunEvals',270000);
tic;
z = fmincon(fun,z0,[],[],Aeq,beq,vlb,vub,@nonlcon,options);
% [z,lambda] = quadprog(Q,c,[],[],Aeq,beq,vlb,vub); % hint: quadprog
t1=toc;

% Calculate objective value
phi1 = 0.0;
PhiOut = zeros(N*mx+M*mu,1);
for i=1:N*mx+M*mu
  phi1=phi1+Q(i,i)*z(i)*z(i);
  PhiOut(i) = phi1;
end

%% Extract control inputs and states
u1  = [z(N*mx+1:2:N*mx+M*mu-1);z(N*mx+M*mu-1)]; % Control input from solution
u2  = [z(N*mx+2:2:N*mx+M*mu);z(N*mx+M*mu)]; % Control input from solution

x1 = [x0(1);z(1:mx:N*mx)];              % State x1 from solution
x2 = [x0(2);z(2:mx:N*mx)];              % State x2 from solution
x3 = [x0(3);z(3:mx:N*mx)];              % State x3 from solution
x4 = [x0(4);z(4:mx:N*mx)];              % State x4 from solution
x5 = [x0(5);z(5:mx:N*mx)];              % State x5 from solution
x6 = [x0(6);z(6:mx:N*mx)];              % State x6 from solution

num_variables = 5/delta_t;
zero_padding = zeros(num_variables,1);
unit_padding  = ones(num_variables,1);

u1  = [zero_padding; u1; zero_padding];
u2  = [zero_padding; u2; zero_padding];
x1  = [pi*unit_padding; x1; zero_padding] - pi;
x2  = [zero_padding; x2; zero_padding];
x3  = [zero_padding; x3; zero_padding];
x4  = [zero_padding; x4; zero_padding];
x5  = [zero_padding; x5; zero_padding];
x6  = [zero_padding; x6; zero_padding];
x = [x1 x2 x3 x4 x5 x6];

%% Plotting
t = 0:delta_t:delta_t*(length(u1)-1);

subplot(711)
stairs(t,u1),grid
title(['q1=' num2str(q1)]);
ylabel('u')
subplot(712)
plot(t,x1,'m',t,x1,'mo'),grid
ylabel('lambda')
subplot(713)
plot(t,x2,'m',t,x2','mo'),grid
ylabel('r')
subplot(714)
plot(t,x3,'m',t,x3,'mo'),grid
ylabel('p')
subplot(715)
plot(t,x4,'m',t,x4','mo'),grid
xlabel('tid (s)'),ylabel('pdot')
subplot(716)
plot(t,x5,'m',t,x5','mo'),grid
xlabel('tid (s)'),ylabel('e')
subplot(717)
plot(t,x6,'m',t,x6','mo'),grid
xlabel('tid (s)'),ylabel('e_dot')
figure;
stairs(t,u2),grid

seconds = (N+2*num_variables)*delta_t;
u = [linspace(0,seconds,size(u2,1))' u1 u2];
u1 = [linspace(0,seconds,size(u2,1))' u1];
u2 = [linspace(0,seconds,size(u2,1))' u2];
x = [linspace(0,seconds,size(x,1))' x];

compute_lq4;
