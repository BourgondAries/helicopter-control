% TTK4135 - Helicopter lab
% Hints/template for problem 2.
% Changed for problem 4.

%% Initialization and model definition
init07; % NB: Change this to the init file corresponding to your helicopter
delta_t = 0.25; % sampling time
h = delta_t;
q1 = 1;
q2 = 1;

alpha1 = K_1*K_pp;
alpha2 = K_1*K_pd;
alpha3 = K_3*K_ep;
alpha4 = K_3*K_ed;

% Discrete time system model. x = [lambda r p p_dot e e_dot]'
A1 = [0   1   0       0        0         0;
      0   0   -K_2    0        0         0;
      0   0   0       1        0         0;
      0   0   -alpha1 -alpha2  0         0;
      0   0   0       0        0         1;
      0   0   0       0        -alpha3   -alpha4]*h + eye(6);
B1 = [0   0   0       alpha1   0         0;
      0   0   0       0        0         alpha3]'*h;

% Number of states and inputs
mx = size(A1,2); % Number of states (number of columns in A)
mu = size(B1,2); % Number of inputs (number of columns in B)

% Time horizon and initialization
N  = 15;                                % Time horizon for states
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
x0 = [x1_0 x2_0 x3_0 x4_0 x5_0 x6_0]'; % Initial values


% Bounds
ul      = -30*pi/180;                   % Lower bound on control -- u1
uu      = 30*pi/180;                    % Upper bound on control -- u1

xl      = -Inf*ones(mx,1);              % Lower bound on states (no bound)
xu      = Inf*ones(mx,1);               % Upper bound on states (no bound)
xl(3)   = ul;                           % Lower bound on state x3
xu(3)   = uu;                           % Upper bound on state x3

% Generate constraints on measurements and inputs
[vlb, vub]      = genbegr2(N,M,xl,xu,ul,uu);  % hint: genbegr2
vlb(N*mx+M*mu)  = 0;                          % We want the last input to be zero
vub(N*mx+M*mu)  = 0;                          % We want the last input to be zero

% Generate the matrix Q and the vector c (objective function weights in the QP problem)
Q1 = zeros(mx,mx);
Q1(1,1) = 1;                            % Weight on state x1
Q1(2,2) = 0;                            % Weight on state x2
Q1(3,3) = 0;                            % Weight on state x3
Q1(4,4) = 0;                            % Weight on state x4
Q1(5,5) = 0;                            % Weight on state x5
Q1(6,6) = 0;                            % Weight on state x6
P1 = diag([q1 q2]);                      % Weight on input
Q = genq2(Q1,P1,N,M,mu);                % Generate Q
c = zeros(N*mx+M*mu,1);                 % Generate c

fun = @(x) x'*Q*x;

%% Generate system matrixes for linear model
Aeq = gena2(A1,B1,N,mx,mu);           % Generate A, hint: gena2
beq = zeros(mx*N,1);        	      % Generate b
beq(1:mx) = A1*x0;                    % Initial value

%% Solve QP problem with linear model
tic;
z = fmincon(fun,z0,[],[],Aeq,beq,vlb,vub,@nonlcon);
t1=toc;

% Calculate objective value
phi1 = 0.0;
PhiOut = zeros(N*mx+M*mu,1);
for i=1:N*mx+M*mu
  phi1=phi1+Q(i,i)*z(i)*z(i);
  PhiOut(i) = phi1;
end

%% Extract control inputs and states
p_c  = [z(N*mx+1:mu:N*mx+M*mu); z(N*mx+M*mu)]; % Control input from solution
e_c  = [z(N*mx+2:mu:N*mx+M*mu); z(N*mx+M*mu)]; % Control input from solution

x1 = [x0(1);z(1:mx:N*mx)];              % State x1 from solution
x2 = [x0(2);z(2:mx:N*mx)];              % State x2 from solution
x3 = [x0(3);z(3:mx:N*mx)];              % State x3 from solution
x4 = [x0(4);z(4:mx:N*mx)];              % State x4 from solution
x5 = [x0(5);z(5:mx:N*mx)];              % State x4 from solution
x6 = [x0(6);z(6:mx:N*mx)];              % State x4 from solution

num_variables = 5/delta_t;
zero_padding = zeros(num_variables,1);
unit_padding  = ones(num_variables,1);

p_c   = [zero_padding; p_c; zero_padding];
e_c   = [zero_padding; e_c; zero_padding];
x1  = [pi*unit_padding; x1; zero_padding] - pi;
x2  = [zero_padding; x2; zero_padding];
x3  = [zero_padding; x3; zero_padding];
x4  = [zero_padding; x4; zero_padding];
x5  = [zero_padding; x5; zero_padding];
x6  = [zero_padding; x6; zero_padding];
x = [x1 x2 x3 x4 x5 x6];

%% Plotting
t = 0:delta_t:delta_t*(length(p_c)-1);

subplot(611)
stairs(t,p_c),grid
title(['q1=' num2str(q1)]);
ylabel('p_c')
subplot(612)
stairs(t,e_c),grid
title(['q2=' num2str(q2)]);
ylabel('e_c')
subplot(613)
plot(t,x1,'m',t,x1,'mo'),grid
ylabel('lambda')
subplot(614)
plot(t,x2,'m',t,x2','mo'),grid
ylabel('r')
subplot(615)
plot(t,x3,'m',t,x3,'mo'),grid
ylabel('p')
subplot(616)
plot(t,x4,'m',t,x4','mo'),grid
xlabel('tid (s)'),ylabel('pdot')

seconds = (N+2*num_variables)*delta_t;
p_c = [linspace(0,seconds,size(p_c,1))' p_c];
e_c = [linspace(0,seconds,size(e_c,1))' e_c];
x = [linspace(0,seconds,size(x,1))' x];