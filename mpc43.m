function main
clc
clear all

addpath('D:/casadi-windows-matlabR2016a-v3.4.5') 
import casadi.*

mpc=50;
T = 3; % Time horizon 3sT3V1B1P05K09XY13g09
N = 50; % number of control intervals
gamma = 0.5;
nu = 0.9;
b = 1;
rho = 0.6;
kappa = 0.5;
epslonn = 0;
y0=3;
x0=1;
z0 = [x0\y0];
n=1;
r=1;
% Declare model variables
t = SX.sym('t');
x1 = SX.sym('x1');

u = SX.sym('u');

% Model equations
xdot = [(x1 + gamma)*u - nu*x1];

% Objective term
f0 = -exp(-rho*t)*(kappa*log(x1)+log(b-u));%???^rho

% Formulate discrete time dynamics
% Fixed step Runge-Kutta 4 integrator
   
if false
% Создаем интегратор с помощью CVODES из SUNDIALS Suite
% использование: x(delta | 0, x0, u) = F('x0', x0, 'p', u)
    ode = struct('x', x, 'p', u, 'ode', xdot, 'quad', f);
    F = integrator('F', 'cvodes', ode, struct('tf', T/N));
else
   DT = T/N;
   f = Function('f', {x1, u, t}, {xdot, f0});
   X0 = MX.sym('X0', n);
   U = MX.sym('U');
   tt = MX.sym('tt');
   X = X0;
   Q = 0;
       [k1, k1_q] = f(X, U, tt);
       [k2, k2_q] = f(X + DT/2 * k1, U, tt+DT/2);
       [k3, k3_q] = f(X + DT/2 * k2, U,tt+DT/2);
       [k4, k4_q] = f(X + DT * k3, U,tt+DT);
       X=X+DT/6*(k1 +2*k2 +2*k3 +k4);
       Q = Q + DT/6*(k1_q + 2*k2_q + 2*k3_q + k4_q);
   
    F = Function('F', {X0, U, tt}, {X, Q}, {'x0','p','t'}, {'xf', 'qf'});
end   
% Собираем переменные и ограничения задачи
u = MX.sym('u',N);
z = MX.sym('z',n);
J = 0;
G = [];
x = z;
psis=0;%не знаю какое оно должно быть
for k = 1:N
    res = F('x0', x, 'p', u(k), 't', (k-1)*DT);    % x(delta | 0, x, u(k))
    x   = res.xf;
    J   = J + res.qf;
end
J=J-exp(-rho*T)*psis*x(mpc);
% Создаем NLP-solver
nlp = struct('f', J, 'x', [z(:); u(:)], 'g', []);
solver = nlpsol('solver', 'ipopt', nlp);

lbw = zeros(N,r);
ubw = (b-epslonn)*ones(N,r);
w0 = zeros(N,r);


%MPC

Nmpc = mpc;
xtau = z0;
X = z0;
U = [];
dv=[];
for tau = 0:Nmpc

% Solve the NLP
sol = solver('x0', [xtau; w0],  'lbg', [], 'ubg', [],...
             'lbx', [xtau; lbw], 'ubx', [xtau; ubw]);
xu = full(sol.x);
u_opt = reshape(xu(n+1:end), r, N);
    
 % выделить решение
 
  J_opt = full(sol.f);
  dv=full(sol.lam_x);
 % dvtau=full(sol.lam_x);
% найти следующее состояние
    res = F('x0', xtau, 'p', u_opt(1),'t', T/N);
    xtau = full(res.xf);
 % запомнить текущее состояние
    U = [U u_opt(1)];
    X = [X xtau];
    %dv=[dv dvtau(2)];
    % подготовить приближение
   w0 = [u_opt(2:end) zeros(1)]';
  % results(T/N*(0:tau+1),X,U,J_opt,t_Elapsed,1)
end
 
% Plot the solution
dv=[dv; 0];
x1_opt = X(1,:);

plot(x1_opt,dv)
% x1_opt = X(1,:);
% 
% U;
% tgrid = linspace(0, T, Nmpc+2);
% clf;
% hold on
% 
% plot(tgrid, x1_opt, '--', 'LineWidth',2)
% stairs(tgrid, [U'; nan], '-.','LineWidth',2)
% xlabel('t')
% legend('z','u')

end