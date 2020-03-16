clc
clear all
addpath('D:/casadi-windows-matlabR2016a-v3.4.5') 
import casadi.*

T = 10; % Time horizon 1sT3V06B1P05K09XY13g09
N = 50; % number of control intervals
gamma = 0.6;
nu = 0.7;
b = 1;
rho = 0.3;
kappa = 0.3;
epslonn=0;
z00 = [1; 3];
z0 = [z00(1)/z00(2)];

% Declare model variables
z1 = SX.sym('z1');
t = SX.sym('t');
u = SX.sym('u');

% Model equations
zdot = [(z1 + gamma)*u - nu*z1];

% Objective term
f0 = -exp(-rho*t)*(kappa*log(z1)+log(b-u));

% Formulate discrete time dynamics
% Fixed step Runge-Kutta 4 integrator
   M = 4; 
   DT = T/N/M;
   f = Function('f', {z1, u, t}, {zdot, f0});
   X0 = MX.sym('X0', 1);
   U = MX.sym('U');
   tt = MX.sym('tt');
   X = X0;
   Q = 0;
   for j=1:M
       [k1, k1_q] = f(X, U, tt);
       [k2, k2_q] = f(X + DT/2 * k1, U, tt+DT/2);
       [k3, k3_q] = f(X + DT/2 * k2, U,tt+DT/2);
       [k4, k4_q] = f(X + DT * k3, U,tt+DT);
       X=X+DT/6*(k1 +2*k2 +2*k3 +k4);
       Q = Q + DT/6*(k1_q + 2*k2_q + 2*k3_q + k4_q);
    end
    F = Function('F', {X0, U, tt}, {X, Q}, {'x0','p','t'}, {'xf', 'qf'});

  
% Start with an empty NLP
w={};
w0 = [];
lbw = [];
ubw = [];
J = 0;
g={};
lbg = [];
ubg = [];

% Formulate the NLP

Xk = z0;

for k=0:N-1
    % New NLP variable for the control
    Uk = MX.sym(['U_' num2str(k)]);
    w = {w{:}, Uk};
    lbw = [lbw, 0];
    ubw = [ubw,  b-epslonn];
    w0 = [w0,  0];

    % Integrate till the end of the interval
    Fk = F('x0',Xk,'p', Uk);
    Xk = Fk.xf;
    J = J+Fk.qf;

    % Add inequality constraint
    g = {g{:}, Xk(1)};%
    lbg = [lbg; 0];
    ubg = [ubg;  inf];
end

% Create an NLP solver
prob = struct('f', J, 'x', vertcat(w{:}), 'g', vertcat(g{:}));
solver = nlpsol('solver', 'ipopt', prob);

% Solve the NLP
sol = solver('x0', w0, 'lbg', lbg, 'ubg', ubg, ...
             'lbx', lbw, 'ubx', ubw);
w_opt = full(sol.x);
dv=full(sol.lam_x);%двойственная
dv=[dv;0];
% Plot the solution
u_opt = w_opt;
x_opt = z0;%
for k=0:N-1
   Fk = F('x0', x_opt(:,end), 'p', u_opt(k+1));
   x_opt = [x_opt, full(Fk.xf)];
end
z1_opt = x_opt(1,:);
plot(z1_opt,dv)
% tgrid = linspace(0, T, N+1);
% clf;
% hold on
% plot(tgrid, z1_opt, '--', 'LineWidth',2)
% 
% stairs(tgrid, [u_opt; nan], '-.','LineWidth',2)
% xlabel('t')
% legend('z','u')
