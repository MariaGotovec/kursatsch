function mpc41
clc
close all

odeopt = odeset('RelTol',1e-8,'AbsTol',1e-8); 

addpath('D:/casadi-windows-matlabR2016a-v3.4.5') 
%addpath('C:\Program Files\MATLAB\R2019b\casadi-windows-matlabR2016a-v3.5.1')
import casadi.*

% �����������
n = 1;
r = 1;
t0 = 0;
T = 10;  % ��������
N = 50; % ���-�� ���������� �����������
Nmpc = 200;

% ��������� ������ ��  CASE: 
gamma = 0.1;
nu = 0.5;
b = 0.55;
rho = 0.15;
kappa =1.5;
eps = 0;
pmax = 40;
zmax = 0.8;

% ��������� �������
y0=1;
x0=0.1;
z0 = x0/y0;

% ���������, ����� ������
[z_s, p_s] = SteadyState;

% ���������� � ��������� ������ 
t = SX.sym('t');
x1 = SX.sym('x1');
u = SX.sym('u');

% �������� ������� (������ ����� �.�.)
xdot = (x1 + gamma)*u - nu*x1;

% �������� ��������
f0 = -exp(-rho*t)*(kappa*log(x1)+log(b-u));%???^rho

% ���������� ��������
% �� 4 ������� � ������������� �����
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

% �������� ���������� � ����������� ������
u = MX.sym('u',N);
z = MX.sym('z',n);
J = 0;
G = [];
x = z;
for k = 1:N
    res = F('x0', x, 'p', u(k), 't', (k-1)*DT);    % x(delta | 0, x, u(k))
    x   = res.xf;
    J   = J + res.qf;
end

% ������� NLP-solver
nlp = struct('f', J, 'x', [z(:); u(:)], 'g', []);
solver = nlpsol('solver', 'ipopt', nlp);

lbw = zeros(N,r);
ubw = (b-eps)*ones(N,r);
w0 = zeros(N,r);


% % ������ ������ ����������� ����������������
% xtau = z0;
% tic
% sol = solver('x0', [xtau; w0],  'lbg', [], 'ubg', [],...
%              'lbx', [xtau; lbw], 'ubx', [xtau; ubw]);
% t_Elapsed = toc;         
% % ��������������� �������: ���������� � ������������ ����������         
% J_opt = full(sol.f);
% primal = full(sol.x);
% 
% u_opt = reshape(primal(n+1:end), r, N);
% [x_opt, p_opt] = PIsystem(0:DT:T, xtau, 0, u_opt);
%     
%results(0:DT:T, x_opt, u_opt,p_opt, J_opt, t_Elapsed, 1)

% ������  mpc
xtau = z0;

z_opt = z0;
w0 =  zeros(N,r);
Uz = [];
XX = nan*zeros(Nmpc,Nmpc+N+2);
p_opt=[];
ptau=[nan];
tic
Tmpc=DT*Nmpc
for tau = t0:DT:Tmpc

% Solve the NLP
    
    sol = solver('x0', [xtau; w0],  'lbg', [], 'ubg', [],...
             'lbx', [xtau; lbw], 'ubx', [xtau; ubw]);
    t_Elapsed = toc;   
    xu = full(sol.x);
    u_opt = reshape(xu(n+1:end), r, N);
    % �������� �������
    J_opt = full(sol.f);
    
        if tau>0   && tau <Tmpc 
            [XX(fix(tau/DT),fix(tau/DT):N+fix(tau/DT)), ptau] = PIsystem(0:DT:T, xtau, 0, u_opt);   
        end
    % ����� ��������� ���������
    res = F('x0', xtau, 'p', u_opt(1),'t', T/N);
    xtau = full(res.xf);
    % ��������� ������� ���������
    Uz = [Uz u_opt(1)];
     z_opt = [z_opt xtau];
     p_opt=[p_opt ptau(1)];
    % ����������� �����������
    w0 = [u_opt(2:end) zeros(1)]';
end

u_opt = Uz;
x_opt = z_opt; 
XX(1,1:Nmpc+2)=x_opt;

%p_opt = mpcPIsystem(0:T/Nmpc:T,0,x_opt,u_opt);

% resultF(0:T/Nmpc:T, x_opt, u_opt,p_opt, J_opt, t_Elapsed, 1)
% 
% [z_s, p_s] = SteadyState;
results2(0:(T+Tmpc)/(N+Nmpc):T+Tmpc, x_opt,u_opt, XX , J_opt, t_Elapsed, 1,p_opt)
%-------------------------------------------------------------
% FUNCTIONS
%-------------------------------------------------------------
function h = h(z)
    h = 1/( b*(z+gamma) );
end
function phi = phi_1(z,p)
    if p>=h(z)
        phi = (b-nu)*z + b*gamma -1/p;
    else
        phi = -nu*z;
    end
end
function phi = phi_2(z,p)
    if p>=h(z)
        phi = -(b-nu-rho)*p -( gamma*kappa + (kappa-1)*z)/( (z+gamma)*z ) ;
    else
        phi = (nu+rho)*p - kappa/z;
    end
end

%-------------------------------------------------------------
function [z_s,p_s] = SteadyState
    % �������� ������� 
    beta0 = rho + nu - b;
    beta1 = rho + nu - kappa*b;
    beta2 = rho + kappa*nu - kappa*b;
    
    fprintf('beta0 = %+8.4f; beta1 = %+8.4f; beta2 = %+8.4f;  CASE: ', ...
        beta0, beta1, beta2);

    if (beta0>0)  && (beta1>0) && (beta2>0) && (nu<b)
        fprintf('A1B1C1\n');  
    end
    if (beta0>0)  && (beta1<=0) && (beta2>0) && (nu<b)
        fprintf('A1B1C2\n');  
    end
    if (beta0>0)  && (beta1>0) && (beta2>0) && (nu>=b)
        fprintf('A1B1C3\n');  
    end
    if (beta0>0)  && (beta1<=0) && (beta2>0) && (nu>=b)
        fprintf('A1B1C4\n');  
    end
    if (beta0<0)  && (beta1>0) && (beta2>0)
        fprintf('A1B1C7\n');  
    end
    sol = fsolve(@(x) [phi_1(x(1),x(2)); phi_2(x(1),x(2))], [1;1]);
    z_s = sol(1); p_s = sol(2);
    fprintf('z_s = %+8.4f; p_s = %+8.4f \n', z_s, p_s);
    
    % ������� �������
%     fig = figure(100);
%     set(fig,'Position', [200    40   400   400]);
%     hold on; grid on; 
%     warning('off')
%     H = [];
%     for zi = 0:0.01:zmax
%         H = [H h(zi)];
%     end
%     plot(0:0.01:zmax,H,'k'); % ������ ������� h(z)
%    
%     fill([0 0:0.01:zmax zmax], [0 H 0],[0.9 0.9 0.9]);
%         
%     fimplicit( @phi_1,[0.01 zmax 0.01 pmax],'--r'); % ������ ������� V1
%     fimplicit( @phi_2,[0.01 zmax 0.01 pmax],'--r'); % ������ ������� V1
%     warning('on')
end                
%-------------------------------------------------------------
function [X,P] = PIsystem(T,x0,p0,ud)
    ud = [ud ud(end)];
    u = @(t) ud(floor(t/DT)+1);
    sys_z = @(t,z) u(t)*(z+gamma) - nu*z;
    sol_z = ode45(sys_z, T, x0, odeopt);
    X = deval(sol_z,T);
    
    sys_p = @(t,p) (rho+nu-u(t))*p - kappa/deval(sol_z,t);
    sol_p = ode45(sys_p, flip(T), 0, odeopt);
    P = deval(sol_p,T);

end

function P = mpcPIsystem(T,p0,xd,ud)
    ud = [ud ud(end)];
    u = @(t) ud(floor(t/DT)+1);
    xd = [xd xd(end)];
    xm = @(t) xd(floor(t/DT)+1);
    sys_p = @(t,p) (rho+nu-u(t))*p - kappa/xm(t);
    sol_p = ode45(sys_p, flip(T), 0, odeopt);
    P = deval(sol_p,T);

end

%-------------------------------------------------------------
function results(T, X, U, Psi, J_opt, time, fNum)
    fprintf('   k  |      u        x        psi      \n');
    fprintf('----------------------------------------\n');
    for k = 1:length(T)-1
        fprintf(' %3.1f  | %+11.6f %+11.6f %+11.6f \n', T(k), U(k), X(k), Psi(k));
    end
    fprintf('������.���������  = %+11.6f \n\n', X(1,end));
    fprintf('�����. ��������   = %6.4f \n\n', J_opt);
    fprintf('����� �������     = %6.4f \n\n', time);

    fig = figure(fNum);
    set(fig,'Position', [0    40   800   600]);

    subplot(2,2,1); hold on; grid on;
    stairs(T,[U nan], '--', 'Linewidth', 1);
    ylim([0,b-eps])
    xlabel('t')
    ylabel('u')
    title('����������� ����������')
    
    subplot(2,2,2); hold on; grid on; 
    plot(T, X, '-b', 'Linewidth', 1)
    xlabel('t')
    ylabel('z')
    title('����������� ����������')
    
    subplot(2,2,3); hold on; grid on; 
    plot(T, Psi, '-b', 'Linewidth', 1)
    xlabel('t')
    ylabel('p')
    title('����������� ������������')
    
    subplot(2,2,4); hold on; grid on; 
    % ��������������� ����������
    warning('off')
    H = [];
    for zi = 0:0.01:zmax
        H = [H h(zi)];
    end 
    plot(0:0.01:zmax,H,'k'); % ������ ������� h(z)
    fill([0 0:0.01:zmax zmax], [0 H 0],[0.9 0.9 0.9]);
        
    fimplicit( @phi_1,[0.01 zmax 0.01 pmax],'--r'); % ������ ������� V1
    fimplicit( @phi_2,[0.01 zmax 0.01 pmax],'--r'); % ������ ������� V1
    warning('on')
    % ������� �������
    plot(X, Psi, '-b', 'Linewidth', 1)
    xlabel('z')
    ylabel('p')
    title('������� �������')
    
end

function results2(T, X, U,xx, J_opt, time, fNum,Psi)
%    fprintf('   k  |      u        x       \n');
%    fprintf('----------------------------------------\n');
%     for k = 1:length(T)-1-Nmpc
%         fprintf(' %3.1f  | %+11.6f %+11.6f \n', T(k), U(k), X(k));
%     end
    fprintf('������.���������  = %+11.6f \n\n', X(1,end));
    fprintf('�����. ��������   = %6.4f \n\n', J_opt);
    fprintf('����� �������     = %6.4f \n\n', time);

    fig = figure(fNum);
    set(fig,'Position', [0    40   1200   300]);

    subplot(1,3,1); hold on; grid on; 
    plot([T nan], xx(1,:), '-b', 'Linewidth', 1)
    for k = 2:4:Nmpc
        plot([T nan], xx(k,:), '--b', 'Linewidth', 1)
    end
    xlabel('t')
    ylabel('z')
    title('a)')
    
    subplot(1,3,2); hold on; grid on;
    stairs(T,[U nan*zeros(1,N)], '--b', 'Linewidth', 1);
    ylim([0,b-eps])
    xlabel('t')
    ylabel('u')
    title('�)')
    
    subplot(1,3,3); hold on; grid on; 
        warning('off')
    H = [];
    for zi = 0:0.01:zmax
        H = [H h(zi)];
    end 
    plot(0:0.01:zmax,H,'k'); % ������ ������� h(z)
    fill([0 0:0.01:zmax zmax], [0 H 0],[0.9 0.9 0.9]);
        
    fimplicit( @phi_1,[0.01 zmax 0.01 pmax],'--r'); % ������ ������� V1
    fimplicit( @phi_2,[0.01 zmax 0.01 pmax],'--r'); % ������ ������� V1
    warning('on')
    % ������� �������
    plot(X, [Psi nan], '-b', 'Linewidth', 1)
    xlabel('z')
    title('�)')

end


end