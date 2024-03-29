%% DMP_discrete: following the studywolf implementation
% Natasha Rouse (April 2021)
% Including more functions from original Schaal code
% Init & Batch_Fit based on Schaal "learn_dcp_batch"
close all
clear all
clc
tic

%% Set parameters
dt = 0.0001; %time increment
bfs = 10;
n_dmps = 1;
tau = 1;
% timesteps = 100;
% n_steps = 2*tau/dt+2;
n_steps = 30000;

%DESIRED TRAJECTORY
l = 1;

y = sin(2.8142*[0:dt:dt*(n_steps-1)]);
yd = [0 diff(y)]/dt;
yd(1) = yd(2);
ydd = [0 diff(yd)]/dt;
ydd(1) = ydd(2);
Y  =[y;yd;ydd];

% figure(4)
% line(desired_traj(:,1),desired_traj(:,2))

% D = dmp_init(n_dmps,bfs,dt,desired_traj(1,:),desired_traj(end,:));

%CENTERS FOR BFS
%From Schaal
% t = (0:1/(D.bfs-1):1)*0.5;
% D.c = exp(-D.cs.ax(1)*t);
% D.cd = (-D.cs.ax(1))*D.c;
%From studywolf
% c = linspace(0,D.cs.runtime,D.bfs);
% for i = 1:length(c)
%     D.c(i) = exp(-D.cs.ax*c(i));
% end
% %From SMP paper
% D.c = logspace(log10(1), log10(.01), bfs);

%VARIANCE OF BFS (TRIAL AND ERROR)
%From studywolf
% D.h = ones(1,D.bfs)*D.bfs^1.5 ./ (D.c'./D.cs.ax)';
%From Schaal
% D.h = (diff(D.c*0.5)).^2;
% D.h = 1./[D.h,D.h(end)];
% %From SMP paper
% D.h = linspace(log10(.3),log10(.002), bfs);

% %Check Offset
% if abs(D.y0-D.goal)<1e-4
%     D.goal = D.goal + 1e-4;
% end
% [D,y_des] = imitate_path(D,desired_traj);
% [D,y,yd,ydd] = dmp_rollout(D,tau);

%% Schaal "learn_dcp_batch" process
%Variables for plotting
% Z = zeros(2,floor(2*tau/D.dt+1),D.n_dmps);
% X = Z;
% V=Z;
% T = zeros(3,floor(2*tau/D.dt+1),D.n_dmps);
% Y=T;
% PSI = zeros(D.bfs,(2*tau/D.dt+1),D.n_dmps);
% W = zeros(D.bfs,(2*tau/D.dt+1),D.n_dmps);
% 
% %generate the minimum jerk trajectory
% t = 0; % time
% td = 0; % time step
% tdd = 0; % time acceleration
% for i = 0:2*tau/dt
%     [t,td,tdd] = min_jerk_step(t,td,tdd,desired_traj(i+2,:),tau-i*dt,dt);
%     T(:,i+1,:) = [t; td; tdd]; % store time, rate of change of time, and temporal acceleration for each degree of freedom in operation space
% end
% 
% %use batch_fit to initialize with minjerk
% [D,~,~,~] = dmp_batch_fit(D,tau,T(1,:,:),T(2,:,:),T(3,:,:));
% 
% %test the fit
% % D = dmp_reset(D);
% D = dmp_set_goal(D,desired_traj(end,:),1);
% 
% for i = 0:2*tau/dt
%     [D,y,yd,ydd,weight] = dmp_run(D,tau); % y is the actual trajectory of the end effector
%     
%     Z(:,i+1,:) = [D.z; D.zd]; 
%     Y(:,i+1,:) = [y; yd; ydd];
%     X(:,i+1,:) = [D.x; D.xd];
%     V(:,i+1,:) = [D.v; D.vd];
%     PSI(:,i+1,:) = D.psi';
% %     W(:,i+1,:) = D.w';
%     W(:,i+1,:) = weight;
% end
% 
% %PLOTTING
% time = (0:dt:tau*2)';
% 
% figure(5)
% clf(5)
% figure(5)
% 
% for i = 1:D.n_dmps
%     hold on
%     
% %position, velocity, acceleration vs. target
% subplot(4,3,1)
% plot(time,[Y(1,:,i)' T(1,:,i)']);
% title('y')
% 
% subplot(4,3,2)
% plot(time,[Y(2,:,i)' T(2,:,i)']);
% title('yd')
% 
% subplot(4,3,3)
% plot(time,[Y(3,:,i)' T(3,:,i)']);
% title('ydd')
% 
% %internal states
% subplot(4,3,4)
% plot(time,Z(1,:,i))
% title('z')
% 
% subplot(4,3,5)
% plot(time,Z(2,:,i))
% title('zd')
% 
% subplot(4,3,6)
% plot(time,PSI(:,:,i))
% title('{Psi activations/},{Basis functions}')
% 
% subplot(4,3,7)
% plot(time,V(1,:,i))
% title('v')
% 
% subplot(4,3,8)
% plot(time,V(2,:,i))
% title('vd')
% 
% subplot(4,3,9)
% plot(time,W(:,:,i))
% title('Linear Model Weights over Time')
% 
% subplot(4,3,10)
% plot(time,X(1,:,i))
% title('x')
% 
% subplot(4,3,11)
% plot(time,X(2,:,i))
% title('xd')
% 
% subplot(4,3,12)
% plot(W(end,:,i))
% title('Weights')
% xlabel(sprintf('tau=%f',tau))
% 
% drawnow
% end
% 
% figure(6)
% hold on
% traj_time = 0:dt:dt*(size(desired_traj, 1)-1);
% subplot(2,1,1)
% plot(traj_time, desired_traj(:,1),time, Y(1,:,1))
% title('Desired vs DMP Y Trajectory')
% % subplot(2,1,2)
% % plot(traj_time, desired_traj(:,2), time,Y(1,:,2))
% % title('Desired vs DMP Z Trajectory')
% 
% % subplot(3,1,3)
% % plot(traj_time, desired_traj(:,3),time,Y(1,:,3))
% % title('Desired vs. Produced Z Trajectory')
% % legend('Desired Z', 'Produced Z')
% toc
%% Dynamics
m_s = 1; % sampler mass
m_c = 1; % cart mass
ndofs = 2; 
g = 9.81;

% initial conditions
q = nan(ndofs, size(Y,2));
q(:,1) = zeros(ndofs, 1);

q_dot = nan(ndofs, size(Y,2));
q_dot(:,1) = [Y(2,1),0]';

% q_ddot = nan(ndofs, size(Y,2));
% q_ddot(:,1) = zeros(ndofs, 1);

u_hist = nan(ndofs, size(Y,2));

figure(7)
%     v = VideoWriter('Fixed Length Underactuated without DMP');
%     v.FrameRate = 10;
%     open(v)
for i = 1:size(Y, 2)-1
    y_ddot =  permute(Y(3,i,:), [3,2,1]);
    q_out = rk4(dt,@dynamics, (i-1)*dt, [q(:,i); q_dot(:,i)], y_ddot, m_c, m_s, l, g);
    
    q(:,i+1) = q_out(1:ndofs);
    q_dot(:,i+1) = q_out(ndofs+1:end);
    
    if mod(i,100)==0
       figure(7);
       plot([q(1,i+1),q(1,i+1)+l*sin(q(2,i+1))],[0,-l*cos(q(2,i+1))],'-b',...
           q(1,i+1)+l*sin(q(2,i+1)),-l*cos(q(2,i+1)),'.k','MarkerSize',20);
       yline(0);
       rectangle('Position',[q(1,i+1)-(.25/2) -(.1/2) 0.25 0.1]);
       axis([-3 3 -5.5 .5]);
       drawnow

%        frame = getframe(gcf);
%        writeVideo(v, frame);
   end


end
%    close(v)

%%
figure(6);
% subplot(2,1,1)
hold on
y0 = q(1,1)+l.*sin(q(2,1));
time = 0:dt:dt*(n_steps-1);
plot(time, Y(1,:), time, q(1,:)+l.*sin(q(2,:))-y0)
legend('Desired Y',  'End Effector Y')

hold off
figure(8);
plot(time, q(1,:),time, q(2,:))
legend('y_c', 'theta')
%%
function [dydt, U] = dynamics(t, Y, y_ddot, m_c, m_s, l, g)
    
    y_ddot_c = (l*sin(Y(2))*Y(4)^2+y_ddot+g*cos(Y(2))*sin(Y(2)))/(0.0001+sin(Y(2))^2);
    theta_ddot = -(y_ddot_c*cos(Y(2))+g*sin(Y(2)))/l;
    
    F_yc = (m_s+m_c)*y_ddot_c+m_s*(l*Y(4)*cos(Y(2))-l*Y(4)^2*sin(Y(2)));
    M_theta = 0;
    
    U = [F_yc; M_theta];
    
    ndofs = length(Y)/2;
    dydt = [Y(ndofs+1:end); y_ddot_c; theta_ddot];
end
%{
function [u, q_ddot, q_dot, q] = robot_dynamics(y_ddot, q, q_dot,m_s, m_c, dt)

    y_c = q(1);
    l = 1;
    theta = q(2);
    
    y_dot_c = q_dot(1);
    l_dot = 0;
    theta_dot = q_dot(2);

    g = 9.81;
                        
    y_ddot_c = (l*sin(theta)*theta_dot^2+y_ddot+g*cos(theta)*sin(theta))/(0.0001+sin(theta)^2);
    theta_ddot = -(y_ddot_c*cos(theta)+g*sin(theta))/l;
    F_yc = (m_s+m_c)*y_ddot_c+m_s*(l*theta_ddot*cos(theta)-l*theta_dot^2*sin(theta));
        

    q_ddot = [y_ddot_c; theta_ddot];
    
    Q_old = [y_c; y_dot_c; theta; theta_dot];
    
    M_theta = 0;
%     Q_new = rk4(dt, 
    Q_new = rk4(Q_old, dt, F_yc, M_theta, m_s, m_c, g, l);
    
    
    u = [F_yc; 0]; % control forces
    
    q = [Q_new(1); Q_new(3)];
    q_dot = [Q_new(2); Q_new(4)];
end
%}
%%

function qqdot_new = rk4(h, fn, t, qqdot, y_ddot, m_c, m_s, l, g)
    % qqdot = [q; q_dot], where q and q_dot are vertical vectors
    % representing the state variables and their derivatives, respectively
    k1 = h*fn(t, qqdot, y_ddot, m_c, m_s, l, g);
    k2 = h*fn(t + h/2, qqdot + k1/2, y_ddot, m_c, m_s, l, g);
    k3 = h*fn(t + h/2, qqdot + k2/2, y_ddot, m_c, m_s, l, g);
    k4 = h*fn(t + h, qqdot + k3, y_ddot, m_c, m_s, l, g);
    
    qqdot_new = qqdot + k1/6 + k2/3 + k3/3 + k4/6;
    
end
%% DMP initialization
function dmp = dmp_init(n_dmps,bfs,dt,y0,goal,w,ay,by)

dmp.n_dmps = n_dmps;
dmp.bfs = bfs;
dmp.dt = dt;
dmp.y0 = ones(1,n_dmps).*y0;

%state variables
dmp.z       = zeros(1,dmp.n_dmps);
dmp.y       = zeros(1,dmp.n_dmps); 
dmp.x       = zeros(1,dmp.n_dmps);
dmp.v       = zeros(1,dmp.n_dmps);
dmp.zd      = zeros(1,dmp.n_dmps);
dmp.yd      = zeros(1,dmp.n_dmps); 
dmp.xd      = zeros(1,dmp.n_dmps);
dmp.vd      = zeros(1,dmp.n_dmps);
dmp.ydd     = zeros(1,dmp.n_dmps);

dmp.g = 0;
dmp.gd = 0;
dmp.G = 0;

if nargin<6
    dmp.w = zeros(n_dmps,bfs);
else
    dmp.w = w;
end
if nargin<7
%     dmp.ay = ones(1,n_dmps)*25;
    dmp.ay = 25;
else
    dmp.ay = ay;
end
if nargin<8
    dmp.by = dmp.ay / 4;
else
    dmp.by = by;
end
dmp.ag = dmp.ay/2;

dmp.A = 0; %original amplitude (max(y)-min(y)) when the primitive was fit
dmp.s = 1; %scale factor for the nonlinear function
dmp.dG = 0; %original goal amplitude (G-y0) when the primitve was fit
dmp.sx2 = zeros(1,bfs);
dmp.sxtd = zeros(1,bfs);
dmp.lambda = 1;
%Canonical state variables
dmp.cs.ax = dmp.ay/3;
dmp.cs.runtime = 1;
dmp.cs.dt = dt;
dmp.cs.timesteps = dmp.cs.runtime / dmp.cs.dt;
% dmp.cs.x = 1;

dmp = dmp_reset(dmp);
end

%% DMP goal
function dmp = dmp_set_goal(dmp,goal,flag)
dmp.G = goal;
dmp.g = goal;
%flag: update x0 with current state
if (flag)
    dmp.x = ones(1,dmp.n_dmps);
    dmp.y0 = dmp.y;
end
if dmp.A ~= 0
    if dmp.A/(abs(dmp.dG)+1e-10)>2.0
        'scaling needs to be set explicitly'
    else
        dmp.s = (dmp.G-dmp.y0)/dmp.dG;
    end
end
end

%% DMP reset
function dmp = dmp_reset(dmp)
dmp.y = dmp.y0;
dmp.yd = zeros(1,dmp.n_dmps);
dmp.ydd = zeros(1,dmp.n_dmps);
dmp.x = ones(1,dmp.n_dmps);
end

%% DMP run
function [dmp,y,yd,ydd,weights] = dmp_run(dmp,tau)
%weighted sum of locally weighted regression models
dmp.psi = exp(-0.5*((dmp.x'-dmp.c).^2).*dmp.h);
amp = dmp.s;
in = dmp.x;
f = sum(in*dmp.w.*dmp.psi)/sum(dmp.psi+1e-10)*amp;

dmp.vd = zeros(1,dmp.n_dmps);
dmp.xd = dmp.cs.ax*(0-dmp.x)*tau;
dmp.zd = (dmp.ay.*(dmp.by.*(dmp.g-dmp.y)-dmp.z)+f)*tau;
dmp.yd = dmp.z*tau;
dmp.ydd = dmp.zd*tau;

dmp.gd = dmp.ag.*(dmp.G-dmp.g);

dmp.x = dmp.xd*dmp.dt + dmp.x;
dmp.v = dmp.vd*dmp.dt + dmp.v;

dmp.z = dmp.zd*dmp.dt + dmp.z;
dmp.y = dmp.yd*dmp.dt + dmp.y;

dmp.g = dmp.gd*dmp.dt + dmp.g;

y = dmp.y;
yd = dmp.yd;
ydd = dmp.ydd;
weights = in*dmp.psi/sum(dmp.psi+1e-10)*amp;
end

%% DMP rollout
function [dmp,y_track,dy_track,ddy_track] = dmp_rollout(dmp,tau)
dmp = dmp_reset(dmp);

if (tau)
    dmp.cs.timesteps = dmp.cs.timesteps/tau;
end

%set up tracking vectors
y_track = zeros(dmp.cs.timesteps,dmp.n_dmps);
dy_track = zeros(dmp.cs.timesteps,dmp.n_dmps);
ddy_track = zeros(dmp.cs.timesteps,dmp.n_dmps);


for t = 1:dmp.cs.timesteps
    %run and record timesteps
    [dmp,y_track(t,:),dy_track(t,:),ddy_track(t,:)] = dmp_step(dmp,0,0);
end
end

%% DMP step
function [dmp,y,yd,ydd] = dmp_step(dmp,error,external_force)
%run the DMP system for a single timestep
tau = 1; %scales the timestep. Increase tau to make the system execute faster
error_coupling = 1/(1+error); %error is optional system feedback

%run canonical system
% x = dmp.cs.x + (-dmp.cs.ax*dmp.cs.x*error_coupling)*tau*dmp.cs.dt;
% dmp.cs.x = x;
[dmp,x] = cs_step(dmp,tau,error_coupling);

%generate basis function activation
psi = exp(-dmp.h .* (x-dmp.c).^2);
% dmp.w
for i = 1:dmp.n_dmps
    %generate forcing term
    f = (x*(dmp.goal(i)-dmp.y0(i)))*dot(psi,dmp.w(i,:));
    if abs(sum(psi)) > 1e-6
        f = f/sum(psi);
    end
    %dmp acceleration
    dmp.ydd(i) = dmp.ay(i)*(dmp.by(i)*dmp.goal(i) - dmp.y(i) - dmp.yd(i)) + f;
    if (external_force)
        dmp.ydd(i) = dmp.ydd(i) + external_force;
    end
    dmp.yd(i) = dmp.yd(i) + dmp.ydd(i)*tau*dmp.dt*error_coupling;
    dmp.y(i) = dmp.y(i) + dmp.yd(i)*tau*dmp.dt*error_coupling;
end
y = dmp.y;
yd = dmp.yd;
ydd = dmp.ydd;
end

%% Canonical System step
function [dmp,x_state] = cs_step(dmp,tau,error_coupling)
% generate a single step of x for discrete (potentially closed) loop
% movements. Decaying from 1 to 0 according to dx = -ax*x.

%tau: gain on execution time. increase to make system execute faster.
%error_coupling: slow down if the error is > 1
% (-dmp.cs.ax*dmp.cs.x*error_coupling)*tau*dmp.cs.dt
x_state = dmp.cs.x + (-dmp.cs.ax*dmp.cs.x*error_coupling)*tau*dmp.cs.dt;
dmp.cs.x = x_state;
end

%% CS rollout
function [dmp,x_track] = cs_rollout(dmp,tau,error_coupling)
if nargin == 1
    tau = 1;
    error_coupling = 1;
end
x_track = zeros(dmp.cs.timesteps,1);
%reset CS
dmp.cs.x = 1;
for t = 1:dmp.cs.timesteps
    [dmp,dmp.cs.x] = cs_step(dmp,tau,error_coupling);
    x_track(t) = dmp.cs.x;
end
end

%% Imitate path
function [dmp,y_des] = imitate_path(dmp,y_des)
%take in a desired trajectory and generate the set of system parameters to
%best realize this path

%set initial state and goal
dmp.y0 = y_des(1,:);
dmp.y_des = y_des;
dmp.goal = y_des(end,:);
if abs(dmp.y0-dmp.goal)<1e-4
    dmp.goal = dmp.goal + 1e-4;
end

path = zeros(dmp.n_dmps,dmp.cs.timesteps);
x = linspace(0,dmp.cs.runtime,length(y_des));
for i = 1:dmp.n_dmps
    for t = 1:dmp.cs.timesteps
        path(i,t) = interp1(x,y_des(:,i),t*dmp.cs.dt);
    end
end
y_des = path';

%calculate velocity of y_des with central difference
dy_des = (gradient(y_des))/dmp.cs.dt;
%calculate acceleration of y_des
ddy_des = gradient(dy_des)/dmp.cs.dt;

f_target = zeros(length(y_des),dmp.n_dmps);
%find the force required to move along this trajectory
for k = 1:dmp.n_dmps
    f_target(:,k) = ddy_des(:,k) - dmp.ay(k)*(dmp.by(k)*(dmp.goal(k)-y_des(:,k))...
        -dy_des(:,k));
end

%efficiently generate weights to realize f_target
dmp = generate_weights(dmp,f_target);

[dmp,x_track] = cs_rollout(dmp);
psi_track = exp(-dmp.h .* (x_track-dmp.c).^2);

%Show Plot
figure(1)
 clf(1)
 figure(1)
%Basis function activations
subplot(2,2,1)
title('Basis Functions')
% plot(psi_track)
hold on
for num = 1:size(psi_track,2)
    plot(psi_track(:,num))
end

%Weighted basis functions
subplot(2,2,2)
title('Weighted Basis Functions')
hold on
for bf = 1:size(psi_track,2)
    plot(psi_track(:,bf)*dmp.w(1,bf),'-')
    plot(psi_track(:,bf)*dmp.w(2,bf),'.')
end

%Plot desired forcing function vs approx
subplot(2,2,3)
title('Desired DMP forcing function')
hold on
for n = 1:dmp.n_dmps
    plot(f_target(:,n),'--')
end
% legend('f_{target}`','w*psi')

%approx forcing function
subplot(2,2,4)
title('Approx DMP forcing function')
hold on
for m = 1:dmp.n_dmps
    plot((psi_track*dmp.w(m,:)')*dmp.cs.dt)
end

psi_plot = psi_track*dmp.w'*dmp.cs.dt;

%Spatial forcing function
figure(2)
 clf(2)
 figure(2)
subplot(1,2,1)
title('Spatial desired forcing function')
plot(f_target(:,1),f_target(:,2),'--')

subplot(1,2,2)
title('Spatial approx. forcing function')
plot(psi_plot(:,1),psi_plot(:,2),'-')
end

%% Generate Weights
function dmp = generate_weights(dmp,f_target)

%calculate x & psi
[dmp,x_track] = cs_rollout(dmp);
psi_track = exp(-dmp.h .* (x_track-dmp.c).^2);
%efficiently calculate bf weights using weighted linear regression
dmp.w = zeros(dmp.n_dmps,dmp.bfs);

for i = 1:dmp.n_dmps
    %spatial scaling term
    k = dmp.goal(i) - dmp.y0(i);
    for j = 1:dmp.bfs
        numer = sum(x_track.*psi_track(:,j).*f_target(:,i));
        denom = sum(x_track.^2 .* psi_track(:,j));
        dmp.w(i,j) = numer/denom;
        if abs(k)>1e-5
            dmp.w(i,j) = dmp.w(i,j)/k;
        end
    end
end
dmp.w(isnan(dmp.w))=0;
end

%% DMP run & fit
function [y,yd,ydd] = dmp_run_fit(dmp,tau,t,td,tdd)

%first check whether this is the first time the primitive is fit, and
%record amplitude & dG info
if dmp.A == 0
    dmp.dG = dmp.G - dmp.y0;
    if dmp.cs.x == 1
        min_y = 1e10;
        max_y = -1e10;
        dmp.s = 1;
    end
end

%the regression target
amp = dmp.s;
ft = (tdd/tau^2 - dmp.ay*(dmp.by*(dmp.g-t)-td/tau))/amp;

%the weighted sum of the locally weighted regression models
dmp.psi = exp(-0.5*((dmp.x-dmpc).^2).*dmp.h);

%update the regression
dmp.sx2 = dmp.sx2*dmp.lambda + dmp.psi*dmp.x^2;
dmp.sxtd = dmp.sxtd*dmp.lambda + dmp.psi*dmp.x*ft;

%compute nonlinearity
in = dmp.x;
f = sum(in*dmp.w*dmp.psi)/sum(dmp.psi+1e-10) * amp;

%integrate
dmp.vd = 0;
dmp.xd = dmp.ax*(0-dmp.cs.x)*tau;

dmp.zd = (dmp.zy*(dmp.by*(dmp.g-dmp.y)-dmp.z)+f)*tau;
dmp.yd = dmp.z*tau;
dmp.ydd = dmp.zd*tau;
dmp.gd = dmp.ag*(dmp.G-dmp.g);
dmp.x = dmp.xd*dt + dmp.x;
dmp.v = dmp.vd*dt + dmp.v;
dmp.z = dmp.zd*dt + dmp.z;
dmp.y = dmp.yd*dt + dmp.y;
dmp.g = dmp.gd*dt + dmp.g;

y = dmp.y;
yd = dmp.yd;
ydd = dmp.ydd;

if dmp.A == 0
    max_y = max(max_y,dmp.y);
    min_y = min(min_y,dmp.y);
    if dmp.cs.x<0.0001
        dmp.A = max_y - min_y;
    end
end
end

%% DMP Batch Fit
function [dmp,Y,Yd,Ydd] = dmp_batch_fit(dmp,tau,traj,trajd,trajdd)
if ismatrix(traj) % 
    traj = traj';
    trajd=trajd';
    trajdd = trajdd';
else
traj = squeeze(traj);
trajd = squeeze(trajd);
trajdd = squeeze(trajdd);
end

if nargin<4
trajd = gradient(traj,dmp.cs.dt);
trajdd = gradient(trajd,dmp.cs.dt);
end

%start state & goal state
y0 = traj(1,:);
goal = traj(end,:);
g = goal;

%amplitude is the max(traj)-min(traj)
A = max(traj) - min(traj);

%hidden states
X = zeros(size(traj));
V = zeros(size(traj));
G = zeros(size(traj));
x = ones(1,dmp.n_dmps);
v = zeros(1,dmp.n_dmps);

for i = 1:length(traj)
    X(i,:) = x;
    V(i,:) = v;
    G(i,:) = g;
    vd = 0;
    xd = dmp.cs.ax.*(0-x)*tau;
    gd = (goal-g).*dmp.ag;
    
    x = xd*dmp.cs.dt + x;
    v = vd*dmp.cs.dt + v;
    g = gd*dmp.cs.dt + g;
end

%regression target
dmp.dG = goal - y0;
dmp.A = max(traj)-min(traj);
dmp.s = 1; %for fitting a new primitive, the scale is always one

amp = dmp.s;
Ft = (trajdd/tau^2 - dmp.ay.*(dmp.by.*(G-traj)-trajd/tau)) /amp;

%weights for each local model
PSI = exp(-0.5*((X*ones(dmp.n_dmps,length(dmp.c))-ones(length(traj),1)*dmp.c).^2)...
    .*(ones(length(traj),1)*dmp.h));
%compute regression
for a = 1:dmp.n_dmps
    dmp.sx2(a,:) = sum(((X(:,a).^2)*ones(1,length(dmp.c))).*PSI,1)';
    dmp.sxtd(a,:) = sum(((X(:,a).*Ft(:,a))*ones(1,length(dmp.c))).*PSI,1)';
    dmp.w(a,:) = (dmp.sxtd(a,:)./(dmp.sx2(a,:)+1e-10))';
end
%compute the prediction
F = sum((X*dmp.w).*PSI,2)./sum(PSI,2) * amp;

z = zeros(1,dmp.n_dmps);
zd = zeros(1,dmp.n_dmps);
y = y0;
Y = zeros(size(traj));
Yd = zeros(size(traj));
Ydd = zeros(size(traj));

for i = 1:length(traj)
    Ydd(i,:) = zd*tau;
    Yd(i,:) = z;
    Y(i,:) = y;
    zd = (dmp.ay.*(dmp.by.*(G(i,:)-y)-z)+F(i,:))*tau;
    yd = z;
    z = zd*dmp.cs.dt + z;
    y = yd*dmp.cs.dt + y;
end

end

%% Minimum Jerk Trajectory
function [x,xd,xdd] = min_jerk_step(x,xd,xdd,goal,tau,dt)
if tau<dt
	return
end

dist = goal - x;

a1   = 0;
a0   = xdd * tau^2;
v1   = 0;
v0   = xd * tau;

t1=dt;
t2=dt^2;
t3=dt^3;
t4=dt^4;
t5=dt^5;

c1 = (6.*dist + (a1 - a0)/2. - 3.*(v0 + v1))/tau^5;
c2 = (-15.*dist + (3.*a0 - 2.*a1)/2. + 8.*v0 + 7.*v1)/tau^4;
c3 = (10.*dist+ (a1 - 3.*a0)/2. - 6.*v0 - 4.*v1)/tau^3;
c4 = xdd/2.;
c5 = xd;
c6 = x;

x   = c1*t5 + c2*t4 + c3*t3 + c4*t2 + c5*t1 + c6;
xd  = 5.*c1*t4 + 4*c2*t3 + 3*c3*t2 + 2*c4*t1 + c5;
xdd = 20.*c1*t3 + 12.*c2*t2 + 6.*c3*t1 + 2.*c4;

end
