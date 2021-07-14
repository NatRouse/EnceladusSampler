%% DMP_discrete: following the studywolf implementation
% Natasha Rouse (Edited June 2021)
% Including more functions from original Schaal code
% Init & Batch_Fit based on Schaal "learn_dcp_batch"
% PLEASE REPLACE DYNAMICS AS YOU SEE FIT!!!! ALL OTHER DMP CODE SHOULD BE
% FUNCTIONAL!!!
close all
clear all
clc

tic

%% Set parameters
dt = 0.01; %time increment
bfs = 10;
n_dmps = 2;
tau = 1;
traj_length = 2*tau/dt+1;
%DESIRED TRAJECTORY
% desired_traj = [-5:0,0:-1:-5;-2*(-5:0),2*(0:-1:-5)]'; %zig zag
% desired_traj = [sin(-2*pi:0.1:2*pi); -2*pi:0.1:2*pi]'; %Upward sine wave
% n_traj = 2*tau/dt+2;
% desired_traj = [0:dt:dt*(traj_length-1)]'; %zeros(1, 126); zeros(1, 126);zeros(1, 126)]'; %Downward sine wave
desired_traj = [0:0.0001:0.0001*(traj_length-1);zeros(1,traj_length)]'; %
% figure(1)
% line(desired_traj(:,1),desired_traj(:,2))

D = dmp_init(n_dmps,bfs,dt,desired_traj(1,:),desired_traj(end,:));

%CENTERS FOR BFS
t = (0:1/(D.bfs-1):1)*0.5;
D.c = exp(-D.cs.ax(1)*t);
D.cd = (-D.cs.ax(1))*D.c;

%VARIANCE OF BFS (TRIAL AND ERROR)
% D.h = (diff(D.c*0.5)).^2;
% D.h = 1./[D.h,D.h(end)];
D.h = exp(-D.ay/3.*t/tau);

[Y,y,cost] = learn_dmp(D,tau,desired_traj);

% dyn.u = u_hist;
% dyn.q = q_hist;
% dyn.qd = q_dot_hist;
% dyn.qdd = q_ddot_hist;
% dyn.ee = end_effector_hist;
% dyn.ee_accel = EE_accel_hist;
% close(v)

%% Cost function evaluation
% cost = Eval(D,dyn,desired_traj);
d = 1;
allcosts = [cost];
for passes = 1:50
%     'enter cost evaluation'
    passes
    bestw = D.w;
    bestcost = cost;
    counter1 = inf;
    counter2 = inf;
    while counter1+counter2 >=2
        counter1 = 0;
        counter2 = 0;
        for i = 1:size(D.w,1)
            for j = size(D.w,2)
                wtrial = bestw;
                wnew = bestw;
                if rem(i,2)==1 || rem(j,2)==1
                    delta = d*(1-2*(rand(1)));
                else
                    delta = d*(1+2*rand(1));
                end
                wtrial(i,j) = bestw(i,j) + delta;
                [~,~,costtrial] = learn_dmp(D,tau,desired_traj,wtrial);
                wnew(i,j) = bestw(i,j) + max(min(bestcost*delta/(bestcost-costtrial), d*5), -d*5);
                [~,~,costnew] = learn_dmp(D,tau,desired_traj,wnew);
                if costnew < bestcost
                    bestw = wnew;
                    bestcost = costnew;
                    counter1 = counter1 + 1;
                end
                if costtrial < bestcost
                    bestw = wtrial;
                    bestcost = costtrial;
                    counter2 = counter2 + 1;
                end
            end
        end
    end
    allcosts = [allcosts,bestcost];
    D.w = bestw;
    cost = bestcost;
end






%% plot
time = (0:D.dt:tau*2)';
figure(3);
hold on
% time = 0:dt:dt*(size(desired_traj, 1)-1);
subplot(2,1,1)
plot(time, squeeze(Y(1,:,1)), time, end_effector_hist(1,:)-end_effector_hist(1,1))
title('Desired vs  Y Trajectory')
legend('Desired Y', 'End Effector Y', 'Location', 'southeast')

subplot(2,1,2)
hold on
plot(time, squeeze(Y(1,:,2)), time, end_effector_hist(2,:)-end_effector_hist(2,1))
title('Desired vs Z Trajectory')
legend('Desired Z', 'End Effector Z')

% plot configuration space
hold off
figure(4);
plot(time, q_hist(1,:), time, q_hist(2,:), time, q_hist(3,:))%, time, q_hist(3,:))
legend('y_c', 'l', 'theta')
title('configuration space')

% plot actuator forces
figure(5);
hold on
subplot(3,1,1)
plot(time, u_hist(1,:))
title('y_c actuator')
subplot(3,1,2)
plot(time, u_hist(2,:))
title('l actuator')
subplot(3,1,3)
plot(time, u_hist(3,:))
title('theta actuator')

hold off

% plot desired and actual accelerations
figure(6); 
subplot(2,1,1)
plot(time, squeeze(Y(3,:,1)), time, EE_accel_hist(1,:))
title('desired vs actual acceleration y')
legend('desired', 'actual')
subplot(2,1,2)
plot(time, squeeze(Y(3,:,2)), time, EE_accel_hist(2,:))
title('desired vs actual acceleration z')
legend('desired', 'actual')

figure(7)
hold on
title('All Costs')
plot(allcosts,'-k.','Linewidth',1.5)
xlabel('Iterations')
ylabel('Cost')

toc

%%
function [u, q_ddot, q_dot, q, trajectory, end_effector_accel] = robot_dynamics(x_ddot, q, q_dot,m_s, m_c, dt)
          
    y_c = q(1);
    l = q(2);
    theta = q(3);
    
    y_c_dot = q_dot(1);
    l_dot = q_dot(2);
    theta_dot = q_dot(3);
    
    %Jacobian
    J = [1 sin(theta) l*cos(theta);
         0  -cos(theta) l*sin(theta)];
%         0  0 1
%         0  0 0
%         0  0 0];
    % time rate of change of Jacobian
    J_dot = [0 theta_dot*cos(theta) l_dot*cos(theta)-l*theta_dot*sin(theta); 
             0 theta_dot*sin(theta) l_dot*sin(theta)+l*theta_dot*cos(theta)];
    
%              0 0 0;
%              0 0 0;
%              0 0 0];
    % mass matrix from equations of motion
    g = 9.81; % 

    M = [m_s+m_c                m_s*sin(theta)   m_s*l*cos(theta);
         m_s*sin(theta)         m_s              0
         m_s*l*cos(theta)        0               m_s*l^2];
    % damping matrix from EOMs
    C = [0                      2*m_s*theta_dot*cos(theta)   -m_s*l*theta_dot*sin(theta); 
        0                             0                       -m_s*l*theta_dot
        0                       2*m_s*l*theta_dot            0];
    K = [0 0 0;
        0 0 0;
        0 m_s*g*sin(theta) 0];

    % stiffness matrix from EOMS
    Mx = inv(J*inv(M)*pinv(J))*(x_ddot-J_dot*q_dot+J*inv(M)*C*q_dot+J*inv(M)*K*q)*pinv(x_ddot);
%     Mx2 = J*M*pinv(J)-J*M*pinv(J)*J_dot*q_dot*pinv(x_ddot)+J*C*q_dot*pinv(x_ddot)+J*K*q*pinv(x_ddot); 

    % total force vector in workspace coordinates
    F_x = Mx*x_ddot;%+m_s*l*g*sin(theta);
    
    % Actuation forces (configuration space)
    F_q = pinv(J)*F_x; % actuation force required
%     disp(F_q(end))
    u = F_q - [0; -g*cos(theta); 0];
    % joint accelerations
        q_ddot = inv(M)*(F_q-C*q_dot-K*q);
    
   

    % update configuration 
    q_dot = q_ddot*dt + q_dot;
    q = q_dot*dt + q; %1/2*q_ddot*dt^2 + 
    
     end_effector_accel = [q_ddot(1)+q_ddot(2)*sin(q(3))+2*q_dot(2)*q_dot(3)*cos(q(3)) + q(2)*q_ddot(3)*cos(q(3))-q(2)*q_dot(3)^2*sin(q(3));
    -q_ddot(2)*cos(q(3)) + 2*q_dot(2)*q_dot(3)*sin(q(3)) + q(2)*q_ddot(3)*sin(q(3))+q(2)*q_dot(3)^2*cos(q(3))];
    
%     if any(abs(end_effector_accel-x_ddot)>=[.1;.1])
%         disp('mismatch')
%         
%     end
    
    y = q(1)+q(2)*sin(q(3));
    z = -q(2)*cos(q(3));
    trajectory = [y;z];
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
traj = squeeze(traj);
trajd = squeeze(trajd);
trajdd = squeeze(trajdd);

if size(traj,1) == 1
    traj = traj';
    trajd = trajd';
    trajdd = trajdd';
end
if nargin<4
trajd = gradient(traj,dmp.cs.dt);
trajdd = gradient(trajd,dmp.cs.dt);
end

%start state & goal state
y0 = traj(:,1);
goal = traj(:,end);
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
    X(:,i) = x;
    V(:,i) = v;
    G(:,i) = g;
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
PSI = exp(-0.5*((X'*ones(dmp.n_dmps,length(dmp.c))-ones(length(traj),1)*dmp.c).^2)...
    .*(ones(length(traj),1)*dmp.h));
%compute regression
for a = 1:dmp.n_dmps
    dmp.sx2(a,:) = sum(((X(a,:)'.^2)*ones(1,length(dmp.c))).*PSI,1)';
    dmp.sxtd(a,:) = sum(((X(a,:).*Ft(a,:))'*ones(1,length(dmp.c))).*PSI,1)';
    dmp.w(a,:) = (dmp.sxtd(a,:)./(dmp.sx2(a,:)+1e-10))';
end
%compute the prediction
F = sum((X'*dmp.w).*PSI,2)./sum(PSI,2) * amp;

z = zeros(dmp.n_dmps,1);
zd = zeros(dmp.n_dmps,1);
y = y0;
Y = zeros(size(traj));
Yd = zeros(size(traj));
Ydd = zeros(size(traj));

for i = 1:length(traj)
    Ydd(:,i) = zd*tau;
    Yd(:,i) = z;
    Y(:,i) = y;
    zd = (dmp.ay.*(dmp.by.*(G(:,i)-y)-z)+F(i,:))*tau;
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

% %% Cost Evaluation Function
% function cost = Eval(dmp,dynamics,desiredtraj,w)
% 'inside eval function'
% if nargin == 3
%     w = dmp.w;
% end
% 
% goal_tolerance = 0.01;
% desiredtraj = desiredtraj(1:end-1,:)';
% 
% size(dynamics.ee)
% size(desiredtraj)
% % cart_position 
% distance = 0;
% timetogoal = nan;
% 
% % theta_deviation = ;
% % endpoint_error = ;
% % velocity = ;
% % acceleration = ;
% end


%% Schaal "learn_dcp_batch" process
function [Y,y,cost] = learn_dmp(D,tau,desired_traj,w)
if nargin==4
    D.w = w;
end
%Variables for plotting
Z = zeros(2,D.n_dmps,floor(2*tau/D.dt+1));
X = Z;
V=Z;
T = zeros(3,D.n_dmps,floor(2*tau/D.dt+1));
Y=T;
PSI = zeros(D.bfs,D.n_dmps,(2*tau/D.dt+1));
W = zeros(D.bfs,D.n_dmps,(2*tau/D.dt+1));
    
% Dynamics
% initial configuration
q = [0,1,0]'; % make sure the desired trajectory and initial configuration are compatible
q_dot = [0, 0, 0]'; % make sure this is is kinematically feasible given the desired trajectory

u_hist = nan(length(q), size(Y,2));
q_ddot_hist = nan(length(q), size(Y,2));
q_ddot_hist(:,1) = zeros(size(q));
q_hist = nan(length(q), size(Y,2));
q_hist(:,1) = q;

q_dot_hist = nan(length(q), size(Y,2));
q_dot_hist(:,1) = q_dot;

% figure(2)
%     v = VideoWriter('DMP Variable Length Fully Actuated');
%     v.FrameRate = 10;
%     open(v)
end_effector_hist = nan(2, size(Y,2));
end_effector_hist(:,1) = [q(1)+q(2)*sin(q(3)); -q(2)*cos(q(3))];
m_s = 10;
m_c = 10;

%generate the minimum jerk trajectory
t = 0; % 
td = 0; % 
tdd = 0; % 
for i = 0:2*tau/D.dt
    [t,td,tdd] = min_jerk_step(t,td,tdd,desired_traj(i+1,:),tau-i*D.dt,D.dt);
    T(:,:,i+1) = [t; td; tdd]; 
end

% Constant velocity T
% desiredt = linspace(0,10,length(desired_traj))';% assumes constant velocity
% desiredydot = [0,0;diff(desired_traj)./(diff(desiredt)*[1,1])]; % desired speed
% desiredydotdot = [0,0;diff(desiredydot)./(diff(desiredt)*[1,1])];
% T = cat(3,desired_traj,desiredydot,desiredydotdot);
% T = permute(T,[3 1 2]);
% T = T(:,1:end-1,:);

%use batch_fit to initialize with minjerk
[D,~,~,~] = dmp_batch_fit(D,tau,T(1,:,:),T(2,:,:),T(3,:,:));

%test the fit
% D = dmp_reset(D);
D = dmp_set_goal(D,desired_traj(end,:),1);

%cost variables
trajpointerror = inf*ones(size(desired_traj,1),1);
qdotdotsum = 0;
for i = 0:2*tau/D.dt
    [D,y,yd,ydd,weight] = dmp_run(D,tau); % y is the actual trajectory of the end effector

    Z(:,:,i+1) = [D.z; D.zd]; 
    Y(:,:,i+1) = [y; yd; ydd];
    X(:,:,i+1) = [D.x; D.xd];
    V(:,:,i+1) = [D.v; D.vd];
    PSI(:,:,i+1) = D.psi';
%     W(:,i+1,:) = D.w';
    W(:,i+1,:) = weight;
    
    ydd_i = squeeze(Y(3,:,i+1)');
    [u,q_ddot, q_dot, q, end_effector, end_effector_accel] = robot_dynamics(ydd_i, q, q_dot, m_s, m_c, D.dt);
    u_hist(:,i+1) = u;
    q_hist(:,i+2) = q;
    q_dot_hist(:,i+2) = q_dot;
    q_ddot_hist(:,i+2) = q_ddot;
    qdotdotsum = rssq(q_ddot) + qdotdotsum;
    end_effector_hist(:,i+2) = end_effector;
    EE_accel_hist(:,i+2) = end_effector_accel; 


%    if mod(i,10)==0
%        figure(2);
%        plot([q(1), end_effector(1)], [0,end_effector(2)]);
%        axis([-6 6 -6 6]);
% 
% %     frame = getframe(gcf);
% 
%        drawnow
% %     writeVideo(v, frame);
%    end
%     pause(0.1)
    trajdistancessq = min((end_effector(1)-desired_traj(i+1,1)).^2 + ...
        (end_effector(2)-desired_traj(i+1,2)).^2, trajpointerror(i+1));
    trajpointerror(i+1) = trajdistancessq;

end
% cost_components = [1000*sum(sqrt(trajpointerror))
%     100*sum(sqrt((end_effector-D.g).^2))
%     sum(rssq(q_dot_hist,1))*D.dt
%     qdotdotsum*D.dt];
cost = 1000*sum(sqrt(trajpointerror)) + 100*sum(sqrt((end_effector'-D.g).^2)) + ...
    sum(rssq(q_dot_hist,2))*D.dt + qdotdotsum*D.dt;

% PLOTTING
time = (0:D.dt:tau*2)';

figure(1)
clf(1)
figure(1)

for i = 1:D.n_dmps 
    hold on

    % position, velocity, acceleration vs. target
    subplot(4,3,1)
    plot(time,[squeeze(Y(1,i,:)) squeeze(T(1,i,:))]); % DMP trained position
    title('y')

    subplot(4,3,2)
    plot(time,[squeeze(Y(2,i,:)) squeeze(T(2,i,:))]); % DMP trained velocity
    title('yd')

    subplot(4,3,3)
    plot(time,[squeeze(Y(3,i,:)) squeeze(T(3,i,:))]); % DMP trained acceleration
    title('ydd')

%     internal states
    subplot(4,3,4)
    plot(time,squeeze(Z(1,i,:)))
    title('z')

    subplot(4,3,5)
    plot(time,squeeze(Z(2,i,:)))
    title('zd')

%     %Updated plot for this slot
%     subplot(4,3,6)
%     hold on
%     for j = 1:D.bfs
%         plot(time,squeeze(PSI(j,:,i))*squeeze(W(j,:,i))')
%         title('{Psi activations/},{Basis functions}')
%     end

    %The original plot for this slot
    subplot(4,3,6)
    plot(time,squeeze(PSI(:,i,:)))
    title('{Psi activations/},{Basis functions}')

    subplot(4,3,7)
    plot(time,squeeze(V(1,i,:)))
    title('v')

    subplot(4,3,8)
    plot(time,squeeze(V(2,i,:)))
    title('vd')

    subplot(4,3,9)
    hold on
    for j = 1:D.bfs
        plot(time,squeeze(W(j,i,:)))
        title('Linear Model Weights over Time')
    end

    subplot(4,3,10)
    plot(time,squeeze(X(1,i,:)))
    title('x')

    subplot(4,3,11)
    plot(time,squeeze(X(2,i,:)))
    title('xd')

    subplot(4,3,12)
    plot(squeeze(W(end,i,:)))
    title('Weights')
    xlabel(sprintf('tau=%f',tau))

    drawnow
end

end

