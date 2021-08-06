% Pendulum trajectory defined by DMP
% Natasha Rouse (August 2021)
% Trying this all again for variable l.

%% GOALS
% 1. Input a lemniscate (infinty sign) as the desired trajectory
% 2. Use DMP process to get initial guess for correct weights (consider an
% open loop system)
% 3. Run through Kati's dynamics and observe pendulum behaviour with Nate's
% plot

% After that...
% 4. Use a cost function to evaluate the dynamics performance ONE TIME.


%% SET-UP
close all
clear all
%Initialize quick-change variables
n_dmps = 2; %# of dmps (& degrees of freedom)
bfs = 10; %# of basis functions / kernels
dt = 0.01; %timestep size
runtime = 10;


%Trajectory (parametric equations for a lemniscate) -- these are the
%desired movements of the sampler
L = 1;
index = linspace(0,2*pi,100);
width = L;
height = .05*L;
y = width/2*sin(index)+.75*width/2;
z = height/2*sin(index*2)-L+height/2;
%--------------------at fixed length
theta = acos(-z/L);
yc = y-L*sin(theta);

traj = [yc; theta]';

figure(1)
hold on
plot(y,z)


figure(3)
hold on
plot(y,z)


%% Initialize DMP
% dmp_init(dmps,bfs,dt,start,goal,...)
D = dmp_init(n_dmps,bfs,dt,runtime,traj,[y(end) z(end)]);

%% Run DMP
[Y,y,cost,dynamics] = learn_dmp(D,0,traj);

%Use backwards kinematics to find the required cart movement and
%lengthening/theta

%For a variable length pendulum





%% DMP initialization
function dmp = dmp_init(n_dmps,bfs,dt,runtime,y,goal,w,ay,by)

dmp.n_dmps = n_dmps;
dmp.bfs = bfs;
dmp.y0 = ones(1,n_dmps).*y(1,:);

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

if nargin<7
    dmp.w = zeros(n_dmps,bfs);
else
    dmp.w = w;
end
if nargin<8
%     dmp.ay = ones(1,n_dmps)*25;
    dmp.ay = 25;
else
    dmp.ay = ay;
end
if nargin<9
    dmp.by = dmp.ay / 4;
else
    dmp.by = by;
end
dmp.ag = dmp.ay/2;

dmp.A = max(y)-min(y); %original amplitude (max(y)-min(y)) when the primitive was fit
dmp.s = 1; %scale factor for the nonlinear function
dmp.dG = dmp.G-dmp.y0; %original goal amplitude (G-y0) when the primitve was fit
dmp.sx2 = zeros(1,bfs);
dmp.sxtd = zeros(1,bfs);
dmp.lambda = 1;

%Canonical state variables
dmp.ax = dmp.ay/3;
dmp.runtime = runtime;
dmp.dt = dt;
dmp.timesteps = dmp.runtime / dmp.dt;
dmp.tau = .5;
% dmp.x = 1;

t = (0:1/(bfs-1):1)*0.5;
dmp.c = exp(-dmp.ax*t)';
% dmp.c = logspace(log10(1), log10(.01), bfs)';
dmp.h = (diff(dmp.c)*0.55).^2;
dmp.h = 1./[dmp.h;dmp.h(end)];
% dmp.h = logspace(log10(.3),log10(.002), bfs)';

dmp = dmp_reset(dmp);
end


%% Minimum Jerk Trajectory
function [x,xd,xdd] = min_jerk_step(x,xd,xdd,goal,tau,dt)
% function [x,xd,xdd] = min_jerk_step(x,xd,xdd,goal,tau, dt) computes
% the update of x,xd,xdd for the next time step dt given that we are
% currently at x,xd,xdd, and that we have tau until we want to reach
% the goal

if tau<dt
    'tau is less than dt'
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

%% DMP Batch Fit
function [dmp,Y,Yd,Ydd] = dmp_batch_fit(dmp,traj,trajd,trajdd)

if nargin<4
trajd = gradient(traj,dmp.dt);
trajdd = gradient(trajd,dmp.dt);
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
    xd = dmp.ax.*(0-x)*dmp.tau;
    if any(goal-g)
        'goal different from g'
    end
    gd = (goal-g).*dmp.ag;
    if any(gd)
        'gd is nonzero'
    end
    x = xd*dmp.dt + x;
    v = vd*dmp.dt + v;
    g = gd*dmp.dt + g;
end

%regression target
dmp.dG = goal - y0;
dmp.A = max(traj)-min(traj);
dmp.s = 1; %for fitting a new primitive, the scale is always one

amp = dmp.s;
Ft = (trajdd/dmp.tau^2 - dmp.ay.*(dmp.by.*(G-traj)-trajd/dmp.tau)) /amp;

%weights for each local model
PSI = exp(-0.5*((X*ones(dmp.n_dmps,length(dmp.c))-ones(length(traj),1)*dmp.c').^2)...
    .*(ones(length(traj),1)*dmp.h'));
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
    Ydd(i,:) = zd*dmp.tau;
    Yd(i,:) = z;
    Y(i,:) = y;
    zd = (dmp.ay.*(dmp.by.*(G(i,:)-y)-z)+F(i,:))*dmp.tau;
    yd = z;
    z = zd*dmp.dt + z;
    y = yd*dmp.dt + y;
end

end

%% DMP set goal
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
        'scaling is set according to amplitude'
        dmp.s = (dmp.G-dmp.y0)/dmp.dG;
    end
end
end

%% DMP reset
function dmp = dmp_reset(dmp,y)
if nargin < 2
    y = zeros(1,dmp.n_dmps);
end
% initialize the state variables
dmp.z       = zeros(1,dmp.n_dmps);
dmp.y       = y; 
dmp.x       = zeros(1,dmp.n_dmps);
dmp.v       = zeros(1,dmp.n_dmps);
dmp.zd      = zeros(1,dmp.n_dmps);
dmp.yd      = zeros(1,dmp.n_dmps);
dmp.xd      = zeros(1,dmp.n_dmps);
dmp.vd      = zeros(1,dmp.n_dmps);
dmp.ydd     = zeros(1,dmp.n_dmps);
% the goal state
dmp.G       = y;
dmp.g       = y;
dmp.gd      = zeros(1,dmp.n_dmps);
dmp.y0      = y;
dmp.s       = ones(1,dmp.n_dmps);

end


%% DMP run
function [dmp,y,yd,ydd,weights] = dmp_run(dmp)
%weighted sum of locally weighted regression models
dmp.psi = exp(-0.5*((dmp.x-dmp.c).^2).*dmp.h);
amp = dmp.s;
in = dmp.x;

f = sum(in*dmp.w.*dmp.psi',2)'/sum(dmp.psi+1e-10)*amp;

dmp.vd = zeros(1,dmp.n_dmps);
dmp.xd = dmp.ax*(0-dmp.x)*dmp.tau;
dmp.zd = (dmp.ay.*(dmp.by*(dmp.g-dmp.y)-dmp.z)+f)*dmp.tau;
dmp.yd = dmp.z*dmp.tau;
dmp.ydd = dmp.zd*dmp.tau;

dmp.gd = dmp.ag.*(dmp.G-dmp.g);

dmp.x = dmp.xd*dmp.dt + dmp.x;
dmp.v = dmp.vd*dmp.dt + dmp.v;
dmp.z = dmp.zd*dmp.dt + dmp.z;
dmp.y = dmp.yd*dmp.dt + dmp.y;

dmp.g = dmp.gd*dmp.dt + dmp.g;

y = dmp.y;
yd = dmp.yd;
ydd = dmp.ydd;
weights = dmp.psi.*in./sum(dmp.psi+1e-10).*amp;
end


%% Schaal "learn_dcp_batch" process
function [Y,y,cost,dynamics] = learn_dmp(dmp,flag,desired_traj,w)
if nargin==4
    dmp.w = w;
end

%Variables for plotting
Z = zeros(2*dmp.tau/dmp.dt+1,dmp.n_dmps);
Zd=Z;
X=Z;
Xd=X;
V=Z;
Vd=V;
T = zeros(2*dmp.tau/dmp.dt+1,dmp.n_dmps);
Td=T;
Tdd=T;
Y=T;
Yd=Y;
Ydd=Y;
PSI = zeros(dmp.bfs,dmp.n_dmps,(2*dmp.tau/dmp.dt+1));
W = zeros(dmp.bfs,dmp.n_dmps,(2*dmp.tau/dmp.dt+1));

if true(flag)
    %generate the minimum jerk trajectory
    t = desired_traj(1,:);
    td = 0;
    tdd = 0;
    for i = 0:2*dmp.tau/dmp.dt
        if i < length(desired_traj)
            [t,td,tdd] = min_jerk_step(t,td,tdd,desired_traj(i+1,:),dmp.tau-i*dmp.dt,dmp.dt);
        else
            [t,td,tdd] = min_jerk_step(t,td,tdd,desired_traj(end,:),dmp.tau-i*dmp.dt,dmp.dt);
        end    
        T(i+1,:) = t;
        Td(i+1,:) = td;
        Tdd(i+1,:) = tdd;
    end

    figure(1)
     clf(1)
     figure(1)
    hold on
    plot(T(:,1),T(:,2))
    plot(Td(:,1),Td(:,2),'.-')
    plot(Tdd(:,1),Tdd(:,2))
    axis equal
    title('Min Jerk Trajectory')

end


T = desired_traj;
[~,Td] = gradient(desired_traj,dmp.dt);
[~,Tdd] = gradient(Td,dmp.dt);

figure(1)
%  clf(1)
%  figure(1)
hold on
plot(desired_traj(:,1),desired_traj(:,2))
% plot(Td(:,1),Td(:,2))
% plot(Tdd(:,1),Tdd(:,2))
title('Gradient Differentiation Trajectory')
legend('desired trajectory (yc vs theta)')


%use batch_fit to initialize with minjerk
[dmp,~,~,~] = dmp_batch_fit(dmp,T,Td,Tdd);

%test the fit
dmp = dmp_reset(dmp);
dmp = dmp_set_goal(dmp,desired_traj(end,:),1);

%cost variables
% trajpointerror = inf*ones(size(desired_traj,1),1);
qdotdotsum = 0;
for i = 0:2*dmp.tau/dmp.dt
    [dmp,y,yd,ydd,weight] = dmp_run(dmp); % y is the actual trajectory of the end effector

    Z(i+1,:) = dmp.z; 
    Zd(i+1,:) = dmp.zd;
    Y(i+1,:) = y;
    Yd(i+1,:) = yd;
    Ydd(i+1,:) = ydd;
    X(i+1,:) = dmp.x;
    Xd(i+1,:) = dmp.xd;
    V(i+1,:) = dmp.v;
    Vd(i+1,:) = dmp.vd;
    PSI(:,:,i+1) = dmp.psi;
%     W(:,i+1,:) = D.w';
    W(:,:,i+1) = weight;
    
    l=1;
ys = y(1) + l*sin(y(2));
zs = -l*cos(y(2));


if i < length(desired_traj)
    trajdistancessq = ((ys-desired_traj(i+1,1)).^2 + ...
    (zs-desired_traj(i+1,2)).^2);
else
    trajdistancessq = ((ys-desired_traj(end,1)).^2 + ...
    (zs-desired_traj(end,2)).^2);
end
trajpointerror(i+1) = trajdistancessq;
end
    
%INSERT DYNAMICS/KINEMATICS HERE
l=1;
ys = Y(:,1) + l*sin(Y(:,2));
zs = -l*cos(Y(:,2));
ycart = Y(:,1);
end_effector = [ys zs];

figure(3)
%  clf(3)
%  figure(3)
hold on 
plot (ys,zs, '.b-') %produced trajectory
plot(desired_traj(:,1),desired_traj(:,2),'m')
% plot(Y(:,1),Y(:,2))
plot (ycart, ycart*0, 'k','Linewidth',1.5) %cart movement (horizontal)
plot([ycart(1),ys(1)], [0,zs(1)], '-k')%initial
[~,i] = max(ys);
plot([ycart(i),ys(i)], [0,zs(i)], '-b')%max 
axis equal
title('sampler trace')
% legend ('produced traj', 'desired traj')

(end_effector(end)-dmp.g).^2
100*sum(sqrt((end_effector-dmp.g).^2))
cost = 1000*sum(sqrt(trajpointerror)) + 100*sum(sqrt((end_effector-dmp.g).^2)) + ...
    sum(rssq(Yd,2))*dmp.dt + qdotdotsum*dmp.dt




%     u_hist(:,i+1) = u;
    q_hist = q;
%     q_dot_hist(:,i+1) = q_dot;
%     q_ddot_hist(:,i+1) = q_ddot;
    qdotdotsum = rssq(q_ddot) + qdotdotsum;
%     end_effector_hist(:,i+1) = end_effector;
%     EE_accel_hist(:,i+1) = end_effector_accel; 


%    if mod(i,10)==0
       figure(2);
       plot([ycart(1),ys(1)], [0,zs(2)]);
       axis([-6 6 -6 6]);

%     frame = getframe(gcf);

       drawnow
%     writeVideo(v, frame);
%    end
    pause(0.1)
% 
%     
% % cost_components = [1000*sum(sqrt(trajpointerror))
% %     100*sum(sqrt((end_effector-D.g).^2))
% %     sum(rssq(q_dot_hist,1))*D.dt
% %     qdotdotsum*D.dt];

% PLOTTING
time = (0:dmp.dt:2*dmp.tau)';

figure(4)
 clf(4)
 figure(4)
hold on
plot(time,Y(:,1))
plot(time,Y(:,2))

figure(5)
 clf(5)
 figure(5)
for i = 1:dmp.n_dmps
    
    % position, velocity, acceleration vs. target
    subplot(4,3,1)
    hold on
    plot(time,Y(:,i)); % DMP trained position
    title('y')

    subplot(4,3,2)
    hold on
    plot(time,Yd(:,i)); % DMP trained velocity
    title('yd')

    subplot(4,3,3)
    hold on
    plot(time,Ydd(:,i)); % DMP trained acceleration
    title('ydd')

%     internal states
    subplot(4,3,4)
    hold on
    plot(time,Z(:,i))
    title('z')

    subplot(4,3,5)
    hold on
    plot(time,Zd(:,i))
    title('zd')

    %Updated plot for this slot
    subplot(4,3,6)
    hold on
    for j = 1:dmp.bfs
        plot(time,squeeze(PSI(j,:,i))*squeeze(W(j,:,i))')
        title('{Psi activations/},{Basis functions}')
    end

    %The original plot for this slot
    subplot(4,3,6)
    hold on
    plot(time,squeeze(PSI(:,i,:)))
    title('{Psi activations/},{Basis functions}')

    subplot(4,3,7)
    hold on
    plot(time,V(:,i))
    title('v')

    subplot(4,3,8)
    hold on
    plot(time,Vd(:,i))
    title('vd')

    subplot(4,3,9)
    hold on
    for j = 1:dmp.bfs
        plot(time,squeeze(W(j,i,:)))
        title('Linear Model Weights over Time')
    end

    subplot(4,3,10)
    hold on
    plot(time,X(:,i))
    title('x')

    subplot(4,3,11)
    hold on
    plot(time,X(:,i))
    title('xd')

    subplot(4,3,12)
    hold on
    plot(squeeze(W(end,i,:)))
    title('Weights')
    xlabel(sprintf('tau=%f',dmp.tau))

    drawnow
end


end


