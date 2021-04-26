%% DMP_discrete: following the studywolf implementation
% Natasha Rouse (April 2021)
close all
clear all
clc
tic

%% Set parameters
dt = 0.01; %time increment
bfs = 10;
n_dmps = 2;
tau = 1; %time scale
% timesteps = 100;

%DESIRED TRAJECTORY
% desired_traj = [-5:0,0:-1:-5;-2*(-5:0),2*(0:-1:-5)]'; %zig zag
% desired_traj = [sin(-2*pi:0.1:2*pi); -2*pi:0.1:2*pi]'; %Upward sine wave
desired_traj = [sin(-2*pi:0.1:2*pi); 0:-0.1:-12.5]'; %Downward sine wave
figure(4)
line(desired_traj(:,1),desired_traj(:,2))

D = dmp_init(n_dmps,bfs,dt,desired_traj(1,:),desired_traj(end,:));

%CENTERS FOR BFS
% %From Schaal
% t = (0:1/(D.bfs-1):1)*0.5;
% D.c = exp(-D.cs.ax*t);
%From studywolf
c = linspace(0,D.cs.runtime,D.bfs);
for i = 1:length(c)
    D.c(i) = exp(-D.cs.ax*c(i));
end
% %From SMP paper
% D.c = logspace(log10(1), log10(.01), bfs);

%VARIANCE OF BFS (TRIAL AND ERROR)
%From studywolf
D.h = ones(1,D.bfs)*D.bfs^1.5 ./ D.c/D.cs.ax;
% %From Schaal
% D.h = (diff(D.c*0.5)).^2;
% D.h = 1./[D.h,D.h(end)];
% %From SMP paper
% D.h = linspace(log10(.3),log10(.002), bfs);

%Check Offset
if abs(D.y0-D.goal)<1e-4
    D.goal = D.goal + 1e-4;
end

[D,y_des] = imitate_path(D,desired_traj);
[D,y,yd,ydd] = dmp_rollout(D,tau);


figure(5)
 clf(5)
 figure(5)
plot(y(:,1),y(:,2))

toc
    
%% DMP initialization
function dmp = dmp_init(n_dmps,bfs,dt,y0,goal,w,ay,by)

%state variables
dmp.z       = 0;
dmp.y       = 0; 
dmp.x       = 0;
dmp.v       = 0;
dmp.zd      = 0;
dmp.yd      = 0; 
dmp.xd      = 0;
dmp.vd      = 0;
dmp.ydd     = 0;

dmp.n_dmps = n_dmps;
dmp.bfs = bfs;
dmp.dt = dt;
dmp.y0 = ones(1,n_dmps).*y0;
dmp.goal = ones(1,n_dmps).*goal;
if nargin<6
    dmp.w = zeros(n_dmps,bfs);
else
    dmp.w = w;
end
if nargin<7
    dmp.ay = ones(1,n_dmps)*25;
else
    dmp.ay = ay;
end
if nargin<8
    dmp.by = dmp.ay / 4;
else
    dmp.by = by;
end

%Canonical state variables
dmp.cs.ax = 1;
dmp.cs.runtime = 1;
dmp.cs.dt = dt;
dmp.cs.timesteps = dmp.cs.runtime / dmp.cs.dt;
dmp.cs.x = 1;

dmp = dmp_reset(dmp);
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

%% DMP reset
function dmp = dmp_reset(dmp)
dmp.y = dmp.y0;
dmp.yd = zeros(1,dmp.n_dmps);
dmp.ydd = zeros(1,dmp.n_dmps);
dmp.cs.x = 1;
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

