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
sigma = linspace(log10(.3),log10(.002), bfs); %variance of basis functions
c = logspace(log10(1), log10(.01), bfs); %spacing of basis functions in time
tau = 0.5; %time scale
timesteps = 100;
%desired trajectory for end-effector: "spiral" downward starting from top
%left
desired_traj = [-5:0,0:-1:-5;-2*(-5:0),2*(0:-1:-5)]';
% line(desired_traj(:,1),desired_traj(:,2))

D = dmp_init(n_dmps,bfs,dt,desired_traj(1,:),desired_traj(end,:));
%generate centres of bfs
des_c = linspace(0,D.cs.runtime,D.bfs);
for i = 1:length(des_c)
    D.c(i) = exp(-D.cs.ax*des_c(i));
end

%set variance of bfs (trial and error to find this)
D.h = ones(1,D.bfs)*D.bfs^1.5 ./ D.c/D.cs.ax;
%check offset
if abs(D.y0-D.goal)<1e-4
    D.goal = D.goal + 1e-4;
end

y_des = imitate_path(D,desired_traj);
[y,yd,ydd] = dmp_rollout(D,timesteps,tau);


toc
    
%% DMP initialization
function dmp = dmp_init(n_dmps,bfs,dt,y0,goal,w,ay,by)

% dmp{50} = struct();
% %time constants for critical damping
% dmp.alpha_z = 25;
% dmp.beta_z  = dmp.alpha_z/4;
% dmp.alpha_g = dmp.alpha_z/2;
% dmp.alpha_x = dmp.alpha_z/3;
% dmp.alpha_v = dmp.alpha_z;
% dmp.beta_v  = dmp.beta_z;
% 
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
% 
% %set goal state
% dmp.g       = 0;
% dmp.gd      = 0;
% dmp.G       = 0;
% 
% %current start state
% dmp.y0 = 0;
% % the orginal amplitude (max(y)-min(y)) when the primitive was fit
% dmp.A       = 0;
% % the original goal amplitude (G-y0) when the primitive was fit
% dmp.dG      = 0;
% % the scale factor for the nonlinear function
% dmp.s       = 1;
% 
% t = (0:1/(bfs-1):1)'*0.5;
% % the local models, spaced on a grid in time by applying the
% % anaytical solutions x(t) = exp(-alpha*t)
% dmp.c       = exp(-dmp.alpha_x*t);
% % we also store the phase velocity at the centers which is used by some
% % applications: xd(t) = x(t)*(-dmp.alpha_x);
% dmp.cd      = dmp.c*(-dmp.alpha_x);
% 
dmp.n_dmps = n_dmps;
dmp.bfs = bfs;
% dmp.psi     = zeros(dmp.bfs,1);
% dmp.w       = zeros(dmp.bfs,1);
% dmp.sx2     = zeros(dmp.bfs,1);
% dmp.sxtd    = zeros(dmp.bfs,1);
% dmp.D       = (diff(dmp.c)*0.55).^2;
% dmp.D       = 1./[dmp.D;dmp.D(end)];
% dmp.lambda  = 1;
dmp.dt = dt;
dmp.y0 = ones(n_dmps,1)*y0;
dmp.goal = ones(n_dmps,1)*goal;
if nargin<6
    dmp.w = zeros(n_dmps,bfs);
else
    dmp.w = w;
end
if nargin<7
    dmp.ay = ones(n_dmps,1)*25;
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
end

%% DMP rollout
function [y_track,dy_track,ddy_track] = dmp_rollout(dmp,timesteps,tau)
dmp = dmp_reset(dmp);

if (tau)
    timesteps = timesteps/tau;
end

%set up tracking vectors
y_track = zeros(timesteps,dmp.bfs);
dy_track = zeros(timesteps,dmp.bfs);
ddy_track = zeros(timesteps,dmp.bfs);

for t = 1:timesteps
    %run and record timesteps
    [y_track(t),dy_track(t),ddy_track(t)] = dmp_step(dmp,0,0);
end
end

%% DMP reset
function dmp = dmp_reset(dmp)
dmp.y = dmp.y0;
dmp.yd = zeros(dmp.n_dmps,1);
dmp.ydd = zeros(dmp.n_dmps,1);
dmp.cs.x = 1;
end

%% DMP step
function [y,yd,ydd] = dmp_step(dmp,error,external_force)
%run the DMP system for a single timestep
tau = 1; %scales the timestep. Increase tau to make the system execute faster
error_coupling = 1/(1+error); %error is optional system feedback

%run canonical system
% x = dmp.cs.x + (-dmp.cs.ax*dmp.cs.x*error_coupling)*tau*dmp.cs.dt;
% dmp.cs.x = x;
x = cs_step(dmp,tau,error_coupling);

%generate basis function activation
psi = exp(-dmp.h .* (x-dmp.c).^2);
for i = 1:dmp.n_dmps
    %generate diminishing front term on the forcing term
    f = (x * (dmp.goal(i) - dmp.y0(i)))*dot(psi,dmp.w(i,:));
    if abs(sum(psi)) > 1e-6
        f = f/sum(psi);
    end
    %dmp acceleration
    dmp.ydd = dmp.ay*(dmp.by*dmp.goal - dmp.y - dmp.yd) + f;
    if (external_force)
        dmp.ydd(i) = dmp.ydd(i) + external_force(i);
    end
    dmp.yd(i) = dmp.yd(i) + dmp.ydd(i)*tau*dmp.dt*error_coupling;
    dmp.y(i) = dmp.y(i) + dmp.yd(i)*tau*dmp.dt*error_coupling;
end
y = dmp.y(i);
yd = dmp.yd(i);
ydd = dmp.ydd(i);
end

%% Canonical System step
function [x_state] = cs_step(dmp,tau,error_coupling)
% generate a single step of x for discrete (potentially closed) loop
% movements. Decaying from 1 to 0 according to dx = -ax*x.

%tau: gain on execution time. increase to make system execute faster.
%error_coupling: slow down if the error is > 1
% (-dmp.cs.ax*dmp.cs.x*error_coupling)*tau*dmp.cs.dt
x_state = dmp.cs.x + (-dmp.cs.ax*dmp.cs.x*error_coupling)*tau*dmp.cs.dt;
dmp.cs.x = x_state;
end

%% CS rollout
function [x_track] = cs_rollout(dmp,tau,error_coupling)
if nargin == 1
    tau = 1;
    error_coupling = 1;
end
x_track = zeros(dmp.cs.timesteps,1);
%reset CS
dmp.cs.x = 1;
for t = 1:dmp.cs.timesteps
    x_track(t) = dmp_step(dmp,tau,error_coupling);
end
end

%% Imitate path
function y_des = imitate_path(dmp,y_des)
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
    for t = 1:length(dmp.cs.timesteps)
        path(i,t) = interp1(x,y_des(:,i),t*dmp.cs.dt);
    end
end
y_des = path;

%calculate velocity of y_des with central difference
dy_des = (gradient(y_des))/dmp.cs.dt;
%calculate acceleration of y_des
ddy_des = gradient(dy_des)/dmp.cs.dt;

f_target = zeros(length(y_des),dmp.n_dmps);
%find the force required to move along this trajectory
for k = 1:dmp.n_dmps
    f_target(:,k) = ddy_des(k) - dmp.ay(k)*(dmp.by(k)*(dmp.goal(k)-y_des(k))...
        -dy_des(k));
end

%efficiently generate weights to realize f_target
dmp = generate_weights(dmp,f_target);

x_track = cs_rollout(dmp);
psi_track = exp(-dmp.h .* (x_track-dmp.c).^2);

%Show Plot
figure(1)
 clf(1)
 figure(1)
%Basis function activations
subplot(2,1,1)
subtitle('Basis functions')
plot(psi_track)
% hold on
% for num = 1:size(psi_track,2)
%     plot(psi_track(:,num))
% end

%Plot desired forcing function vs approx
subplot(2,1,2)
subtitle('DMP forcing function')
hold on
for n = 1:dmp.n_dmps
    plot(f_target(:,n),'--')
    plot(sum(psi_track*dmp.w(n))*dmp.cs.dt)
end
legend('f_{target}`','w*psi')
end

%% Generate Weights
function dmp = generate_weights(dmp,f_target)
%calculate x & psi
x_track = cs_rollout(dmp);
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

