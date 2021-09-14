%% SET-UP
close all
clear all
rng(1)
%Initialize quick-change variables
n_dmps = 2; %# of dmps (& degrees of freedom)
bfs = 20; %# of basis functions / kernels
dt = 0.001; %timestep size
runtime = 10;
plotflag = 0;

%Trajectory (parametric equations for a lemniscate) -- these are the
%desired movements of the sampler
L = 1;
s2 = linspace(-pi/4,pi/4,200);
s3 = linspace(pi/4,3*pi/4,100);
s4 = linspace(3*pi/4,7*pi/4,200);
index = [s2,s3,s4];
% index = linspace(-pi/4,7*pi/4,500);

width = L;
height = .05*L;
y = width/2*sin(index)+.75*width/2;
z = height/2*sin(index*2)-L+height/2;
theta = acos(-z/L);
yc = y-L*sin(theta);

lem = [y; z]';
traj = [yc; theta]';
if traj(1,:)==traj(end,:)
    traj(end,:) = traj(end,:)+0.001;
end

c2 = ones(length(s2),1)*[1 0 1];
c3 = ones(length(s3),1)*[0 0 1];
c4 = ones(length(s4),1)*[0 1 1];
trajcolor = [c2;c3;c4];
% trajcolor = ones(length(index),1)*[0 0 1];

%Plot the desired traj (y-z & yc-theta)
figure(1)
 clf(1)
 figure(1)
 subplot(121)
hold on
for i = 1:length(y)
    plot(y(i),z(i),'.-','Color',trajcolor(i,:))
end
% axis equal
title('Trajectory as defined')

 subplot(122)
hold on
for i = 1:length(y)
    plot(yc(i),theta(i),'.-','Color',trajcolor(i,:))
end
title('Desired Traj')
xlabel('y_{cart}')
ylabel('theta')

 
%Initialize DMP
D = dmp_init(n_dmps,bfs,dt,runtime,traj);

%------------------------------------RUN DMP
Y=zeros(2*D.tau/dt+1,n_dmps);
Yd=Y;
Ydd=Y;

[~,trajd] = gradient(traj,dt);
trajd = [0,0;trajd(1:end-1,:)];
[~,trajdd] = gradient(trajd,dt);
trajdd = [0,0;trajdd(1:end-1,:)];

%Batch fit
[D,Y,Yd,Ydd] = dmp_batch_fit(D,traj,trajd,trajdd,plotflag);

%test the fit
D = dmp_reset(D);
D = dmp_set_goal(D,traj(end,:),1);

[D,Y,Yd,Ydd] = dmp_run_fit(D,plotflag);

initialcost = check_fit(D,lem,Y,Yd,Ydd);

costplot = [];
costplot = [costplot,initialcost];
cost = initialcost;

delta = 1;
%Let's learn!
for passes = 1:10
    bestw = D.w;
    bestcost = cost;
    %Trying random stuff first
    for iter = 1:50
        if rand(1)>0.5
            w = 1-(0.5*rand(2,bfs));
            [D,y,yd,ydd] = dmp_run_fit(D,plotflag,w);
            cost = check_fit(D,lem,y,yd,ydd);
            if cost < bestcost
                bestw = D.w;
                bestcost = cost;
                delta = 1;
            end
        else
            w = bestw + 3*delta*randn(2,bfs);
            [D,y,yd,ydd] = dmp_run_fit(D,plotflag,w);
            cost = check_fit(D,lem,y,yd,ydd);
            if cost < bestcost
                bestw = D.w;
                bestcost = cost;
            end
        end
    end
    
    D.w = bestw;
    cost = bestcost;
    delta = delta*0.99;
    [D,y,yd,ydd] = dmp_run_fit(D,plotflag,w);
    cost = check_fit(D,lem,y,yd,ydd);
    costplot = [costplot,bestcost];
    
    y
    plot_pendulum(y)
    pause
        
    %gradient descent
    choosenew = inf;
    choosedel = inf;
    while choosenew + choosedel >= 2
        choosenew = 0;
        choosedel = 0;
        for i = 1:size(D.w,1)
            for j = 1:size(D.w,2)
                %Vary weights (plus/minus some delta)
                wtrial = bestw;
                wnew = bestw;
                d = delta*(1-2*rand(1));
                wtrial(i,j) = bestw(i,j) + d;

                %Check cost - higher or lower?
                [D,y,yd,ydd] = dmp_run_fit(D,plotflag,wtrial);
                trialcost = check_fit(D,lem,y,yd,ydd);

                wnew(i,j) = bestw(i,j) + max(min(bestcost*d/(bestcost-trialcost),...
                    delta*5), -delta*5);
                [D,y,yd,ydd] = dmp_run_fit(D,plotflag,wnew);
                newcost = check_fit(D,lem,y,yd,ydd);

                if newcost < bestcost
                    bestw = wnew;
                    bestcost = newcost;
                    choosenew = choosenew + 1;
                end
                if trialcost < bestcost
                    bestw = wtrial;
                    bestcost = trialcost;
                    choosedel = choosedel + 1;
                end
            end
        end
    end
    costplot = [costplot,bestcost];
    plot_pendulum(y)
    pause(1)
end

figure(5)
 clf(5)
 figure(5)
plot(costplot)

costplot

plot_pendulum(y)
%% DMP initialization
function dmp = dmp_init(n_dmps,bfs,dt,runtime,y,w,ay,by)

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
dmp.h = (diff(dmp.c)*0.55).^2;
dmp.h = 1./[dmp.h;dmp.h(end)];

dmp = dmp_reset(dmp);
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

%% DMP Set Goal
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

%% DMP Batch Fit
function [dmp,Y,Yd,Ydd] = dmp_batch_fit(dmp,traj,trajd,trajdd,flag)

%start state & goal state
y0 = traj(1,:);
goal = traj(end,:);
g = y0;

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
    gd = (goal-g).*dmp.ag;
    x = xd*dmp.dt + x;
    v = vd*dmp.dt + v;
    g = gd*dmp.dt + g;
end

%regression target
dmp.dG = goal - y0;
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
    dmp.w(a,:) = (dmp.sxtd(a,:)./(dmp.sx2(a,:)+1e-10));
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

%% DMP run
function [dmp,y,yd,ydd,weights] = dmp_run(dmp,w)
%weighted sum of locally weighted regression models
dmp.psi = exp(-0.5*((dmp.x-dmp.c).^2).*dmp.h);
amp = dmp.s;
in = dmp.x;

if nargin == 2
    dmp.w = w;
end
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

%% DMP run fit
function [dmp,Y,Yd,Ydd] = dmp_run_fit(dmp,flag,w)


if nargin < 3
    w = dmp.w;
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
WEIGHTS = W;

for i = 0:2*dmp.tau/dmp.dt
    [dmp,y,yd,ydd,weights] = dmp_run(dmp,w); % y is the actual trajectory of the end effector
    
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
    W(:,:,i+1) = dmp.w';
    WEIGHTS(:,:,i+1) = weights;
end


if flag == 1
    
    plot_pendulum(Y);
    
    % PLOTTING
    time = (0:dt:2*D.tau)';
    figure(4)
     clf(4)
     figure(4)
    for i = 1:n_dmps

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
    %     subplot(4,3,6)
    %     hold on
    %     for j = 1:D.bfs
    %         plot(time,squeeze(PSI(j,:,i))*squeeze(W(j,:,i))','*')
    %         title('{Psi activations/},{Basis functions}')
    %     end

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
        for j = 1:D.bfs
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
        plot(squeeze(WEIGHTS(end,i,:)))
        title('Weights')
        xlabel(sprintf('tau=%f',D.tau))
        
        drawnow
    end
end


end

%% Plot pendulum moving
function plot_pendulum(traj)

l=1; %length of pendulum
ys = traj(:,1) + l*sin(traj(:,2));
zs = -l*cos(traj(:,2));

figure(4)
 clf(4)
 figure(4)
for i = 1:length(traj)
    plot([traj(i,1) ys(i)],[0 zs(i)],'k')
    axis([-4 4 -4 4])
    drawnow
end
plot(ys,zs,'k')
axis([-4 4 -4 4])


end

%% Check the fit of the current DMP
% W
function cost = check_fit(dmp,traj,y,yd,ydd)

l=1;
ys = y(:,1) + l*sin(y(:,2));
zs = -l*cos(y(:,2));
end_effector = [ys zs];

trajdistancessq = inf*ones(size(y,1),1);
for i = 1:2*dmp.tau/dmp.dt
    trajdistancessq = min( min((ys(i)-traj(:,1)).^2 + (zs(i)-traj(:,2)).^2),...
    trajdistancessq);
end
yddsum = sum(rssq(ydd,2));

cost = 1000*sum(sqrt(trajdistancessq)) + 100*sum(sqrt((end_effector(end)-dmp.g).^2)) + ...
    sum(rssq(yd,2))*dmp.dt + yddsum*dmp.dt;

end



