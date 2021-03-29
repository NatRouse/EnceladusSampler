%%

clear all;
close all;
clc

g = 9.81; % gravity (m/s^2)

% max time, time step, and # of iterations. Using for loop instead of while loo
% makes preallocation easier
t_max = 100;
dt = .001;
iterations = t_max/dt+1;

% mass of car/gondola and point mass (sampler) (kg)
m_c = 50; 

sampler_radius = .05; %0.023127 corresponds to 1kg tungsten
Cd = .3;
A = 4*pi*sampler_radius^2;
rho = 19300; % density of tungsten
V = 4/3*pi*sampler_radius^3;
m_s = rho*V; 
%% initial robot configuration
l = 1; %(m)
theta = 0; %(rad)
R_y = -1; % (m)

R_y_dot = 0; % (m/s)
theta_dot = 0; % (rad/s)
l_dot = 0; % (m/s)

theta_ddot= 0; % rad/s^2 ASSUMPTION

% initial end effector position and velocities
y = R_y+l*sin(theta); % (m)
y_dot = R_y_dot+l_dot*sin(theta)+l*theta_dot*cos(theta); % (m/s) 

z = -l*cos(theta); % (m)
z_dot = -l_dot*cos(theta) + l*theta_dot*sin(theta); %(m/s)

%% histories for degrees of freedom and their derivatives
dof_hist = nan(iterations, 3);
dof_d_hist = nan(iterations, 3);
dof_dd_hist = nan(iterations, 3);

dof_hist(1,:) = [R_y theta l];
dof_d_hist(1,:) = [R_y_dot theta_dot l_dot];
dof_dd_hist(1,:) = [0 0 0];
% history for end effector position and derivatives
y_hist = nan(iterations,1);
y_hist(1) = [y];

y_d_hist = nan(iterations,1);
y_d_hist(1) = [y_dot];
y_dd_hist = nan(iterations,1);
z_hist = nan(iterations,1);
z_hist(1) = [z];
z_d_hist = nan(iterations,1);
z_d_hist(1) = [z_dot];
z_dd_hist = nan(iterations,1);

% history of forces to move car and increase tether length
F_Ry_hist = nan(length(dof_hist), 1);
F_l_hist = nan(length(dof_hist), 1);

% error history
error_hist = nan(iterations, 1);

%% constants for PD control of end effector position
%    K = 100.0*[0.099999999999990   0.039740446974352  -0.070597818316656   0.024439480872370  -0.000000033690876  -0.000000096623315;
%                 0.000000003369090  -0.000000004954016   0.000981083831427   0.000000091669307   1.000000000000002   1.000099995000491];
K_struct = load('K_mc50_rsampler5cm.mat');
K = K_struct.K;

k_Ry_y = K(1,1);
c_Ry_y = K(1,2);

k_theta_y = K(1,3);
c_theta_y = K(1,4);

k_l_y = K(1,5);
c_l_y = K(1,6);

k_Ry_z = K(2,1);
c_Ry_z = K(2,2);

k_theta_z = K(2,3);
c_theta_z = K(2,4);

k_l_z = K(2,5);
c_l_z = K(2,6);

%% Define map in terms of how "good" each position is 
sensor_radius = 2;
rho_air = 1.225;
D = 10;
dz = .1;
z_range = -10:dz:0;
dr = 0.1;
y_range = -D/2:dr: D/2;
mu = 1.81e-5;
Q = .1; % m^3/s
[Y1, Z1, u] = pipe_flow(rho_air, mu, Q, z_range, dr, 'circular', D);
U = repmat(u, [length(z_range), 1]);  
goodness = U.*(-Z1);
surf(Y1,Z1,goodness);
%%
% Create a grid of evenly spaced points in two-dimensional space.
% y1 = -2.5:0.05:2.5;
% z_range = -10:0.1:0;
% 
% [Y1,Z1] = meshgrid(y1,z_range);
% mesh = [Y1(:) Z1(:)];
% 
% goodness = 100*(-(Y1-.5).^2-Z1+1);
% %Evaluate the pdf of the normal distribution at the grid points.
% goodness = reshape(goodness,length(z_range),length(y1));
% 
% % Plot the pdf values.
% 
% surf(y1,z_range,goodness)
% caxis([min(goodness(:))-0.5*range(goodness(:)),max(goodness(:))])
% axis([min(y1) max(y1) min(z_range) max(z_range)])
% xlabel('y1')
% ylabel('z_range')
% zlabel('Probability Density')

%%
explored = zeros(size(goodness));
explored_hist = zeros(size(explored,1), size(explored,2), iterations);
goodness_hist = 0;
max_goodness = max(goodness_hist);

y_desired = y;
z_desired = z;
y_desired_hist = [];
z_desired_hist = [];

error_eqn = 1;
R_y_ddot = 0;
l_ddot = 0;
%%
i = 1;
% R_y_desired = 0;
% l_desired = 10;
theta_desired = 0;
while i < iterations && error_eqn > .2
    %% Define the mapped area
    is_in_range = ((Y1-y).^2+(Z1-z).^2)<=sensor_radius^2;
    explored = or(is_in_range, explored);
    explored_hist(:,:,i) = explored; % update history of mapped region
    
    good_in_range = goodness.*is_in_range; % Value of region in view
    desired_index = find(good_in_range == max(good_in_range, [], 'all'));
    if length(desired_index)>1
        desired_index = desired_index(end);       
    end
    if  goodness(desired_index) > max_goodness
        max_goodness = goodness(desired_index);
        y_desired = Y1(desired_index);
        z_desired = Z1(desired_index);
%         disp('new target: ')
        y_desired_hist = [y_desired_hist; y_desired];
        z_desired_hist = [z_desired_hist; z_desired];
        R_y_desired = y_desired;
        l_desired = -z_desired;
        goodness_hist = [goodness_hist;goodness(desired_index)];
    end    
    
    %% End effector desired accelerations
    mod_theta = sign(theta)*mod(abs(theta), 2*pi);
    if mod_theta>-pi/2 && mod_theta<pi/2
        s = 1;
    elseif mod_theta<-pi/2 || mod_theta>pi/2
        s = -1;
    end
    y_index = find(min(abs(y_range-y), [], 'all')==abs(y_range-y));
    z_index = find(min(abs(z_range-z), [], 'all')==abs(z_range-z));
    
    % find z_dot = l_dot*cos(theta)+l*theta_dot*sin(theta);
    % vz_rel = U-z_dot;
    %vy_rel = R_y_dot+l_dot*sin(theta)+l*theta_dot*cos(theta);
    %v_rel = sqrt(vz_rel^2+vy_rel^2);
    
    % Define location of sphere in the flow
%     y = R_y+l*sin(theta);
%     z = -l*cos(theta);
%     y_dot = R_y_dot + l_dot.*sin(theta) + l.*theta_dot.*cos(theta);
%     z_dot = -l_dot.*cos(theta)+l.*theta_dot.*sin(theta);
    
    y_dot_rel = -y_dot;
    z_dot_rel = u-z_dot;
    
    u_rel = (y_dot_rel^2+z_dot_rel.^2).^0.5;
    y_sampler_max = y_range(min(abs(y+sampler_radius-y_range))==abs(y+sampler_radius-y_range));
    y_sampler_min = y_range(min(abs(y-sampler_radius-y_range))==abs(y-sampler_radius-y_range));
    y_min_index = find(y_range == y_sampler_min);
    y_max_index = find(y_range == y_sampler_max);
    u_rel_sampler = u_rel(y_min_index:y_max_index);
    avg_y_index = round((y_max_index+y_min_index)/2);
    fluid_angle = atan2(z_dot_rel(avg_y_index),y_dot_rel);
    
    F_fluid = mean(1/2*rho_air*Cd*u_rel_sampler.^2*A);
    F_fluid_y = F_fluid*cos(fluid_angle);
    F_fluid_z = F_fluid*sin(fluid_angle);

    F_z = -k_Ry_z*(R_y-R_y_desired) - c_Ry_z*R_y_dot - k_theta_z*(theta-theta_desired) - c_theta_z*theta_dot - k_l_z*(l-l_desired) - c_l_z*l_dot + F_fluid_z;
%     y = R_y + l*sin(theta);
%     y_dot = R_y_dot + l_dot*sin(theta) + l*theta_dot*cos(theta);
%     F_y = -k_y*(y-y_desired) - c_y*y_dot;
    F_y = -k_Ry_y*(R_y-R_y_desired) - c_Ry_y*R_y_dot - k_theta_y*(theta-theta_desired) - c_theta_y*theta_dot - k_l_y*(l-l_desired) - c_l_y*l_dot + F_fluid_y;    
    F_l = F_z/cos(theta);
    F_Ry = F_y-F_l*sin(theta);


    %% accelerations of controlled DOFs needed to achieve desired end
    % effector accelerations
    l_ddot = -(R_y_ddot*sin(theta)+R_y_dot*theta_dot*cos(theta)-l*theta_dot^2-R_y_dot*theta_dot*cos(theta)-g*cos(theta))+F_l/m_s;
    R_y_ddot = (F_Ry-m_s*(l_ddot*sin(theta)+2*l_dot*theta_dot*cos(theta)+theta_ddot*l*cos(theta)-theta_dot^2*l*sin(theta)))/(m_c+m_s);
    % uncontrolled dof acceleration
    theta_ddot= -(2*l_dot*theta_dot+R_y_ddot*cos(theta)+g*sin(theta))/l;

    %% update degrees of freedom
    R_y = R_y + R_y_dot*dt + 1/2*R_y_ddot*dt^2;
    R_y_dot = R_y_dot + R_y_ddot*dt;
    
    l = l + l_dot*dt + 1/2*l_ddot*dt^2;
    l_dot = l_dot + l_ddot*dt;
    
    theta = theta + theta_dot*dt + 1/2*theta_ddot*dt^2;
    theta_dot = theta_dot + theta_ddot*dt;
    dof_hist(i,:) = [R_y theta l];
    dof_d_hist(i,:) = [R_y_dot theta_dot l_dot];
    dof_dd_hist(i,:) = [R_y_ddot theta_ddot l_ddot];
    
    %% update end effector coords and velocities
    y_dot = R_y_dot + l_dot.*sin(theta) + l.*theta_dot.*cos(theta);
    z_dot = -l_dot.*cos(theta)+l.*theta_dot.*sin(theta);
    
    y = R_y+l*sin(theta);
    z = -l*cos(theta);
    
    y_hist(i) = y;
    y_d_hist(i) = y_dot;
%     y_dd_hist(i) = y_ddot;
    
    z_hist(i) = z;
    z_d_hist(i) = z_dot;
%     z_dd_hist(i) = z_ddot;   

    % update error and error history
    error_eqn = sqrt((y_desired-y)^2 + (z_desired-z)^2)+ sqrt((l*theta_dot)^2+l_dot^2+R_y_dot^2);
    error_hist(i) = error_eqn;
        
    % Update force history
    F_Ry_hist(i) = F_Ry;
%     F_l_hist(i) = F_l;
    
    i = i+1;
end

 
%% Plot DOFs vs their respective time derivatives   
figure(2);
subplot(2, 2, 1) 
plot(dof_hist(:,1), dof_d_hist(:,1))
title('R_y_dot vs R_y')

subplot(2, 2, 2);
plot(dof_hist(:,2), dof_d_hist(:,2))
title('theta_dot vs theta')

subplot(2,2, 3);
plot(dof_hist(:,3), dof_d_hist(:,3))
title('l_dot vs l')


%% Plot configuration space vs time
figure(3)

subplot(2,2,1)
plot(dt:dt:dt*length(dof_hist), dof_hist(:,1))
title('Cart Position vs Time')
xlabel('Time (s)')
ylabel('Cart Position (m)')

hold on
subplot(2,2,2)
plot(dt:dt:dt*length(dof_hist), dof_hist(:,2))
title('Pendulum Angle vs Time')
xlabel('Time (s)')
ylabel('Pendulum Angle from Vertical (rad)')

subplot(2,2,3)
plot(dt:dt:dt*length(dof_hist), dof_hist(:,3))
title('Tether Length vs Time')
xlabel('Time (s)')
ylabel('Tether Length (m)')

 %% plot robot behavior
plot_bot = true;
if plot_bot
%         v = VideoWriter('Pendulum on Cart3');%
%         v.FrameRate = 10;%
%         open(v)%
    X1_reshaped = reshape(Y1, [1, numel(Y1)]);
    Z1_reshaped = reshape(Z1, [1, numel(Z1)]);

    for j = 1:50:i
    figure(4);
    goodness_reshaped = reshape(goodness.*explored_hist(:,:,j), [1, numel(goodness)]);
    goodness_reshaped(goodness_reshaped == 0) = nan;
    Ry_gd = dof_hist(j,1);
    theta_gd = dof_hist(j,2);
    l_gd = dof_hist(j,3);
    y = Ry_gd+l_gd*sin(theta_gd);
    z = -l_gd*cos(theta_gd);
    h = scatter(X1_reshaped, Z1_reshaped, [], goodness_reshaped, 'fill');
    hold on
    plot([Ry_gd; y],[0;z], 'r', 'Linewidth', 2)
    axis([min(y_range) max(y_range) min(z_range) 0])
    hold off
%         frame = getframe(gcf);
    drawnow
%         writeVideo(v, frame);
    end 
%         close(v)
end
