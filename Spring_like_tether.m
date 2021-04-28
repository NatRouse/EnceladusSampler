%% fixed length tether 
close all;
clear all;
clc
plot_robot = true;
%% Environmental Parameters
g = 9.81; % gravity (m/s^2)

% max time, time step, and # of iterations. Using for loop instead of while
% loop
% makes preallocation easier
t_max = 10;
dt = .001;
iterations = t_max/dt+1;

% vent geometry
D = 1; % pipe diameter (m)
dz = .1; % used for defining the vent
z_range = -10:dz:0;
dr = 0.01;
y_range = -D/2:dr: D/2;

% Fluid parameters
rho_air = 1.225; % density (kg/m^3)
mu = 1.81e-5; % dynamic viscosity
Q =100; % volumetric flowrate (m^3/s)

% specify the nominal fluid flow
[Y1, Z1, u] = pipe_flow2(rho_air, mu, Q, z_range, dr, 'circular', D);

Y1_reshaped = reshape(Y1, [1, numel(Y1)]); % convert meshgrid to vector form
Z1_reshaped = reshape(Z1, [1, numel(Z1)]);

U0 = repmat(u, [length(z_range), 1]); % velocity 
I = 0.05; % turbulence intensity
u_prime_avg = I*U0; % standard deviation of the turbulence

% How desireable each point in the flow is
goodness = U0.*(-Z1);
surf(Y1,Z1,goodness);

%% Robot physical parameters
sensor_radius = 1; % sampler can detect flowrate up this distance away (m)
% mass of cart/gondola and point mass (sampler) (kg)
m_c = 50; 
sampler_radius = .05; %0.023127 corresponds to 1kg tungsten (m)
Cd = .3; % approximate drag coefficient of sphere in turbulent flow (assuming continuum mechanics)
A = 4*pi*sampler_radius^2; % frontal surface area of sampler (m^2)
rho = 19300; % density of tungsten (kg/m^3)
V = 4/3*pi*sampler_radius^3; % sampler volume (m^3)
m_s = rho*V; % sampler mass (kg)


%% Initial conditions and history
F_l = 0;
F_l_eq = 0; 

l_unstretched = 4;

l = l_unstretched;
l_dot = 0;
l_ddot = 0;

theta = 0;
theta_dot = 0;

y_c = 0;
y_c_dot = 0;
y_c_ddot = 0;

theta_hist = [theta];
theta_dot_hist = [theta_dot];

l_hist = [l];
l_unstretched_hist = [l_unstretched];

y_c_hist = [y_c];
y_c_dot_hist = [y_c_dot];

F_l_sampler_hist = [];
%% Robot control parameters
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

k_l_z = 100;%K(2,5);
c_l_z = K(2,6)/10;

%% Desired End Effector Coordinates

z_desired = -8;
y_desired = 0.25;    
l_desired = -z_desired;
y_c_desired = y_desired;
theta_desired = 0;

%% Force histories
F_l_control_hist = [];
F_z_control_hist = [];
F_y_control_hist = [];

F_l_sampler_control_hist = [];
F_y_sampler_control_hist = [];
F_z_sampler_control_hist = [];


i = 1;

while i<t_max/dt+1
    
    
    %% end effector position
    y = y_c + l*sin(theta);
    y_dot = y_c_dot + l_dot*sin(theta) + l*theta_dot*cos(theta);
    z = -l*cos(theta);
    z_dot = -l_dot*cos(theta) + l*theta_dot*sin(theta);
    
    if z>0
        disp("The sampler is outside the vent!")
       break
    end
    
%% Identify fluid forces
    y_index = find(min(abs(y_range-y), [], 'all')==abs(y_range-y));
    z_index = find(min(abs(z_range-z), [], 'all')==abs(z_range-z));
    
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
%     disp(F_fluid)
    F_fluid_y = F_fluid*cos(fluid_angle);
    F_fluid_z = F_fluid*sin(fluid_angle);
%% Calculate control forces
    F_z_control = -k_Ry_z*(y_c-y_c_desired) - c_Ry_z*y_c_dot - k_theta_z*(theta-theta_desired) - c_theta_z*theta_dot + k_l_z*(l_unstretched-l_desired)*cos(theta) + c_l_z*l_dot*cos(theta);
    F_z_control_hist = [F_z_control_hist; F_z_control];
    
    F_l_control = F_z_control/cos(theta); % 
    F_l_control_hist = [F_l_control_hist; F_l_control];

%     F_l = F_l_eq;
    F_y_cart_control = -k_Ry_y*(y_c-y_c_desired) - c_Ry_y*y_c_dot - k_theta_y*(theta-theta_desired) - c_theta_y*theta_dot + k_l_y*(l_unstretched-l_desired)*sin(theta) + c_l_y*l_dot*sin(theta);    
%     M = 0;
%% Calculate forces experienced by each body
    % if the tether is not in tension
    if l<l_unstretched
        k_l_z = 0;% 0.001; % tether does not resist compression
        % Neither pushing out tether or reeling in tether will affect the
        % sampler.
        F_z_sampler_control = 0; 
        F_l_sampler_control = 0;
        F_y_sampler_control = 0; 
        color = 'r';
        
    else % if the tether is in tension
        k_l_z = 10000; % tether resists tension
        % unraveling more tether will not move the sampler (but it will
        % change the tether length, which will impact the tension on the
        % next time step)
        if F_l_control < 0
%             F_z_sampler_control = 0;
%             F_y_sampler_control = 0;
            F_l_sampler_control = 0; %F_l_control;
            F_z_sampler_control = 0;
            F_y_sampler_control = 0;
            
        % Reeling in tether will apply force directly to the sampler if it is already fully extended    
        elseif F_l_control >= 0
            F_z_sampler_control = F_z_control;
            F_l_sampler_control = F_l_control;
            F_y_sampler_control = -F_l_sampler_control*sin(theta);
        
        end
        
        color = 'b';
%         else 
%            F_l = 0;
%            F_y_c = 0;
%         end
    end
    
    F_y_sampler_net = F_y_sampler_control + F_fluid_y;
    F_z_sampler_net = F_z_sampler_control + F_fluid_z;
    
    F_z_cart = 0; % no net force on the cart in vertical direction because it's on a rigid rail
%     F_y_cart_net = F_y_cart_control + F_y_sampler_control + F_fluid_y; % net force on cart in horizontal depends on control and tension in the tether
    
M_theta = 0 + 0 + F_fluid_y*l*cos(theta) + F_fluid_z*l*sin(theta);
F_y_cart = F_y_cart_control + F_fluid_y - F_y_sampler_control;
F_l = 0-F_l_sampler_control-F_fluid_y*sin(theta)+F_fluid_z*cos(theta); % Upward flow of fluid decreases tether length


    theta_ddot = -(2*l_dot*theta_dot+y_c_ddot*cos(theta)+g*sin(theta))/l + M_theta/l^2;
    y_c_ddot = (F_y_cart - m_s*(l_ddot*sin(theta)+2*l_dot*theta_dot*cos(theta)-l*theta_ddot*cos(theta)-l*theta_dot^2*sin(theta)))/(m_s+m_c);
    l_ddot = (F_l-k_l_z*l+k_l_z*l_unstretched)/m_s+g*cos(theta)+l*theta_dot^2-y_c_ddot*sin(theta);

    theta = 1/2*theta_ddot*dt^2+theta_dot*dt+theta;
    theta_dot = theta_ddot*dt+theta_dot;
    
    l = 1/2*l_ddot*dt^2 + l_dot*dt + l;
    l_dot = l_ddot*dt + l_dot;

    y_c = 1/2*y_c_ddot*dt^2+y_c_dot*dt+y_c;
    y_c_dot = y_c_ddot*dt+y_c_dot;
    

%         l_unstretched = 5;
%     if F_l_control<0   % Only increase tether length if controller applied a force releasing tether. Do not reel in tether
%         l_unstretched = l-F_l_control/k_l_z; %2
%     else
%         F_l_unstretched = 0;
% %         l_unstretched = l-F_l_unstretched/k_l_z;
%     end 

    
    if k_l_z>0   % only adjust tether if l>l_unstretched
        l_unstretched = l-F_l_control/k_l_z; %2
    else
        F_l_unstretched = 0;
%         l_unstretched = l-F_l_unstretched/k_l_z;
    end 
    
    l_hist = [l_hist; l];
    l_unstretched_hist = [l_unstretched_hist; l_unstretched];
    
    theta_hist = [theta_hist; theta];
    theta_dot_hist = [theta_dot_hist; theta_dot];

    y_c_hist = [y_c_hist; y_c];
    y_c_dot_hist = [y_c_dot_hist; y_c_dot];
    
    if mod(i-1, 20)==0 && plot_robot
        plot([y_c y_c+l*sin(theta)], [0 -l*cos(theta)], color)
        axis ([-4 4 -10 1])
        drawnow

    end
    F_l_sampler_hist = [F_l_sampler_hist;F_l_sampler_control];
    
    i = i+1;
end

%%
figure(2);
plot(l_hist)
hold on 

plot(l_unstretched_hist)

%%
figure(3); 
plot(y_c_hist)
%%
figure(4)
plot(F_l_sampler_hist(10:end),'.')
hold on
plot(F_l_control_hist(10:end),'-')
