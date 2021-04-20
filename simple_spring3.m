%% unforced spring
close all;
clear all;
clc

m_c = 50; 
sampler_radius = .05; %0.023127 corresponds to 1kg tungsten (m)
rho = 19300; % density of tungsten (kg/m^3)
V = 4/3*pi*sampler_radius^3; % sampler volume (m^3)
m_s = rho*V; % sampler mass (kg) 
g = 9.81;
dt = 0.001;

F_l = 0;
F_l_eq = 0; 

l_eq = 1;

l = l_eq;
l_dot = 0;
l_ddot = 0;

theta = pi/4;
theta_dot = 0;

y_c = 0;
y_c_dot = 0;
y_c_ddot = 0;

theta_hist = [theta];
theta_dot_hist = [theta_dot];

l_hist = [l];
l_eq_hist = [l_eq];

y_c_hist = [y_c];
y_c_dot_hist = [y_c_dot];


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

z_desired = -2;
y_desired = 0;    
l_desired = -z_desired;
y_c_desired = y_desired;
theta_desired = 0;

i = 1;
t_max = 10;

while i<t_max/dt+1

    F_z = -k_Ry_z*(y_c-y_c_desired) - c_Ry_z*y_c_dot - k_theta_z*(theta-theta_desired) - c_theta_z*theta_dot - k_l_z*(l-l_desired) - c_l_z*l_dot;
    F_l_eq = F_z;%/cos(theta);
%     F_l = F_l_eq;
    F_y = -k_Ry_y*(y_c-y_c_desired) - c_Ry_y*y_c_dot - k_theta_y*(theta-theta_desired) - c_theta_y*theta_dot - k_l_y*(l-l_desired) - c_l_y*l_dot;% + F_fluid_y;    
%     M = 0;

    if l<l_eq
        k = 0;
        F_l = 0;
        F_y_c = 0;
    else
        k = 10000;
        if F_l_eq < 0
            F_l = F_l_eq;
            F_y_c = F_y-F_l*sin(theta);

        else 
           F_l = 0;
           F_y_c = 0;
        end
    end


    theta_ddot = -(2*l_dot*theta_dot+y_c_ddot*cos(theta)+g*sin(theta))/l ;
    y_c_ddot = (F_y_c - m_s*(l_ddot*sin(theta)+2*l_dot*theta_dot*cos(theta)-l*theta_ddot*cos(theta)-l*theta_dot^2*sin(theta)))/(m_s+m_c);
    l_ddot = (F_l-k*l+k*l_eq)/m_s+g*cos(theta)+l*theta_dot^2-y_c_ddot*sin(theta);
    %l = 1/2*l_ddot*dt^2+l_dot*dt+l;
    %l_dot = l_ddot*dt+l_dot;

    theta = 1/2*theta_ddot*dt^2+theta_dot*dt+theta;
    theta_dot = theta_ddot*dt+theta_dot;
    
    l = 1/2*l_ddot*dt^2 + l_dot*dt + l;
    l_dot = l_ddot*dt + l_dot;

    y_c = 1/2*y_c_ddot*dt^2+y_c_dot*dt+y_c;
    y_c_dot = y_c_ddot*dt+y_c_dot;
    
    l_eq = 2; %l+F_l_eq/m_s; %2

    l_hist = [l_hist; l];
    l_eq_hist = [l_eq_hist; l_eq];
    
    theta_hist = [theta_hist; theta];
    theta_dot_hist = [theta_dot_hist; theta_dot];

    y_c_hist = [y_c_hist; y_c];
    y_c_dot_hist = [y_c_dot_hist; y_c_dot];
    
    if mod(i-1, 20)==0
        plot([y_c l*sin(theta)], [0 -l*cos(theta)])
        axis ([-4 4 -4 4])
        drawnow

    end
    i = i+1;
end

plot(l_hist)
hold on 
plot(l_eq_hist)
figure; 
plot(y_c_hist)
