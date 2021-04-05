% rho_air = 1.225;
% D = 1; % pipe diameter (m)
% dz = .1;
% z_range = -10:dz:0;
% dr = 0.1;
% y_range = -D/2:dr: D/2;
% mu = 1.81e-5;
% Q = 10; % m^3/s
% [Y1, Z1, u] = pipe_flow(rho_air, mu, Q, z_range, dr, 'circular', D);
% U = repmat(u, [length(z_range), 1]);
% U0_max = max(U, [], 'all');
% goodness = U.*(-Z1);
% surf(Y1,Z1,goodness);
function [Y, Z, U] = pipe_flow(rho, mu, Q, z_range, dr, varargin)
    % Calculates the fluid flow behavior inside a pipe. Valid for viscous,
    % subsonic laminar flow. 

    % rho =	density of the fluid, kg/m^3
    % u	= flow speed, m/s
    % L	= characteristic linear dimension, m
    % m = dynamic viscosity of the fluid, 
    % D = pipe diameter, m
    % L_pipe = pipe length, m
    % theta = angle of pipe wrt to the horizontal
    % calculate pressure drop throughout the pipe
    % del_p = Q*128*mu*L_pipe/(pi*D^4) + rho*g*sin(theta); 

    % Calculate Reynold's number 
    if strcmp(varargin{1}, 'circular') 
       D = varargin{2};
       Re_crit = 2100;
       A = pi*(D/2)^2;
       u_avg = Q/A; 
    elseif strcmp(varargin{1}, 'rectangular')
       a = varargin{2}(1);
       b = varargin{2}(2);
       A= a*b;
       P = 2*(a+b);
       D = 4*A/P;
       u_avg = Q/A;
       Re_crit = 2800; % Assuming high aspect ratio, the critical Reynold's 
       % number approaches the value for a flat plate (Critical Reynolds Number
       % for Newtonian Flow in Rectangular Ducts, 1988)

    end



    Re = rho*u_avg*D/mu;
    r = -D/2:dr:D/2;
    R = D/2;

    if Re<Re_crit % laminar flow
        V_c = 2*Q/(pi*R^2);
        u = V_c*(1-(2*r./D).^2); % Velocity profile
    elseif Re>Re_crit % turbulent flow
        
%         L = max(abs(z_range));
%         f = (1/(-1.8*log10((epsilon/(3.7*D))^1.11+6.9/Re)))^2; % friction factor from modified Colebrook formula
%         del_p = Q*128*mu*L/D^4; % Eqn 8.9 in Fundamentals of Fluid Mechanics (fully developed laminar flow in pipes)

        
        n = 7; % power law exponent depends on Reynolds number. It's ~7 for many flows 
        V_c = (n+1)*(2*n+1)/(2*n^2)*u_avg;
        u = V_c*(1-abs(r./R)).^(1/n);
%         del_z = range(z_range);
%         
%         del_p = rho*f*L*v^2/(2*D)-rho*g*del_z;
%         
%         tau_w = del_p*D/(4*L); % Eqn 8.5
%         u_star = (tau_w/rho)^.5;
%         
%         % viscuous sublayer
%         nu = mu/rho; % kinematic viscosity
%         y_vs = R-r;
%         u_bar_v = y_vs*u_star^2/nu; 
%         
%         % overlap sublayer
%         u_bar_o = u_star*2.5*ln(y_ts*u_star/nu)+
%         
%         % Turbulent sublayer
%         u_bar_t = Vc*(1-r/R)^1/n
%         
        
       %  solve the friction factor
%       e_abs = 0.002; % lower bound for surface roughness of sea ice (m) 
%        %Global Trends of Sea Ice: Small-Scale Roughness and Refractive Index,
%        %2010
%        syms f
%        e_rel = e_abs/D;
%        f = vpa(solve(1/sqrt(f) == -2.0*log10(e_rel/3.7 + 2.51/(Re*sqrt(f))), f));
%     %    del_p = 
%         u = [];
    end
    
    U = u;
    [Z, Y] = ndgrid(z_range, r);
%     U = repmat(u, [length(z_range), 1]);    
end
