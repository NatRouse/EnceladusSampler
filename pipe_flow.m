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
       %  solve the friction factor
      e_abs = 0.002; % lower bound for surface roughness of sea ice (m) 
       %Global Trends of Sea Ice: Small-Scale Roughness and Refractive Index,
       %2010
       syms f
       e_rel = e_abs/D;
       f = vpa(solve(1/sqrt(f) == -2.0*log10(e_rel/3.7 + 2.51/(Re*sqrt(f))), f));
    %    del_p = 
        u = [];
    end
    
    U = u;
    [Z, Y] = ndgrid(z_range, r);
%     U = repmat(u, [length(z_range), 1]);    
end
