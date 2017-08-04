function [ Xout ] = stateTransition( dynamicsfcn, dt, Xin, p, cur_T )
% State transition model that propagates the state vector (X) from current
% time (cur_T) to the future time (cur_T+dt).
%
% Inputs:
%           dt      - time interval between current time and future time
%           Xin     - input state vector (at current time)
%           p       - parameter struct
%           cur_T   - current time
%
% Outputs:
%           Xout    - output state vector (at future time)
%
% Anshuman Das, Cornell University
% Wednesday, August 2, 2018


tol     = 1e-6;
options = odeset('AbsTol',tol,'RelTol',tol);
t0      = 0;


for i = 1:1:dt
    
    p.xm = p.moon_x(i+cur_T);
    p.ym = p.moon_y(i+cur_T);
    p.zm = p.moon_z(i+cur_T);
    p.xs = p.sun_x(i+cur_T);
    p.ys = p.sun_y(i+cur_T);
    p.zs = p.sun_z(i+cur_T);
    
    [~,Xout] = ode45( dynamicsfcn, [t0,t0+1], Xin, options, p);
    
    t0 = t0+1;
    Xin = Xout(end,:)';
    
end

Xout = Xin;

end

