function [ Z ] = measurementModel( X, p, cur_T)
% Dynamics equations for three-body problem. Three bodies involved are the
% spacecraft, moon, and sun. All positions and velocities are given in ECI
% frame. These are the discretized equations that will be used as the
% state-transition matrix between k-1 and k.
%
% Inputs:
%           x       - state vector for x (6 element)
%           p       - parameters
%
%
% Outputs:
%           z       - measurement vector (5 element)
%
% Anshuman Das, Cornell University
% Wednesday, August 2, 2018


% Cubesat positions
xc = X(1);      yc = X(2);      zc = X(3);
xcdot = X(4);   ycdot = X(5);   zcdot = X(6);

% Ephemerides positions
p.xm = p.moon_x(cur_T);
p.ym = p.moon_y(cur_T);
p.zm = p.moon_z(cur_T);
p.xs = p.sun_x(cur_T);
p.ys = p.sun_y(cur_T);
p.zs = p.sun_z(cur_T);

%% Equations of motion (actual equations)
%
Z(1)    = (p.P/p.THETA)*acos((-xc*p.xm - yc*p.ym - zc*p.zm + xc^2 + yc^2 + zc^2)/(sqrt(xc^2+yc^2+zc^2)*((p.xm-xc)^2+(p.ym-yc)^2+(p.zm-zc)^2)));
Z(2)    = (p.P/p.THETA)*acos((-xc*p.xs - yc*p.ys - zc*p.zs + xc^2 + yc^2 + zc^2)/(sqrt(xc^2+yc^2+zc^2)*((p.xs-xc)^2+(p.ys-yc)^2+(p.zs-zc)^2)));
Z(3)    = (p.P/p.THETA)*acos((p.xm*(p.xs-xc) + p.ym*(p.ys-yc) + p.zm*(p.zs-zc) - xc*p.xs - yc*p.ys - zc*p.zm + xc^2 + yc^2 + zc^2)/(sqrt(xc^2+yc^2+zc^2)*((p.xm-xc)^2+(p.ym-yc)^2+(p.zm-zc)^2)));
Z(4)    = (2*p.P/p.THETA)*atan2(p.re,sqrt(xc^2 + yc^2 + zc^2));
Z(5)    = (2*p.P/p.THETA)*atan2(p.rm,sqrt((p.xm-xc)^2 + (p.ym-yc)^2 + (p.zm-zc)^2));
%}

%% Equations of motion (simplified measurement equations)
%{
Z(1)    = xc;
Z(2)    = yc;
Z(3)    = zc;
Z(4)    = xcdot;
Z(5)    = ycdot;
Z(6)    = zcdot;
%}
Z = Z';

end