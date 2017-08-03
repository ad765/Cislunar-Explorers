function [ Xdot ] = dynamicsModel( ~, X, p )
% Dynamics equations for three-body problem. Three bodies involved are the
% spacecraft, moon, and sun. All positions and velocities are given in ECI
% frame. These are the continuous equations that will be used as the
% state-transition between the a posteriori current state and the a priori
% next state.
% 
% Inputs:
%           dt      - time steptime (unused, but don't delete)
%           z       - state vector for x at k-1 (18 element)
%           p       - parameters (mass; position of sun, moon, earth)
%
%
% Outputs:
%           zdot     - state vector for x at k (18 element)
%
% Anshuman Das, Cornell University
% Wednesday, August 2, 2018


% Cubesat positions and speeds
xc = X(1);
yc = X(2);
zc = X(3);
xcdot = X(4);
ycdot = X(5);
zcdot = X(6);

% Equations of motion
Xdot(1) = xcdot;
Xdot(2) = ycdot;
Xdot(3) = zcdot;
Xdot(4) = -p.muE*xc/(xc^2+yc^2+zc^2)^(3/2) - p.muM*((xc-p.xm)/((xc-p.xm)^2+(yc-p.ym)^2+(zc-p.zm)^2)^(3/2) + p.xm/(p.xm^2+p.ym^2+p.zm^2)^(3/2)) - p.muS*((xc-p.xs)/((xc-p.xs)^2+(yc-p.ys)^2+(zc-p.zs)^2)^(3/2) + p.xs/(p.xs^2+p.ys^2+p.zs^2)^(3/2));
Xdot(5) = -p.muE*yc/(xc^2+yc^2+zc^2)^(3/2) - p.muM*((yc-p.ym)/((xc-p.xm)^2+(yc-p.ym)^2+(zc-p.zm)^2)^(3/2) + p.ym/(p.xm^2+p.ym^2+p.zm^2)^(3/2)) - p.muS*((yc-p.ys)/((xc-p.xs)^2+(yc-p.ys)^2+(zc-p.zs)^2)^(3/2) + p.ys/(p.xs^2+p.ys^2+p.zs^2)^(3/2));
Xdot(6) = -p.muE*zc/(xc^2+yc^2+zc^2)^(3/2) - p.muM*((zc-p.zm)/((xc-p.xm)^2+(yc-p.ym)^2+(zc-p.zm)^2)^(3/2) + p.zm/(p.xm^2+p.ym^2+p.zm^2)^(3/2)) - p.muS*((zc-p.zs)/((xc-p.xs)^2+(yc-p.ys)^2+(zc-p.zs)^2)^(3/2) + p.zs/(p.xs^2+p.ys^2+p.zs^2)^(3/2));

Xdot = Xdot';

end