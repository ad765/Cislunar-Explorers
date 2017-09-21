% Unscented Kalman Filter for Trajectory Determination
%
% The following code simulations the position determination of the lunar
% cubesat given initial conditions and parameters for the
%
% Use https://ssd.jpl.nasa.gov/horizons.cgi to find correct ephemeris.
% Directions for ephemeris tables are given in README.md.

clear,clc

%% PROGRAM INPUTS
moonTable   = 'moon_eph.txt';   % Moon ephemeris data sheet
sunTable    = 'sun_eph.txt';    % Sun ephemeris data sheet
dt_fil      = 10;               % filter update time interval (in seconds)
minutes     = 10;               % end time (in minutes)
p.L         = 6;                % length of state
p.M         = 5;                % length of measurements
xc0     = 408e3+6378.1e3;           
yc0     = 0;
zc0     = 0;
xcdot0  = 0;
ycdot0  = 7.67e3;
zcdot0  = 0;


%% Unpack and spline ephermerides data
% All in meters and in intervals of seconds
[moonX, moonY,  moonZ,  ~,~,~]  = txt2csv(moonTable);
[sunX,  sunY,   sunZ,   ~,~,~]  = txt2csv(sunTable);

delta_t = 1/60;
time    = 0:1:length(moonX)-1;            % In intervals of 1 min
interptime = 0:delta_t:length(moonX)-1;   % In intervals of 1 sec

p.moon_x  = interp1(time,moonX,interptime,'spline')*1000;
p.moon_y  = interp1(time,moonY,interptime,'spline')*1000;
p.moon_z  = interp1(time,moonZ,interptime,'spline')*1000;

p.sun_x  = interp1(time,sunX,interptime,'spline')*1000;
p.sun_y  = interp1(time,sunY,interptime,'spline')*1000;
p.sun_z  = interp1(time,sunZ,interptime,'spline')*1000;


%% Parameters and Initializations

% Time array
t0      = 0;                % start time (in seconds)
secpmin = 60;               % seconds per minute
tf      = minutes*secpmin;  % end time (in seconds)
dt_dyn  = 1;                % dynamics update time interval (in seconds)

% Planetary and camera parameters
G       = 6.67e-11;         % gravitational constant
ms      = 1.98855e30;       % mass of sun
me      = 5.9736e24;        % mass of earth
mm      = 7.348e22;         % mass of moon
p.rs    = 6.957e8;          % radius of sun
p.re  	= 6.371e6;          % radius of earth
p.rm  	= 1.737e6;          % radius of moon
p.P  	= 2464;             % pixels resolution (Raspberry PI Camera v2)
p.THETA = 52.40*pi/180;     % field of view     (CAMERA MODULE datasheet)
p.muM   = G*mm;             % std gravitational parameter of moon
p.muE   = G*me;             % std gravitational parameter of earth
p.muS   = G*ms;             % std gravitational parameter of sun
p.q     = 1e-4;             % std deviation of process noise
p.r     = 0.12*(180/pi)*p.P/p.THETA; % std deviation of measurement noise
p.Q_k   = p.q^2*eye(p.L);     % covariance of process noise
p.R_k   = p.r^2*eye(p.M);     % covariance of measurement noise

% Initial state
X0      = [xc0; yc0; zc0; xcdot0; ycdot0; zcdot0];

%% Dynamics Model Check (Simulation)
X_state = [];

for k = t0+1:dt_dyn:tf+1
    X_temp = stateTransition( @dynamicsModel, dt_dyn, X0, p, k);
    X_state = [X_state, X_temp];
    X0 = X_temp;
end

%% Measurement Model Check
%{
X_state_noisy = zeros(size(X_state));
for i = 1:length(X_state)
    X_state_noisy(:,i) = measurementModel(X_state(:,i),p);
end

x_noisy = X_state_noisy(1,:)';
y_noisy = X_state_noisy(2,:)';
z_noisy = X_state_noisy(3,:)';

disp(mean(x_noisy-x))
disp(mean(y_noisy-y))
disp(mean(z_noisy-z))

disp(std(x_noisy-x))
disp(std(y_noisy-y))
disp(std(z_noisy-z))

figure(3)
hold on
subplot(3,1,1)
plot(t,x_noisy-x,'k')
subplot(3,1,2)
plot(t,y_noisy-y,'k')
subplot(3,1,3)
plot(t,z_noisy-z,'k')
hold off
%}

%% Unscented Kalman Filter

% Initial state estimate
X_kkm1  = [xc0; yc0; zc0; xcdot0; ycdot0; zcdot0];

% Initial state covariance
P_kkm1 = 1e6*eye(p.L);

% UKF parameters
p.alpha   = 1e-3;                     % default, tunable
p.kappa   = 0;                        % default, tunable
p.beta    = 2;                        % default, tunable
p.lambda  = p.alpha^2*(p.L+p.kappa)-p.L;
p.Wm      = [p.lambda/(p.L+p.lambda);
    (1/(2*(p.lambda+p.L)))*ones(2*p.L,1)];
p.Wc      = [p.lambda/(p.L+p.lambda) + (1-p.alpha^2+p.beta);
    (1/(2*(p.lambda+p.L)))*ones(2*p.L,1)];

THEORY      = [];
STATE       = [];
MEASURED    = [];
P11         = [];
P22         = [];
P33         = [];
INNOVATION  = [];

for k = t0+1:dt_fil:tf+1
    
    % Simulated measurement
    Y_meas      = measurementModel( X_state(:,k), p, k) + normrnd(0,p.r,[p.M,1]);
    THEORY      = [THEORY, X_state(:,k)];
    MEASURED    = [MEASURED, Y_meas];
    
    % Unscented Kalman filter
    [ X_kp1k, P_kp1k, diff ] = UKF( @dynamicsModel, ...
        measurementModel( X_state(:,k), p, k),...
        X_kkm1, P_kkm1, Y_meas, p, k, dt_fil);
    
    % Reinitialize
    X_kkm1  = X_kp1k;
    P_kkm1  = P_kp1k;
    
    % Tabulate arrays
    STATE       = [STATE, X_kkm1];
    INNOVATION  = [INNOVATION, diff];
    P11         = [P11, 3*sqrt(P_kkm1(1,1))];
    P22         = [P22, 3*sqrt(P_kkm1(2,2))];
    P33         = [P33, 3*sqrt(P_kkm1(3,3))];
    
end

%% Data Assimilation
% Theoretical
x_theory = THEORY(1,:);
y_theory = THEORY(2,:);
z_theory = THEORY(3,:);

% Estimates
x_est = STATE(1,:);
y_est = STATE(2,:);
z_est = STATE(3,:);

% Measurements
x_mea = MEASURED(1,:);
y_mea = MEASURED(2,:);
z_mea = MEASURED(3,:);

%% Graphing and Simulation
figure(1)
hold on
axis equal
plot3(x_est,y_est,z_est)
plot3(x_theory,y_theory,z_theory)
hold off

figure(2)
hold on
plot(P11)
plot(-P11)
%plot(x_est-x_theory)
title('Position (x) Error vs. Time')
legend('Estimate','Theory')
hold off

figure(3)
hold on
plot(P22)
plot(-P22)
%plot(y_est-y_theory)
title('Position (y) Error vs. Time')
legend('Estimate','Theory')
hold off

figure(4)
hold on
plot(P33)
plot(-P33)
%plot(z_est-z_theory)
title('Position (z) Error vs. Time')
legend('Estimate','Theory')
hold off

figure(5)
hold on
title('Cubesat-Moon Angle Innovation')
plot(INNOVATION(1,:))
hold off

figure(6)
hold on
title('Cubesat-Sun Angle Innovation')
plot(INNOVATION(2,:))
hold off

figure(7)
hold on
title('Moon-Sun Angle Innovation')
plot(INNOVATION(3,:))
hold off

figure(8)
hold on
title('Earth Diameter Innovation')
plot(INNOVATION(4,:))
hold off

figure(9)
hold on
title('Moon Diameter Innovation')
plot(INNOVATION(5,:))
hold off
