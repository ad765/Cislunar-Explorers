% Unscented Kalman Filter for Trajectory Determination
%
% The following code simulations the position determination of the lunar
% cubesat given initial conditions and parameters for the 


clear,clc

%% Unpack and spline ephermerides data
% All in m and m/s
moon    = csvread('moon_30min.csv',1,1)*1000;
sun     = csvread('sun_30min.csv',1,1)*1000;

delta_t = 1/1800;
time    = 0:1:length(moon(:,1))-1;            % In intervals of 30 min
interptime = 0:delta_t:length(moon(:,1))-1;   % In intervals of 1 sec

p.moon_x  = interp1(time,moon(:,1),interptime,'spline');
p.moon_y  = interp1(time,moon(:,2),interptime,'spline');
p.moon_z  = interp1(time,moon(:,3),interptime,'spline');

p.sun_x  = interp1(time,sun(:,1),interptime,'spline');
p.sun_y  = interp1(time,sun(:,2),interptime,'spline');
p.sun_z  = interp1(time,sun(:,3),interptime,'spline');


%% Parameters and Initializations

% Time array
t0      = 0;                % start time (in seconds)
secpmin = 60;               % seconds per minute
minutes = 10;               % end time (in minutes)
tf      = minutes*secpmin;  % end time (in seconds)
dt_dyn  = 1;                % dynamics update time interval (in seconds)
dt_fil  = 10;               % filter update time interval (in seconds)

% State and measurement properties
L       = 6;                % length of state
M       = 5;                % length of measurements

% Planetary and camera parameters
G       = 6.67e-11;         % gravitational constant
ms      = 1.98855e30;       % mass of sun
me      = 5.9736e24;        % mass of earth
mm      = 7.348e22;         % mass of moon
p.rs    = 6.957e8;          % radius of sun
p.re  	= 6.371e6;          % radius of earth
p.rm  	= 1.737e6;          % radius of moon
p.P  	= 2464;             % pixels resolution (Raspberry PI Camera v2)
p.THETA = 0.9337;           % field of view     (CAMERA MODULE datasheet)
p.muM   = G*mm;             % std gravitational parameter of moon
p.muE   = G*me;             % std gravitational parameter of earth
p.muS   = G*ms;             % std gravitational parameter of sun
p.q     = 1e-5;             % std deviation of process noise
p.r     = 0.12*(180/pi)*p.P/p.THETA; % std deviation of measurement noise
Q_k     = p.q^2*eye(L);     % covariance of process noise
R_k     = p.r^2*eye(M);     % covariance of measurement noise

% Initial state
xc0     = 408e3+p.re;
yc0     = 0;
zc0     = 0;
xcdot0  = 0;
ycdot0  = 7.67e3;
zcdot0  = 0;
X0      = [xc0; yc0; zc0; xcdot0; ycdot0; zcdot0];

%% Dynamics Model Check
X_state = [];

for i = t0+1:dt_dyn:tf+1
    X_temp = stateTransition( dt_dyn, X0, p, i);
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
P_kkm1 = 1e6*eye(L);

% UKF parameters
alpha   = 1e-3;                     % default, tunable
kappa   = 0;                        % default, tunable
beta    = 2;                        % default, tunable
lambda  = alpha^2*(L+kappa)-L;
Wm      = [lambda/(L+lambda);
    (1/(2*(lambda+L)))*ones(2*L,1)];
Wc      = [lambda/(L+lambda) + (1-alpha^2+beta);
    (1/(2*(lambda+L)))*ones(2*L,1)];

THEORY      = [];
STATE       = [];
MEASURED    = [];
P11         = [];
P22         = [];
P33         = [];
INNOVATION  = [];

for i = t0+1:dt_fil:tf+1
    
    P11     = [P11, 3*sqrt(P_kkm1(1,1))];
    P22     = [P22, 3*sqrt(P_kkm1(2,2))];
    P33     = [P33, 3*sqrt(P_kkm1(3,3))];
    
    % Simulated measurement
    Y_meas  = measurementModel( X_state(:,i), p, i) + normrnd(0,p.r,[M,1]);
    THEORY = [THEORY, X_state(:,i)];
    MEASURED = [MEASURED, Y_meas];
    
    % Sigma points of state estimate (k-1)
    Xsp_kkm1 = sigmaPoints( X_kkm1, real(P_kkm1), lambda+L);
    lsp = length(Xsp_kkm1);
    
    % Sigma points of state measurement
    Ysp_kkm1 = zeros(M,lsp);
    for j = 1:lsp
        Ysp_kkm1(:,j) = measurementModel( Xsp_kkm1(:,j), p, i);
    end
    
    % State measurement
    Y_k = zeros(M,1);
    for j = 1:lsp
        Y_k = Y_k + Wm(j)*Ysp_kkm1(:,j);
    end
    
    % State measurement covariance
    P_YY = zeros(M,M);
    for j = 1:lsp
        P_YY = P_YY + Wc(j)*(Ysp_kkm1(:,j) - Y_k)*(Ysp_kkm1(:,j) - Y_k)';
    end
    P_YY = P_YY + R_k;
    
    % State estimation-measurement cross-covariance
    P_XY = zeros(L,M);
    for j = 2:lsp
        P_XY = P_XY + (Xsp_kkm1(:,j) - X_kkm1)*(Ysp_kkm1(:,j) - Y_k)';
    end
    P_XY = (1/(2*lambda))*P_XY;
    
    % Kalman gain
    K = P_XY/P_YY;
    
    % State estimation and state estimation error covariance
    INNOVATION = [INNOVATION, (Y_meas - Y_k)];
    X_kk = X_kkm1 + K*(Y_meas - Y_k);
    P_kk = (P_kkm1 - K*P_YY*K');
    
    % Sigma points of state estimate (k)
    Xsp_kk = sigmaPoints( X_kk, P_kk, lambda+L);
    
    % State-transition model
    Xsp_kp1k = zeros(size(Xsp_kkm1));
    for j = 1:lsp
        Xsp_kp1k(:,j) = stateTransition( dt_fil, Xsp_kk(:,j), p, i);
    end
    
    % State estimation
    X_kp1k = zeros(L,1);
    for j = 1:lsp
        X_kp1k = X_kp1k + Wm(j)*Xsp_kp1k(:,j);
    end
    
    % Covariance estimation
    P_kp1k  = zeros(L,L);
    for j = 1:lsp
        P_kp1k = Wc(j)*(Xsp_kp1k(:,j) - X_kp1k)*(Xsp_kp1k(:,j) - X_kp1k)';
    end
    P_kp1k  = P_kp1k + Q_k;
    
    % Reinitialize
    X_kkm1  = X_kp1k;
    P_kkm1  = P_kp1k;
    
    STATE   = [STATE, X_kp1k];
    
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
plot(x_est-x_theory)
title('Position (x) Error vs. Time')
legend('Estimate','Theory')
hold off

figure(3)
hold on
plot(P22)
plot(-P22)
plot(y_est-y_theory)
title('Position (y) Error vs. Time')
legend('Estimate','Theory')
hold off

figure(4)
hold on
plot(P33)
plot(-P33)
plot(z_est-z_theory)
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
