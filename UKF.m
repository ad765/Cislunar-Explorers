function [ X_kp1k, P_kp1k, INNOVATION ] = UKF( dynamicsFCN, measurementFCN,...
    X_kkm1, P_kkm1, Y_meas, p, T_CURRENT, dT_FILTER)
% Unscented Kalman filter algorithm. Uses a proper weightage of the
% estimate and dynamics model 




% Sigma points of state estimate (k-1)
Xsp_kkm1 = sigmaPoints( X_kkm1, P_kkm1, p.lambda+p.L);
lsp = size(Xsp_kkm1,2);

% Sigma points of state measurement
Ysp_kkm1 = zeros(p.M,lsp);
for i = 1:lsp
    Ysp_kkm1(:,i) = measurementFCN;
end

% State measurement
Y_k = zeros(p.M,1);
for i = 1:lsp
    Y_k = Y_k + p.Wm(i)*Ysp_kkm1(:,i);
end

% State measurement covariance
P_YY = zeros(p.M,p.M);
for i = 1:lsp
    P_YY = P_YY + p.Wc(i)*(Ysp_kkm1(:,i) - Y_k)*(Ysp_kkm1(:,i) - Y_k)';
end
P_YY = P_YY + p.R_k;

% State estimation-measurement cross-covariance
P_XY = zeros(p.L,p.M);
for i = 2:lsp
    P_XY = P_XY + (Xsp_kkm1(:,i) - X_kkm1)*(Ysp_kkm1(:,i) - Y_k)';
end
P_XY = (1/(2*(p.lambda+p.L)))*P_XY;

% Kalman gain
K = P_XY*pinv(P_YY);

% State estimation and state estimation error covariance
X_kk = X_kkm1 + K*(Y_meas - Y_k);
P_kk = (P_kkm1 - K*P_YY*K');

% Sigma points of state estimate (k)
Xsp_kk = sigmaPoints( X_kk, P_kk, p.lambda+p.L);

% State-transition model
Xsp_kp1k = zeros(size(Xsp_kkm1));
for i = 1:lsp
    Xsp_kp1k(:,i) = stateTransition( dynamicsFCN, dT_FILTER, Xsp_kk(:,i), p, T_CURRENT);
end

% State estimation
X_kp1k = zeros(p.L,1);
for i = 1:lsp
    X_kp1k = X_kp1k + p.Wm(i)*Xsp_kp1k(:,i);
end

% Covariance estimation
P_kp1k  = zeros(p.L,p.L);
for i = 1:lsp
    P_kp1k = p.Wc(i)*(Xsp_kp1k(:,i) - X_kp1k)*(Xsp_kp1k(:,i) - X_kp1k)';
end
P_kp1k  = P_kp1k + p.Q_k;

INNOVATION = Y_meas - Y_k;


end