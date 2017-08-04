function [ X_kp1k, P_kp1k, INNOVATION ] = UKF( dynamicsFCN, measurementFCN,...
    X_kkm1, P_kkm1, Y_meas, p, T_CURRENT, dT_FILTER)
% Unscented Kalman filter algorithm. Uses a proper weightage of the
% estimate and dynamics model 




% Sigma points of state estimate (k-1)
Xsp_kkm1 = sigmaPoints( X_kkm1, P_kkm1, p.lambda+p.L);
lsp = length(Xsp_kkm1);

% Sigma points of state measurement
Ysp_kkm1 = zeros(p.M,lsp);
for j = 1:lsp
    Ysp_kkm1(:,j) = measurementFCN;
end

% State measurement
Y_k = zeros(p.M,1);
for j = 1:lsp
    Y_k = Y_k + p.Wm(j)*Ysp_kkm1(:,j);
end

% State measurement covariance
P_YY = zeros(p.M,p.M);
for j = 1:lsp
    P_YY = P_YY + p.Wc(j)*(Ysp_kkm1(:,j) - Y_k)*(Ysp_kkm1(:,j) - Y_k)';
end
P_YY = P_YY + p.R_k;

% State estimation-measurement cross-covariance
P_XY = zeros(p.L,p.M);
for j = 2:lsp
    P_XY = P_XY + (Xsp_kkm1(:,j) - X_kkm1)*(Ysp_kkm1(:,j) - Y_k)';
end
P_XY = (1/(2*p.lambda))*P_XY;

% Kalman gain
K = P_XY*pinv(P_YY);

% State estimation and state estimation error covariance
X_kk = X_kkm1 + K*(Y_meas - Y_k);
P_kk = (P_kkm1 - K*P_YY*K');

% Sigma points of state estimate (k)
Xsp_kk = sigmaPoints( X_kk, P_kk, p.lambda+p.L);

% State-transition model
Xsp_kp1k = zeros(size(Xsp_kkm1));
for j = 1:lsp
    Xsp_kp1k(:,j) = stateTransition( dynamicsFCN, dT_FILTER, Xsp_kk(:,j), p, T_CURRENT);
end

% State estimation
X_kp1k = zeros(p.L,1);
for j = 1:lsp
    X_kp1k = X_kp1k + p.Wm(j)*Xsp_kp1k(:,j);
end

% Covariance estimation
P_kp1k  = zeros(p.L,p.L);
for j = 1:lsp
    P_kp1k = p.Wc(j)*(Xsp_kp1k(:,j) - X_kp1k)*(Xsp_kp1k(:,j) - X_kp1k)';
end
P_kp1k  = P_kp1k + p.Q_k;

INNOVATION = Y_meas - Y_k;


end