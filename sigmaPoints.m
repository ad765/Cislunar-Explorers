function [ SP ] = sigmaPoints( X, P, lambda)
% Calculates sigma points of state vector.
%
% Inputs:
%           X       - State vector
%           P       - State covariance matrix
%           lambda  - Scaling factor from relationship
%                       lambda = alpha^2*(L+kappa)
%
% Outputs:
%           SP      - Sigma points (ns  by  2*L+1 array)
%
% Anshuman Das, Cornell University
% Wednesday, August 2, 2018

%% Initialization
L = length(X);
SP = zeros(L,2*L);

%% Cholesky factorization (matrix square-root)
A = chol(P,'upper');

%% Sigma point calculation
for i = 1:L
    SP(:,i)   = X + sqrt(lambda)*A(:,i);
    SP(:,i+L) = X - sqrt(lambda)*A(:,i);
end
SP = [X, SP];

end