function [ SP ] = sigmaPoints( X, P, c)
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

%% Cholesky factorization (matrix square-root)
A = sqrt(c)*chol(P)';

%% Sigma point calculation
Y = X(:,ones(1,numel(X)));
SP = [X, Y+A, Y-A];

end