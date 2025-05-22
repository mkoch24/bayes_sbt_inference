function [phi_y,log_det_C_y_hat_inv] = computeMarLik(k,y_obs,...
    mu_y_hat,C_y_hat_inv,C_theta_hat_inv,HTCy_invH,log_det_Cy_inv)

% [nlogpy,log_det_C_y_hat_inv] = computeMarginalLikelihood(k,theta,
% theta_mean,C_y_hat_inv,C_theta_hat_inv,HTCy_invH,log_det_Cy_inv)
%
% computes the marginal likelihood (upto a normalization constant) of the
% sample theta and the log-determinant of the covariance matrix
%
% -------- INPUT VARIABLES --------
% k                 = number of layers 
% y_obs             = observation data  (d x 1)
% mu_y_hat          = mean of the normal marginal likelihood distribution  (d x 1) 
% C_y_hat_inv       = covariance of the normal marginal likelihood
% distribution (d x d)
% C_theta_hat_inv   = inverse covariance matrix of the general normal
% conditional posterior distribution over theta (k x k)
% HTCy_invH         = H^TCy^(-1)H (k x k)
% log_det_Cy_inv    = log-determinant of the inverse data covariance matrix
% Cy^-1
% -------- OUTPUT VARIABLES --------
% phi_y             = marginal likelihood of current sample upto a
% normalization constant
% log_det_C_y_hat_inv = log-determinant of the covariance of the normal marginal
% likelihood distribution

% marginal likelihood of current sample upto a normalization constant
phi_y           = 0.5*(y_obs-mu_y_hat)'*C_y_hat_inv*(y_obs- mu_y_hat);

% log-determinant of the covariance matrix using Weinstein–Aronszajn identity
log_det_C_y_hat_inv = log(prod(ones(k,1)-...
  (diag(HTCy_invH)).*(1./diag(C_theta_hat_inv))))+log_det_Cy_inv;
% alternatively can use
% log(det(eye(k)-HTCy_invH/C_theta_hat_inv))-log_det_Cy_inv;

end