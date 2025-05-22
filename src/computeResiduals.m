function [res_decorr] = computeResiduals(y_obs,z_ind_y,theta,k,d,sample_ind,c,sigma2)

% [res_decorr] = residuals(y_obs,z_ind_y,theta,k,d,obs_ind,sample_ind,
% c,sigma2)

% computes the decorrelated residuals
%
% -------- INPUT VARIABLES --------
% y_obs     = obs data (d x 1)
% z_ind_y   = ids of observation data intervals where interfaces are present ((n_samp x 1) cell)
% theta     = current (initialized to zero) and previous layer properties (k_max x n_samp)
% k         = number of layers of previous sample
% d         = length of obs data vector
% sample_ind= index of random samples chosen between burn_in_samp and n_samp
% c         = proxy to correlation length ((n_samp x 1) cell)
% sigma2    = variance of the GPS ((n_samp x 1) cell)
%
% -------- OUTPUT VARIABLES --------
% res_decorr= decorrelated residual (d x delta_n)

res         = zeros(d,length(sample_ind)); % residual vector
res_decorr  = zeros(d,length(sample_ind)); % decorrelated residual


for i = 1:length(sample_ind)
  c_tmp         = c{sample_ind(i),1};
  sigma2_tmp    = sigma2{sample_ind(i),1};
  z_ind_y_tmp   = z_ind_y{sample_ind(i),1};
  H             = zeros(d,k(sample_ind(i)));

  for j = 2:k(sample_ind(i))+1    
    cov_mat_sym = zeros(z_ind_y_tmp(j)-z_ind_y_tmp(j-1));     
    for m = 2:z_ind_y_tmp(j)-z_ind_y_tmp(j-1)
      cov_mat_sym(m,1:m) = [c_tmp(j-1)^(m-1),cov_mat_sym(m-1,1:m-1)];
    end
    cov_mat = cov_mat_sym + cov_mat_sym'+diag(ones(z_ind_y_tmp(j)-...
      z_ind_y_tmp(j-1),1));
    cov_mat = cov_mat*sigma2_tmp(j-1);
    U = chol(cov_mat);    

    % compute residual
    H(z_ind_y_tmp(j-1)+1:z_ind_y_tmp(j),j-1) = 1;
    res(z_ind_y_tmp(j-1)+1:z_ind_y_tmp(j),i) = y_obs(z_ind_y_tmp(j-1)+1:...
        z_ind_y_tmp(j))-H(z_ind_y_tmp(j-1)+1:z_ind_y_tmp(j),:)*...
        theta(1:k(sample_ind(i)),sample_ind(i));
    res_decorr(z_ind_y_tmp(j-1)+1:z_ind_y_tmp(j),i) = ...
      (res(z_ind_y_tmp(j-1)+1:z_ind_y_tmp(j),i)'/U)';
  end
end


end