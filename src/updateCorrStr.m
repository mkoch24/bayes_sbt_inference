function [c_tmp,sigma2_tmp] = updateCorrStr(k,id_z,c,sigma2)

% [c_tmp,sigma2_tmp] = updateCorrStr([k,k_dash],id_z,c,sigma2)
%
% calculates the components of the new covariance matrix depending on the
% type of move.
%
% -------- INPUT VARIABLES --------
% k         = matrix of number of layers containing proposed (k_dash) and 
% previous (k) values (2 x 1)
% id_z      = position id of new interface depth in the updated interface
% depth vector
% c         = vector of c, a proxy to corr. length of the GP in each layer (k(1) x 1)
% sigma2    = vector of sigma^2, the variance of the GP in each layer (k(1) x 1)
% -------- OUTPUT VARIABLES --------
% c         = updated c vector  (k(2) x 1)
% sigma2    = updated sigma^2 vector (k(2) x 1)

% birth move
if k(2) > k(1)

  % modify inverse of existing R (correlation) matrix for the split layer
  % i.e. split it into two parts. New layer has the same properties of the
  % original layer being split

  c_tmp                      = zeros(k(2),1);
  sigma2_tmp                 = zeros(k(2),1);
  c_tmp(id_z:id_z+1)         = c(id_z);
  sigma2_tmp(id_z:id_z+1)    = sigma2(id_z);

  c_tmp(1:id_z-1)            = c(1:id_z-1);
  c_tmp(id_z+2:k(2))         = c(id_z+1:k(1));
  sigma2_tmp(1:id_z-1)       = sigma2(1:id_z-1);
  sigma2_tmp(id_z+2:k(2))    = sigma2(id_z+1:k(1));

% death move
elseif k(2) < k(1)

  % modify inverse of existing R (correlation) matrix for the merged layer
  % i.e. merge two layers into one. 

  c_tmp(1:id_z-1,1)          = c(1:id_z-1);
  c_tmp(id_z:k(2),1)         = c(id_z+1:k(1));
  sigma2_tmp(1:id_z-1,1)     = sigma2(1:id_z-1);
  sigma2_tmp(id_z:k(2),1)    = sigma2(id_z+1:k(1));

% perturbation and no perturbation moves
else

  % do nothing and keep the same correlation matrix

  c_tmp                     = c;
  sigma2_tmp                = sigma2;
end

end