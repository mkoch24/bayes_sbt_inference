function [C_y_inv,log_det_Cy_inv,HTC_y_inv,HTC_y_invH] = ...
  computeCyInvMult(k,d,c,sigma2_inv,id_H,H)  

%[C_y_inv,log_det_Cy_inv,HTC_y_inv,HTC_y_invH] = computeCyInvMult(k,d,c,
% 1./sigma2,id_H,H_old)
%
% computes the components of the tridiagonal inverse data covariance matrix
% (C_y^-1), its log-determinant and its multiples
%
% -------- INPUT VARIABLES --------
% k         = number of layers of previous sample 
% d         = length of observation data vector 
% c         = vector of c, a proxy to corr. length of the GP in each layer (k x 1)
% sigma2_inv= vector of (sigma^2)^-1, the variance of the GP in each layer (k x 1)
% id_H      = ids of observation data intervals where interfaces are present (k+1 x 1)
% H         = matrix mapping the parameter theta to the data (d x k) 
% -------- OUTPUT VARIABLES --------
% C_y_inv           = C_y^-1 matrix (d x d)
% log_det_Cy_inv    = log-determinant of C_y^-1
% HTC_y_inv         = H^T(C_y)^(-1) 
% HTC_y_invH        = H^T(C_y)^(-1)H

% calculate Cy^-1
term1 = sigma2_inv./(1-c.^2);
term2 = -c.*term1;
term3 = (1+c.^2).*term1;

diag1 = zeros(d,1);
diag2 = zeros(d,1);
for i = 1:k
  if id_H(i)+1-id_H(i+1) == 0
    diag1(id_H(i)+1:id_H(i+1),1) = [term1(i,1)];
    diag2(id_H(i)+1:id_H(i+1),1) = 0;
    diag3(id_H(i)+1:id_H(i+1),1) = 0;   
  else
    diag1(id_H(i)+1:id_H(i+1),1) = [term1(i,1);repmat(term3(i,1),...
      id_H(i+1)-id_H(i)-2,1);term1(i,1)];
    diag2(id_H(i)+1:id_H(i+1),1) = [repmat(term2(i,1),id_H(i+1)-id_H(i)-1,1);0];
    diag3(id_H(i)+1:id_H(i+1),1) = [0;repmat(term2(i,1),id_H(i+1)-id_H(i)-1,1)];
  end
end

% sparse format
% C_y_inv = spdiags([diag2 diag1 diag3],-1:1,N,N);

% full format
C_y_inv = diag(diag1,0) + diag(diag2(1:end-1),-1) + diag(diag3(2:end),1);

% log-determinant of Cy^-1
d_diff = diff(id_H);
log_det_Cy_inv = sum(-(d_diff-1).*log(1-c.^2)+d_diff.*log(sigma2_inv));

% calculate H'Cy^(-1) and H'Cy^(-1)H
HTC_y_inv = H'*C_y_inv; 
HTC_y_invH = HTC_y_inv*H;

end