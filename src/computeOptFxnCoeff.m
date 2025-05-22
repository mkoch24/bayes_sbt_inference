function [a1,a2,a3] = computeOptFxnCoeff(Y)

% [a1(i),a2(i),a3(i)] = compute_opt_fxn_coeff(y_obs(id_H(i)+1:id_H(i+1))- theta(i))
%
% returns the opt fxn coeffs
%
% -------- INPUT VARIABLES --------
% Y  = (y-H*theta)_i for each layer i (d_i x 1, where d_i is the number of 
% obs. data that lie in the layer i) 
% -------- OUTPUT VARIABLES --------
% a1 = Y_2^2 + Y^3^2 + ... + Y_(d_i)^2
% a2 = -2*(Y_1Y_2 + Y_2Y_3 + ... + Y_(d_i-1)Y_(d_i))
% a3 = Y_1^2 + Y_2^2 + ... + Y_(d_i)^2 

a1  = sum(Y(2:end-1).^2);

a21 = Y(1:end-1);
a22 = Y(2:end);
a2  = -2*sum(a21.*a22);

a3 = sum(Y(1:end).^2);

end