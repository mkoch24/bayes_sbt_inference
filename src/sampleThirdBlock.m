function [state] = sampleThirdBlock(const, state, coeff)

% [struct(state)] = sampleThirdBlock(struct(const), struct(state)
% struct(coeff))
%
% returns a sample of sigma^2 from an inverse gamma distribution
%
% -------- INPUT VARIABLES --------
% alpha0, beta0 = initial coefficients of the inverse gamma distribution
% id_H  = ids of observation data intervals where interfaces are present (k+1 x 1)
% Y     = (y-H*theta)_i for each layer i (d_i x 1, where d_i is the number of 
% obs. data that lie in the layer i)
% a1    = Y_2^2 + Y^3^2 + ... + Y_(d_i)^2
% a2    = -2*(Y_1Y_2 + Y_2Y_3 + ... + Y_(d_i-1)Y_(d_i))
% a3    = Y_1^2 + Y_2^2 + ... + Y_(d_i)^2
%
% -------- OUTPUT VARIABLES --------
% sigma2= updated vector of sigma^2 (k x 1)  

ip      = state.ip;
alpha0  = const.alpha0;
beta0   = const.beta0;
z_ind_y = state.z_ind_y{ip,1};
c       = state.c{ip,1};
a1      = coeff.a1;
a2      = coeff.a2;
a3      = coeff.a3;

f_c = @(c)((a1.*c.^2+a2.*c+a3)./(1-c.^2)); % YT*R_inv*Y in fxn form
YRY = f_c(c);

% generate new inverse gamma random variable
alpha         = alpha0(1) + 0.5*(diff(z_ind_y));
beta          = beta0(1)  + 0.5*YRY;
sigma2_inv    = gamrnd(alpha,1./beta);
sigma2        = 1./sigma2_inv;

state.sigma2{ip,1} = sigma2;

end