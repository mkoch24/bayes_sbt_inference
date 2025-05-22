function [state, coeff] = sampleSecondBlock(data, const, state)

%[struct(state), struct(coeff)] = sampleSecondBlock(struct(data),
% struct(const), struct(state))
%
% solves the optimization problem related to block 2 to give c
%
% -------- INPUT VARIABLES --------
% y_obs     = observation data vector (d x 1)
% z_ind_y   = ids of observation data intervals where interfaces are present (k+1 x 1)
% theta     = current sample of theta (k x 1)
% k         = current number of layers
% sigma2    = vector of sigma^2, the variance of the GP in each layer (k x 1)
% c_min      = lower bound on c (0 < c < 1)
% c_max      = upper bound on c (0 < c < 1)
%
% -------- OUTPUT VARIABLES --------
% c_dash    = updated vector of c (k x 1)
% a1,a2,a3  = coefficients of the opt fxn in block 2

ip          = state.ip;
y_obs       = data.y_obs;
z_ind_y     = state.z_ind_y{ip,1};
theta       = state.theta(1:state.k(ip),ip);
k           = state.k(ip);
sigma2      = state.sigma2{ip,1};
c_min       = const.c_min;
c_max       = const.c_max;

c_dash      = zeros(k,1);
a1          = zeros(k,1);
a2          = zeros(k,1);
a3          = zeros(k,1);

for i = 1:k
  % for each layer, compute the opt fxn coefficients
  [a1(i),a2(i),a3(i)] = computeOptFxnCoeff(y_obs(z_ind_y(i)+1:z_ind_y(i+1))- theta(i));

  % define opt fxn
  f_c = @(c)((z_ind_y(i+1)-z_ind_y(i)-1)*log(1-c^2)/2)+...
    0.5*(a1(i)*c^2+a2(i)*c+a3(i))/(sigma2(i)*(1-c^2));
  
  % perform optimization
  % options = optimset('Display','iter'); % uncomment to display iterations
  c_dash(i)  = fminbnd(f_c,c_min,c_max); % modified fminbnd 
end

state.c{ip,1} = c_dash;
coeff.a1      = a1;
coeff.a2      = a2;
coeff.a3      = a3;
end