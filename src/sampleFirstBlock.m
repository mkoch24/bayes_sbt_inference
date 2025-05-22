function [state,stats] = sampleFirstBlock(data, const, state, stats)

%[struct(state), struct(stats)] = sampleFirstBlock(struct(data),
% struct(const), struct(state), struct(stats))
% 
% samples the first block of the algorithm with RJMCMC to yield k,z,theta,
% other variables related to interface position and number of accepted
% samples
%
% -------- INPUT VARIABLES --------
% ip        = current sample
% rand_k    = indicator vector of type of move  (4 x 1) 
% k         = number of layers of previous sample
% rand_num  = random integer determining type of move
% k_min, k_max = minand max number of layers
% z         = interface depth matrix containing current (initialized to
% zero) and previous sample interface depths (k_max+1 x 2)
% delta_int = interval spacing of the grid in which interfaces can be placed
% z_ind_int = interval ids where interfaces are present (k+1 x 1)
% r         = minimum interval spacing between two interfaces
% y_depth   = observation data depths (d x 1)
% n         = total number of intervals
% z_H       = partitioned obs depths
% z_int     = partitioned interval depths
% H_old     = H matrix of the previous sample
% c         = proxy to correlation length
% sigma2    = sigma^2 (k x 1)
% z_ind_y   = ids of observation data intervals where interfaces are present (k+1 x 1)
% theta     = current (initialized to zero) and previous layer properties (k_max x 2)
% post_prob = conditional posterior density of theta
% mu_theta_prior, std_dev_theta_prior = mean and std dev of the prior
% y_obs     = obs data
% thetamin, thetamax, = min and max bounds of theta
% Z_ratio_prior= prior ratio of normalization constants (Z_k/Z_k')
% n_acc, n_acc_birth, n_acc_death, n_acc_perturb, n_acc_fixed = total
% number of accepted samples along with accepted samples in birth, death,
% perturbation and no perturbation move
%
% -------- OUTPUT VARIABLES --------
% k_new     = new number of layers
% z         = interface depth matrix containing current (updated) and 
% previous sample interface depths (k_max+1 x 2) 
% theta     = current (updated) and previous layer properties (k_max x 2)
% sigma2_new= sigma^2 updated to reflect addition/deletion of layer (k_new x 1)
% id_int_new= new interval ids where interfaces are present (k_new+1 x 1) 
% id_H_new  = new ids of observation data intervals where interfaces are present (k_new+1 x 1)
% H_old     = updated H vector 
% post_prob_new= new conditional posterior density of theta
% phi_y     = marginal likelihood of previous sample
% n_acc, n_acc_birth, n_acc_death, n_acc_perturb,n_acc_fixed = total
% number of accepted samples along with accepted samples in birth, death,
% perturbation and no perturbation move

d                   = data.d;
y_obs               = data.y_obs;
y_depth             = data.y_depth; 
z_y                 = data.z_y;

rand_k              = const.rand_k;
k_min               = const.k_min;
k_max               = const.k_max; 
delta_int           = const.delta_int;
r                   = const.r;
n                   = const.n;
z_int               = const.z_int;
theta_min           = const.theta_min;
theta_max           = const.theta_max;
Z_ratio_prior       = const.Z_ratio_prior;
mu_theta_prior      = const.mu_theta_prior;
std_dev_theta_prior = const.std_dev_theta_prior;

ip                  = state.ip;
rand_num            = state.rand_num;
k                   = state.k(ip-1);
z                   = state.z(:,ip-1:ip);
z_ind_int           = state.z_ind_int{ip-1,1};
theta               = state.theta(:,ip-1:ip);
c                   = state.c{ip-1,1}; 
sigma2              = state.sigma2{ip-1,1};
z_ind_y             = state.z_ind_y{ip-1};
H                   = state.H;

post_prob           = stats.post_prob(ip-1);
n_acc               = stats.n_acc;     
n_acc_birth         = stats.n_acc_birth;         
n_acc_death         = stats.n_acc_death;        
n_acc_perturb       = stats.n_acc_perturb;    
n_acc_fixed         = stats.n_acc_fixed; 

flag_z = 0;
[k_dash,q_k,q_k_dash,rand_num] = proposalDensityk(rand_k,k,rand_num,...
    k_min,k_max);

% k_diff  = k_dash-k;
% fprintf('step = %i k_diff = %i \n',ip,k_diff);


%% proposals for z in birth, death, perturbation and no perturbation cases
[z,id_z,flag_z] = findInterfaceDepth(z,k,k_dash,rand_num,delta_int,...
  z_ind_int,r,y_depth,rand_k,flag_z);

%% new observation matrix
H_new              = zeros(d,k_dash);
z_ind_y_tmp        = zeros(k_dash+1,1);
z_ind_y_tmp(1)     = 0;      z_ind_y_tmp(k_dash+1) = d;
z_ind_int_tmp      = zeros(k_dash+1,1);
z_ind_int_tmp(1)   = 0;      z_ind_int_tmp(k_dash+1) = n+1;

if flag_z == 1
  fprintf('There is an error in setting the number of intervals');
  pause
else
  for i = 2:k_dash
    z_ind_y_tmp(i) = find(z(i,2) > z_y(:,1) & z(i,2) <= z_y(:,2));
    H_new(z_ind_y_tmp(i-1)+1:z_ind_y_tmp(i),i-1) = 1;
  end
end
H_new(z_ind_y_tmp(i)+1:d,i) = 1;

for i = 2:k_dash
  z_ind_int_tmp(i) = find(z(i,2) > z_int(:,1) & z(i,2) <= z_int(:,2));
end

% Covariances of the prior distributions of theta
C_theta           = diag(repmat(std_dev_theta_prior^2,k,1));
C_theta_inv       = eye(k)/C_theta;
C_theta_dash      = diag(repmat(std_dev_theta_prior^2,k_dash,1));
C_theta_dash_inv  = eye(k_dash)/C_theta_dash;

%% LHS stationary distribution

% recompute Cy^(-1) and its multiples because c and sigma2 have been
% updated in block 2 and block 3
[C_y_inv,log_det_Cy_inv,HTC_y_inv,HTC_y_invH] = ...
  computeCyInvMult(k,d,c,1./sigma2,z_ind_y,H);

% recompute covariance matrix of LHS posterior
C_theta_hat_inv = HTC_y_invH + C_theta_inv;

if ip == 2
  mu_theta_hat = C_theta_hat_inv\(HTC_y_inv*y_obs +...
    C_theta_inv*repmat(mu_theta_prior,k,1));
  post_prob = (sqrt(prod(diag(C_theta_hat_inv)))/...
    (2*pi())^(k/2))*...
    exp(-0.5*(theta(1:k,1)-mu_theta_hat)'*...
    C_theta_hat_inv*(theta(1:k,1)-mu_theta_hat));
end

% recompute the mean and covariance of the LHS marginal likelihood
mu_y_hat    = H*repmat(mu_theta_prior,k,1);
C_y_hat_inv = C_y_inv - HTC_y_inv'/C_theta_hat_inv*HTC_y_inv;   % faster
% version of C_y_hat_inv = (C_y+H*C_theta*H')^-1;

% recompute marginal likelihood and log-determinant of the covariance
% matrix of the LHS marginal likelihood
[phi_y,log_det_C_y_hat_inv]  = computeMarLik(k,y_obs,...
  mu_y_hat,C_y_hat_inv,C_theta_hat_inv,HTC_y_invH,log_det_Cy_inv);

%% RHS stationary distribution

% first update correlation structure; only components required to define
% the tridiagonal inverse covariance matrix are stored
[c_tmp,sigma2_tmp] = updateCorrStr([k,k_dash],id_z,c,sigma2);

% compute Cy^(-1) and its multiples because a birth/ death/
% perturbation/ no perturbation move has taken place
[C_y_dash_inv,log_det_Cy_dash_inv,HTC_y_dash_inv,HTC_y_dash_invH] = ...
  computeCyInvMult(k_dash,d,c_tmp,1./sigma2_tmp,z_ind_y_tmp,H_new);

% compute mean and covariance of the RHS normal conditional
% posterior distribution over theta
C_theta_hat_dash_inv = HTC_y_dash_invH + C_theta_dash_inv;
mu_theta_hat_dash(1:k_dash,1) = C_theta_hat_dash_inv\(HTC_y_dash_inv*y_obs +...
  C_theta_dash_inv*repmat(mu_theta_prior,k_dash,1));
L = chol(C_theta_hat_dash_inv);

% compute mean and the covariance of the RHS marginal likelihood
mu_y_hat_dash = H_new*repmat(mu_theta_prior,k_dash,1);
C_y_hat_dash_inv = C_y_dash_inv - HTC_y_dash_inv'/C_theta_hat_dash_inv*...
  HTC_y_dash_inv; % faster version of C_y_hat_inv = (C_y+H*C_theta*H')^-1

% compute marginal likelihood and log-determinant of the covariance
% matrix of the RHS marginal likelihood
[phi_y_dash,log_det_C_y_hat_dash_inv]  = computeMarLik(k_dash,...
  y_obs,mu_y_hat_dash,C_y_hat_dash_inv,C_theta_hat_dash_inv,...
  HTC_y_dash_invH,log_det_Cy_dash_inv);

%% Proposal of new sample of theta
% first, need to find the bounds in standard normal space
Xmin   = L*(theta_min(1:k_dash,1)-mu_theta_hat_dash(1:k_dash,1));
Xmax   = L*(theta_max(1:k_dash,1)-mu_theta_hat_dash(1:k_dash,1));

% propose new sample from the truncated standard normal distribution
X_star = trandn(Xmin,Xmax);

% get back the sample in the original space
theta_hat_dash(1:k_dash,1) = L\X_star + mu_theta_hat_dash(1:k_dash,1);

%% compute acceptance probability
% birth move
if rand_num <= rand_k(1)

  mu_theta_hat = C_theta_hat_inv\(HTC_y_inv*y_obs +...
    C_theta_inv*repmat(mu_theta_prior,k,1));

  mu_theta_cdf_tmp = [mu_theta_hat_dash(id_z,1),...
    mu_theta_hat_dash(id_z+1,1),mu_theta_hat(id_z)];
  C_theta_cdf_tmp = [1/sqrt(C_theta_hat_dash_inv(id_z,id_z)),...
    1/sqrt(C_theta_hat_dash_inv(id_z+1,id_z+1)),...
    1/sqrt(C_theta_hat_inv(id_z,id_z))];

  cdf_post(1:3) = normcdf(theta_max(1:3)',mu_theta_cdf_tmp,...
    C_theta_cdf_tmp)-normcdf(theta_min(1:3)',mu_theta_cdf_tmp,...
    C_theta_cdf_tmp);

  Z_hat_ratio           = cdf_post(1)*cdf_post(2)/cdf_post(3);

  p_z_ratio             = priorProbzk(n,k-1,r);
  n_avlbl               = computeAvlblInterval(z_ind_int,r);
  prior_prop_ratio      = p_z_ratio*n_avlbl;

  A = min(1,exp(-phi_y_dash+phi_y+0.5*(log_det_C_y_hat_dash_inv-...
    log_det_C_y_hat_inv))*q_k*Z_hat_ratio/...
    (Z_ratio_prior*q_k_dash)*prior_prop_ratio);

  % death move
elseif rand_num > rand_k(1) && rand_num <= rand_k(2)

  mu_theta_hat = C_theta_hat_inv\(HTC_y_inv*y_obs +...
    C_theta_inv*repmat(mu_theta_prior,k,1));

  mu_theta_cdf_tmp = [mu_theta_hat_dash(id_z-1,1),...
    mu_theta_hat(id_z-1),mu_theta_hat(id_z)];
  C_theta_cdf_tmp = [1/sqrt(C_theta_hat_dash_inv(id_z-1,id_z-1)),...
    1/sqrt(C_theta_hat_inv(id_z-1,id_z-1)),...
    1/sqrt(C_theta_hat_inv(id_z,id_z))];

  cdf_post(1:3) = normcdf(theta_max(1:3)',mu_theta_cdf_tmp,...
    C_theta_cdf_tmp)-normcdf(theta_min(1:3)',mu_theta_cdf_tmp,...
    C_theta_cdf_tmp);

  Z_hat_ratio           = cdf_post(1)/(cdf_post(2)*cdf_post(3));

  p_z_ratio             = priorProbzk(n,k_dash-1,r);
  n_avlbl               = computeAvlblInterval(z_ind_int_tmp,r);
  prior_prop_ratio      = 1/(p_z_ratio*n_avlbl);

  A = min(1,exp(-phi_y_dash+phi_y+0.5*(log_det_C_y_hat_dash_inv-...
    log_det_C_y_hat_inv))*q_k*Z_hat_ratio*Z_ratio_prior/...
    q_k_dash*prior_prop_ratio);

  % perturbation mmove
elseif rand_num == rand_k(3)

  mu_theta_hat = C_theta_hat_inv\(HTC_y_inv*y_obs +...
    C_theta_inv*repmat(mu_theta_prior,k,1));

  mu_theta_cdf_tmp = [mu_theta_hat_dash(id_z-1,1),...
    mu_theta_hat_dash(id_z,1),mu_theta_hat(id_z-1),mu_theta_hat(id_z)];
  C_theta_cdf_tmp = [1/sqrt(C_theta_hat_dash_inv(id_z-1,id_z-1)),...
    1/sqrt(C_theta_hat_dash_inv(id_z,id_z)),...
    1/sqrt(C_theta_hat_inv(id_z-1,id_z-1)),...
    1/sqrt(C_theta_hat_inv(id_z,id_z))];

  if mu_theta_cdf_tmp(1:2) == mu_theta_cdf_tmp(3:4)

    Z_hat_ratio = 1;

  else

    cdf_post(1:4) = normcdf(theta_max(1:4)',mu_theta_cdf_tmp,...
      C_theta_cdf_tmp)-normcdf(theta_min(1:4)',mu_theta_cdf_tmp,...
      C_theta_cdf_tmp);
    Z_hat_ratio = cdf_post(1)*cdf_post(2)/...
      (cdf_post(3)*cdf_post(4));

  end

  A = min(1,exp(-phi_y_dash+phi_y+0.5*(log_det_C_y_hat_dash_inv-...
    log_det_C_y_hat_inv))*q_k*Z_hat_ratio/q_k_dash);

  % no perturbation move
else

  A = min(1,exp(-phi_y_dash+phi_y+0.5*(log_det_C_y_hat_dash_inv-...
    log_det_C_y_hat_inv))*q_k/q_k_dash);

end

%% Metropolis accept-reject criterion
U = rand;

if U <= A
  theta(1:k_dash,2)    = theta_hat_dash(1:k_dash,1);
  n_acc                = n_acc+1;

  post_prob_new = (sqrt(prod(diag(C_theta_hat_dash_inv)))/...
    (2*pi())^(k_dash/2))*...
    exp(-0.5*(theta(1:k_dash,2)-mu_theta_hat_dash(1:k_dash,1))'*...
    C_theta_hat_dash_inv*(theta(1:k_dash,2)-mu_theta_hat_dash(1:k_dash,1)));

  if rand_num <= rand_k(1)
    n_acc_birth       = n_acc_birth + 1;
  elseif rand_num > rand_k(1) && rand_num <= rand_k(2)
    n_acc_death       = n_acc_death + 1;
  elseif rand_num == rand_k(3)
    n_acc_perturb     = n_acc_perturb + 1;
  else
    n_acc_fixed       = n_acc_fixed + 1;
  end
  k_new               = k_dash;
  H                   = H_new;
  z_ind_y_new         = z_ind_y_tmp;
  sigma2_new          = sigma2_tmp;
  z_ind_int_new       = z_ind_int_tmp;
else
  theta(1:k,2)        = theta(1:k,1);
  post_prob_new       = post_prob;
  z(:,2)              = zeros (k_max+1,1);
  z(1:k+1,2)          = z(1:k+1,1);
  k_new               = k;
  z_ind_y_new         = z_ind_y;
  sigma2_new          = sigma2;
  z_ind_int_new       = z_ind_int;
end

state.k(ip)               = k_new;
state.z(:,ip-1:ip)        = z;
state.theta(:,ip-1:ip)    = theta;
state.sigma2{ip,1}        = sigma2_new;
state.z_ind_int{ip,1}     = z_ind_int_new;
state.z_ind_y{ip,1}       = z_ind_y_new;
state.H                   = H;

stats.phi_y(ip)            = phi_y;  
stats.post_prob(ip)        = post_prob_new;
stats.n_acc                = n_acc;
stats.n_acc_birth          = n_acc_birth;
stats.n_acc_death          = n_acc_death;        
stats.n_acc_perturb        = n_acc_perturb;
stats.n_acc_fixed          = n_acc_fixed;

end
