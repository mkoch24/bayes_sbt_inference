function [data, const, state, stats] = initializeState(obsFile, n_samp,...
    burn_in_samp, delta_n, k_min, k_max, k_init, rand_k, delta_int, r,...
    z_min, z_max, z_init, theta_min,theta_max, theta_init, ...
    mu_theta_prior, std_dev_theta_prior,c_min, c_max, cor_len_init,...
    flag_true, flag_sbt)

%
% [data,hyper,state] = initializeState(obsFile, n_samp,burn_in_samp, 
% delta_n, k_min, k_max, k_init, [rand_k], delta_int, r, z_min, z_max,
% [z_init], theta_min,theta_max, theta_init, mu_theta_prior,
% std_dev_theta_prior,c_min, c_max, cor_len_init,flag,flag_sbt)
% 
% load input variables related to 3-block MCMC and build initial structs
%
% -------- INPUT VARIABLES --------
% data fields: y_obs, y_depth, z_H, z_int, d, n
% 
% constant fields: n_samp, burn_in_samp, delta_n, k_min, k_max, rand_k, 
% delta_int, n, z_min, z_max, z_int, r, mu_theta_prior,std_dev_theta_prior,
% theta_min, theta_max, ,Z_ratio_prior, c_min, c_max, alpha0, beta0, flag
% 
% state fields: ip, k, z, theta, c, sigma2, z_ind_int, z_ind_y, H, 
% 
% statistics fields: post_prob, phi_y, n_acc, n_acc_birth, n_acc_death,
% n_acc_perturb, n_acc_fixed
%
% -------- OUTPUT VARIABLES --------
% structs: data, const, state and stats used in 3-block MCMC

% parameters related to data
S                         = load(obsFile,'y_obs','y_depth');
data.y_obs                = S.y_obs; % observation data
data.y_depth              = S.y_depth; % depth of obs data points
data.d                    = numel(data.y_depth); % number of obs data pts
data.z_y                  = [data.y_depth(1:end-1), data.y_depth(2:end)]; % partitioned obs depths

% fixed parameters for the analysis
const.n_samp              = n_samp; % total number of MCMC samples
const.burn_in_samp        = burn_in_samp; % burn-in samples
const.delta_n             = delta_n; % number of posterior samples used to calculate residuals 
const.k_min               = k_min; % minimum number of layers
const.k_max               = k_max; % maximum number of layers
const.rand_k              = rand_k; % indicator vector of type of move  (4 x 1)
const.delta_int           = delta_int; % interval spacing of the grid in which interfaces can be placed
const.n                   = floor(data.y_depth(end)/const.delta_int); % total number of intervals
const.z_min               = z_min; % min interface depth
const.z_max               = z_max; % max interface depth
const.z_int               = [(0:const.n-1)',(1:const.n)']*const.delta_int; % partitioned interval depths
const.r                   = r; % minimum interval spacing between two interfaces
if isempty(mu_theta_prior)
    const.mu_theta_prior  = mean(data.y_obs); % mean of prior distbtn of theta
else
    const.mu_theta_prior  = const.mu_theta_prior;
end
const.std_dev_theta_prior = std_dev_theta_prior; % std dev of the prior distbtn of theta
const.theta_min           = theta_min * ones(k_max,1); % min bound on theta
const.theta_max           = theta_max * ones(k_max,1); % max bound on theta
const.Z_ratio_prior       = normcdf(const.theta_max(1), ...
    const.mu_theta_prior, const.std_dev_theta_prior) ...
    - normcdf(const.theta_min(1), ...
    const.mu_theta_prior, const.std_dev_theta_prior); % prior ratio of normalization constants (Z_k/Z_k')
const.c_min               = c_min; % lower bound on c (0 < c_min < c_max < 1)
const.c_max               = c_max; % upper bound on c (0 < c_max < c_max < 1)
const.alpha0              = 3;    % Inv-Gamma distbtn parameter
const.beta0               = 0.25; % Inv-Gamma distbtn parameter
const.flag                = flag_true; % flag to determine if true GP parameters are available

% define and initialize state arrays related to 3-block MCMC
state.ip                  = 1; % current sample number
state.k                   = zeros(1,n_samp); % number of layers
state.k(1)                = k_init;   % initial number of layers
state.z                   = zeros(k_max+1,n_samp); % layer interface depths 
state.z(2:k_init,1)       = z_init;   % z(2:k_init+1) = [0; z_init; z_max], (k_init+1) x 1 vector           
state.z(k_init+1,1)       = z_max;    % max. interface depth          
state.theta               = zeros(k_max,n_samp); 
state.theta(1:k_init,1)   = theta_init*ones(k_init,1); % initial theta
state.c                   = cell(n_samp,1); % proxy to correlation length
state.c{1}                = exp(-delta_int/cor_len_init).*ones(k_init,1); % initial c
state.sigma2              = cell(n_samp,1); % variance of the GPs
state.sigma2{1}           = std(data.y_obs)^2 * ones(k_init,1); %initial variance
state.z_ind_int           = cell(n_samp,1); % interval ids where interfaces are present (k+1 x 1)
state.z_ind_y             = cell(n_samp,1); % ids of obs data intervals where interfaces are present (k+1 x 1)
state.H                   = zeros(data.d,k_init); % H matrix 

% variables related to statistics in 3-block MCMC
stats.post_prob           = zeros(n_samp,1); % conditional posterior density for each MCMC sample
stats.phi_y               = zeros(1,n_samp); % marginal likelihood for each MCMC sample
stats.n_acc               = 0; % number of accepted samples
stats.n_acc_birth         = 0; % number of accepted samples (birth) 
stats.n_acc_death         = 0; % number of accepted samples (death)
stats.n_acc_perturb       = 0; % number of accepted samples (perturb)
stats.n_acc_fixed         = 0; % number of accepted samples (no pertubation)
stats.plot_sbt            = flag_sbt; % 1 if sbt plots are needed, 0 otherwise

% build initial z_ind_int vector
z_ind_int_tmp      = zeros(k_init+1,1);
z_ind_int_tmp(1)   = 0; z_ind_int_tmp(k_init+1) = const.n+1;
for i=2:3
    z_ind_int_tmp(i) = find(state.z(i,1)>const.z_int(:,1) & ...
        state.z(i,1)<=const.z_int(:,2));
end
state.z_ind_int{1} = z_ind_int_tmp;

% build initial z_ind_y vector
z_ind_y_tmp        = zeros(k_init+1,1);
z_ind_y_tmp(1)     = 0; z_ind_y_tmp(k_init+1) = data.d;
for i=2:3
    z_ind_y_tmp(i) = find(state.z(i,1)>data.z_y(:,1) & ...
        state.z(i,1)<=data.z_y(:,2));
    state.H(z_ind_y_tmp(i-1)+1:z_ind_y_tmp(i),i-1) = 1;
end
state.z_ind_y{1}   = z_ind_y_tmp;

% load true parameters of GPs used to create synthetic data (if available)
if const.flag == 1
    S_true              = load(obsFile,'k_true','z_true','theta_true',...
        'sigma2_true','l_true');
    const.k_true        = S_true.k_true;
    const.z_true        = S_true.z_true;
    const.theta_true    = S_true.theta_true;
    const.sigma2_true   = S_true.sigma2_true;
    const.l_true        = S_true.l_true;
else
    const.k_true        = [];
    const.z_true        = [];
    const.theta_true    = [];
    const.sigma2_true   = [];
    const.l_true        = [];
end