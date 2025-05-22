function [stats] = postProcessing(data, const, state, stats)

% [struct(state), struct(stats)] = postProcessing(struct(data),
% struct(const), struct(state), struct(stats))
%
% post processing to plot all graphs
%
% -------- INPUT VARIABLES --------
% structs: data, const, state, stats. Please see initializeState.m
%
% -------- OUTPUT VARIABLES --------
% updated struct: stats. Please see initializeState.m


d                   = data.d;
y_obs               = data.y_obs;
y_depth             = data.y_depth;

n_samp              = const.n_samp;
burn_in_samp        = const.burn_in_samp;
delta_n             = const.delta_n;
k_min               = const.k_min;
k_max               = const.k_max;
delta_int           = const.delta_int;
r                   = const.r;
z_min               = const.z_min;
z_max               = const.z_max;
theta_min           = const.theta_min;
theta_max           = const.theta_max;
plot_true           = const.flag;
k_true              = const.k_true;
z_true              = const.z_true;
theta_true          = const.theta_true;
sigma2_true         = const.sigma2_true;
l_true              = const.l_true;

k                   = state.k;
z                   = state.z;
theta               = state.theta;
c                   = state.c;
sigma2              = state.sigma2;
z_ind_y             = state.z_ind_y;

post_prob           = stats.post_prob;
phi_y               = stats.phi_y;
n_acc               = stats.n_acc;
n_acc_birth         = stats.n_acc_birth;
n_acc_death         = stats.n_acc_death;
n_acc_perturb       = stats.n_acc_perturb;
n_acc_fixed         = stats.n_acc_fixed;
plot_sbt            = stats.plot_sbt;

%% ========================================================================
% plot log-likelihood

phi_y_mean     = zeros(n_samp,1);
phi_y_mean(1)  = phi_y(1);
for i = 2:n_samp
    phi_y_mean(i) = (phi_y_mean(i-1)*(i-1)+phi_y(i))/i;
end

plotLogMarginalLikelihood(phi_y,n_samp,phi_y_mean) 

%{
fname = sprintf('fig_marginal_log_likelihood.svg');
set(gcf, 'Renderer', 'painters')
print(gcf,'-dsvg',fname);
%}

%% ========================================================================
% compute acceptance probabilities
stats.acc_prob_total    = n_acc/n_samp*100;
stats.acc_prob_birth    = n_acc_birth/n_samp*100;
stats.acc_prob_death    = n_acc_death/n_samp*100;
stats.acc_prob_perturb  = n_acc_perturb/n_samp*100;
stats.acc_prob_fixed    = n_acc_fixed/n_samp*100;

fprintf(['Acceptance probabilities:\n ' ...
    'Total = %0.2f Birth = %0.2f Death = %0.2f Perturb = %0.2f Fixed = %0.2f\n'],...
    stats.acc_prob_total,stats.acc_prob_birth,stats.acc_prob_death,...
    stats.acc_prob_perturb,stats.acc_prob_fixed);

%% ========================================================================
% compute and plot residuals
rng(25,'twister')
sample_ind = randi([burn_in_samp,n_samp],delta_n,1);
[res_decorr] = computeResiduals(y_obs,z_ind_y,theta,k,d,sample_ind,c,sigma2);

res_unc_vec = reshape(res_decorr',delta_n*d,1);
y_depth_rep = repmat(y_depth,1,delta_n);
z_vec       = reshape(y_depth_rep',delta_n*d,1);

x_lab       = 'Residuals';
y_lab       = 'Depth (m)';
add_hist    = 1; % add histogram of average decorrelated residuals
nbins       = 10;
nx          = 50;
ny          = 320;
nor         = 'pdf';
x0          = 25; % x position of fig
y0          = 12;  % y position of fig
width       = 14; % width of fig
height      = 14; % height of fig
lw          = 0.75; % linewidth
col_prior   = [0.65 0.82 0.92]; % colour of the prior
col_hist    = [0 0.1 0.5]; % colour of the histogram

input_param = {add_hist,res_decorr,nbins,col_prior,col_hist,lw*2,x0,y0,...
    width,height};
plotResiduals(nor,res_unc_vec,z_vec,nx,ny,x_lab,y_lab,-3.5,3.5,z_min,...
    z_max,input_param);

%{
fname = sprintf('fig_residuals.svg');
set(gcf, 'Renderer', 'painters')
print(gcf,'-dsvg',fname);
%}

%% ========================================================================
% plot histogram for k and find optimum k
clear pdf;
[Nk,edgek]    = histcounts(k(1,burn_in_samp+1:n_samp));
I_k           = find(Nk == max(Nk));
k_opt         = (edgek(I_k)+edgek(I_k+1))/2; % optimum k
k_prior_y     = 1/(k_max-k_min+1);

x0            = 5;
y0            = 6.9;
width         = 25;
height        = 7;
input_param   = {plot_true,k_true,lw,x0,y0,width,height,col_hist};

plotChainHistk(burn_in_samp,n_samp,k,k_prior_y,k_min,k_max,input_param);

%{
fname = sprintf('fig_k.svg');
set(gcf, 'Renderer', 'painters')
print(gcf,'-dsvg',fname);
%}

%% ========================================================================
% find optimum z based on optimum k and plot conditional posterior pdfs

int_depth_col       = reshape(z(2:k_max,burn_in_samp+1:n_samp),...
    (n_samp-burn_in_samp)*(k_max-1),1);
int_depth_col(int_depth_col==z_min)=[];
int_depth_col(int_depth_col==z_max)=[];
int_depth_col(int_depth_col==0)    =[];

% find position of interfaces (mean and max)
k_index     = k_min:k_max;
z_ind       = cell(15,1);
z_ind_temp  = cell(15,1);

for i = 1:length(k_index)
    z_ind_temp{i,:} = burn_in_samp+find(k(burn_in_samp+1:n_samp)==k_index(i));
    z_ind{i,:} = z(k_min:k_index(i),z_ind_temp{i,:});
end

txt         = 'z';
bw_z        = 0.1;   % bandwidth for kernel density estimate
x0          = 1;
y0          = 2;
width       = k_opt*4;
height      = 8;
col_dens    = col_hist;

input_param = {plot_true,z_true,txt,bw_z,lw,x0,y0,width,height,col_dens};

% find and plot the pdf of interface depths corresponding to the sample
% with max conditional posterior density
[z_max_pos_prob,k_max_pos_prob] = ...
    postProcessMaxDensz(post_prob,z_ind,z_ind_temp,k_index,k_opt,input_param);

%{
fname = sprintf('fig_conditional_post_density.svg');
set(gcf, 'Renderer', 'painters')
print(gcf,'-dsvg',fname);
%}

%% ========================================================================
% plot 2d contours of theta
thin     = 10;  % thin samples of theta at this frequency to plot contours
nx       = 100; %
ny       = d;

theta_N  = zeros(d,n_samp);
for i = 1:n_samp
    for j = 1:k(i)
        theta_N(z_ind_y{i,1}(j)+1:z_ind_y{i,1}(j+1),i) = theta(j,i)*...
            ones(z_ind_y{i,1}(j+1)-z_ind_y{i,1}(j),1);
    end
end

theta_N_vec     = reshape(theta_N(:,burn_in_samp+1:thin:n_samp)',...
    (n_samp-burn_in_samp)/thin*d,1);
y_depth_rep     = repmat(y_depth,1,(n_samp-burn_in_samp)/thin);
z_vec           = reshape(y_depth_rep',(n_samp-burn_in_samp)/thin*d,1);

x_lab           = '\it I_{c}';
y_lab           = 'Depth (m)';
add_ci          = 1;
nor             = 'pdf';
x0              = 32;
y0              = 5;
width           = 12;
height          = 12;

input_param     = {add_ci,d,theta_N,burn_in_samp,n_samp,k_max_pos_prob,...
    z_max_pos_prob,y_obs,lw*3,x_lab,y_lab,z_min,z_max,x0,y0,width,...
    height,theta_min,theta_max,int_depth_col,col_dens,delta_int,r,...
    plot_true,k_true,z_true,theta_true};

plotContourThetaPdfz(nor,theta_N_vec,z_vec,nx,ny,y_depth,input_param);

%{
fname = sprintf('fig_theta_contour_z_pdf.svg');
set(gcf, 'Renderer', 'painters')
print(gcf,'-dsvg',fname);
%}

%% ========================================================================
% determine c_opt and sigma2_opt and plot their posterior marginal pdfs

% returns the mean of the samples of l and sigma^2 corresponding to k_opt
% and the value of l and sigma^2 corresponding to the max. posterior
% density

c_ind_temp      = find(k(burn_in_samp+1:n_samp)==k_opt);
n_k_opt         = length(c_ind_temp);
c_opt           = zeros(k_opt,n_k_opt);
sigma2_opt      = zeros(k_opt,n_k_opt);
for i = 1:length(c_ind_temp)
    c_opt(:,i)      = c{burn_in_samp+c_ind_temp(i),1};
    sigma2_opt(:,i) = sigma2{burn_in_samp+c_ind_temp(i),1};
end
l_opt           = -delta_int./log(c_opt);
sigma2_l        = 1e-4;
sigma2_u        = 1e4;
l_l             = 1e-4;
l_u             = 1e4;

txt             = '\sigma';
bw_sigma2       = 0.1*ones(k_opt,1); % bandwidth for kernel density estimate
log_scale       = 1;
x0              = 1;
y0              = 25;
width           = k_opt*5;
height          = 7;

input_param     = {txt,bw_sigma2,lw*2,x0,y0,width,height,col_dens};
[stats.sigma2_mean,stats.sigma2_max_den] = ...
    plotMarginalsSigmal(sigma2_opt,n_k_opt,sigma2_l,sigma2_u,k_opt,...
    log_scale,plot_true,sigma2_true,input_param);
%{
fname = sprintf('fig_sigma2');
set(gcf, 'Renderer', 'painters')
print(gcf,'-vector', '-dsvg',fname);
%}

x0              = 1;
y0              = 15;
txt             = '\it l';
bw_l            = 0.1*ones(k_opt,1);

input_param     = {txt,bw_l,lw*2,x0,y0,width,height,col_dens};
[stats.l_mean,stats.l_max_den] = ...
    plotMarginalsSigmal(l_opt,n_k_opt,l_l,l_u,k_opt,log_scale,plot_true,...
    l_true,input_param);

%{
fname = sprintf('fig_l');
set(gcf, 'Renderer', 'painters')
print(gcf,'-vector', '-dsvg',fname);
%}
%% ========================================================================
% compute the soil behavior type (sbt) and plot (sbt) figs
if plot_sbt == 1
    z_y                   = [y_depth(1:d-1)  y_depth(2:d)];
    z_ind_y_opt           = zeros(k_opt+1,1);
    z_ind_y_opt(1)        = 0;
    z_ind_y_opt(k_opt+1)  = d;
    theta_mean_avg        = zeros(k_opt,1);
    stats.sbt             = zeros(k_opt,1);

    for i = 1:k_opt-1
        z_ind_y_opt(i+1)    = find(z_max_pos_prob(i,1) > z_y(:,1) &...
            z_max_pos_prob(i,1) <= z_y(:,2));
    end

    x0          = 45;
    y0          = 5;
    width       = 5;
    height      = 12;
    input_param = {burn_in_samp,n_samp,plot_true,k_true,z_true,theta_true,...
        lw*2,x_lab,y_lab,x0,y0,width,height};
    theta_mean  = plotSbt(theta_N,y_obs,y_depth,z_max_pos_prob,k_opt,z_min,...
        z_max,input_param);

    for i = 1:k_opt
        theta_mean_avg(i,1) = mean(theta_mean(z_ind_y_opt(i)+1:...
            z_ind_y_opt(i+1)));
        if theta_mean_avg(i,1) < 1.31
            stats.sbt(i,1) = 2;
        elseif theta_mean_avg(i,1) < 3.6 && theta_mean_avg(i,1) > 2.95
            stats.sbt(i,1) = 3;
        elseif theta_mean_avg(i,1) < 2.95 && theta_mean_avg(i,1) > 2.6
            stats.sbt(i,1) = 4;
        elseif theta_mean_avg(i,1) < 2.6 && theta_mean_avg(i,1) > 2.05
            stats.sbt(i,1) = 5;
        elseif theta_mean_avg(i,1) < 2.05 && theta_mean_avg(i,1) > 1.31
            stats.sbt(i,1) = 6;
        else
            stats.sbt(i,1) = 7;
        end
        fprintf(['Soil Behaviour Type corresponding to layer %i: \n' ...
            '%i \n'],i,stats.sbt(i))
    end


    %{
fname = sprintf('fig_sbt');
set(gcf, 'Renderer', 'painters')
print(gcf,'-vector', '-dsvg',fname);
    %}
end