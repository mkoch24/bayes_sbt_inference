%% run_example.m
% example for “Bayesian trans-dimensional soil behaviour type inference”
clc
clear

% 1) add directory to path 
addpath(genpath(fullfile(pwd, 'src')))
addpath(genpath(fullfile(pwd, 'plotting')))
addpath(genpath(fullfile(pwd, 'data')))

rng(25,'twister')

% 2) initialize
[data, const, state, stats] = initializeState( ...
    fullfile(pwd,'data','obs_data_real.mat'), ...data_file
    1e5,...n_samp
    1e4,...burn_in_samp
    1e4,...delta_n % number of posterior samples used to calculate residuals
    2,...  k_min
    20,... k_max 
    3, ... k_init
    [4,8,9,10],...rand_k
    0.05,...delta_int 
    2,...r
    0.05,...z_min
    14.9,...z_max
    [2;5],...z_init
    0.52,...theta_min
    4.12,...theta_max
    2,...theta_init
    [],...mu_theta_prior % initialized to mean of observations by default
    8,...std_dev_theta_prior
    0,...c_min
    1,...c_max
    1,...cor_len_init
    0,...flag indicating if true GP parameters are available (only for synthetic data)
    0); % flag indicating if sbt analysis and plots are required or not

% 3) run the sampler
tic
[state,stats] = runMCMC(data, const, state, stats);
toc

% 4) post processing
stats = postProcessing(data, const, state, stats);

