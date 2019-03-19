%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function: prepare all system parameters for BASIS stored in 'sys_para'
%           and prepare parallel sections (open pool or cluster)
%           created by Stephen Wu, 2015 April 28
%           edited by Daniel Waelchli, Georgios Arampatzis 2016 spring 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Last update: 2016 Jul 25, by Daniel Waelchli
%   2016-07-25: (Daniel)
%       1. update to mTMCMC
%   2015-06-02: (Stephen)
%       1. add N_burn for burn-in period for each MCMC chain
%   2015-06-09: (Stephen)
%       1. add Gen_burn for controling stages of burn-in period
%   2015-09-18: (Stephen)
%       1. add min. ini. value for optimization of p (annealing power)
%       2. add tmp_tol for function tolerance of the optimization
%   2015-09-25: (Stephen)
%       1. change matlabpool to parpool
%   2015-10-24: (Stephen)
%       1. add description on N_s,N_burn,minP,maxP
%   2015-11-03: (Stephen)
%       1. correct variable name for opt_setup in else-statement

function sys_para = Input_SimpleSystem(method, param)

%% Parallelization settings

% Turn on/off parallelization
sys_para.enable_para = true;
% Number of available cpu for parallelization
sys_para.N_cpu = 2;
sys_para.spmd_par = 0;
% Max. length of a chain to control parallel efficiency
% Default: 1 - no bias case (see BASIS paper by Stephen Wu)
sys_para.max_cLen = 1; 

%% Set paths
addpath(genpath('../functions/tools'))
sys_para.model_folder = '../functions/simple_system';
sys_para.data_folder = '../../data/simple_system/';


%% Load ODE synthetic data
load([ sys_para.data_folder 'data_simple_ode.mat' ] );

%% Load the ODE system and solver information
d = load( '../functions/simple_system/simple_system_ode/ode_simple.mat' );
data.ode = d.ode;

%% Set the parallel pool starter
addpath( '../functions/tools' );


%% Set parameters according to method
sys_para.method = method;

switch method
    case 'BASIS'
        loglike_function = 'simple_ode_loglike';
        sys_para.cov_check.method   = 'BASIS';
           
    case 'smMALA'
        loglike_function = 'simple_ode_loglike_smMALA';
        
        switch param.cov_check

                case {'NONE','EIG'}
                    sys_para.cov_check.method   = param.cov_check;

                otherwise
                    error('Unknown cov. check method. Choose NONE, EIG') 
        end
    
    case {'pMALA'}
        loglike_function = 'simple_ode_loglike_pMALA';
        
        switch param.cov_check

                case {'NONE','EIG'}
                    sys_para.cov_check.method   = param.cov_check;

                otherwise
                    error('Unknown cov. check method. Choose NONE, EIG') 
        end

    case 'SN'
        loglike_function = 'simple_ode_loglike_SN';
        
            switch param.cov_check

                case {'NONE','EIG'}
                    sys_para.cov_check.method   = param.cov_check;

                otherwise
                    error('Unknown cov. check method. Choose NONE, EIG') 
            end

    otherwise
       error('Unknown method. Choose BASIS, smMALA, pMALA or SN ') 
       
end

%% Restart settings
% Continue from last generation file
sys_para.restart = false;
% Load file
sys_para.gen_file = '';


%% System settings
% Number of samples at each intermediate stage
sys_para.N_s = 1000;
% Number of samples at the last stage (default: same as N_s)
sys_para.N_s_fin = sys_para.N_s;
% Max. number of stages (default: 1000)
sys_para.max_gen = 10000; 
% Scaling parameter for covariance matrix in MCMC proposal (Gaussian)
% (beta^2 in Ching&Chen 2007, suggested 0.04 in the paper) 
sys_para.beta2 = 0.04;
% scaling parameter of proposal
sys_para.eps = 0.04;
% Boolean for automatic tuning of the scaling param eps during burn-in 
% Note: if set TRUE standard burn-in skipped
sys_para.auto_rescale = false;
% Discard steps durin rescaling, yes/no
sys_para.discard_rescaling = true;
% Max number of repetitions used during auto_rescaling
sys_para.max_rep = 5;
% Target acceptance rate during rescaling
sys_para.tacr = 0.574;
% Target acceptance rate bounds to rescale sys_para.c
sys_para.tacr_bds = [0.4, 0.6];
% Boolean for turning burn-in on/off
sys_para.burnin = true;
% Number of stages for using burn-in in MCMC (no burn-in = 0, always = inf)
% Note: = 1 also gives no burn-in because first stage has no MCMC
sys_para.Gen_burn = 20;
% Length of burn-in period for each MCMC at each stage (a vector)
% Note: length of N_burn = Gen_burn or else last element of N_burn used
%   when current generation > length of N_burn
sys_para.N_burn = [2];
% Number of stages for doing rescaling instead of burn-in in MCMC
% Note: = 1 also gives no rescaling because first stage has no MCMC
sys_para.Gen_rescale = 20;% Length of burn-in period for each MCMC at each stage (a vector)
% Note: length of N_rescale = Gen_burn or else last element of N_rescale used
%   when current generation > length of N_burn
sys_para.N_rescale = 20;
% Boolean for turning system messages on/off
sys_para.TF_msg = true;
% Boolean for saving intermediate results
sys_para.TF_save = true;

% File name for saving output of BASIS_Master
sys_para.save_file = [sys_para.data_folder, method, '_', ...
    param.cov_check, '_', sprintf('%02d',sys_para.N_s)];

%% Model and likelihood settings
% change warning of ode15s on Tolerance to erro for catching
warnId = 'MATLAB:ode15s:IntegrationTolNotMet';
warnstate = warning('error', warnId);

warning('off','MATLAB:nearlySingularMatrix')

% Dimension of theta (main model parameter)
sys_para.N_dim = 3;
% Hard bounds on the parameter space (theta)
% row 1: min, row 2: max, #column = N_dim
sys_para.hard_bds = [ -2 -2  0;
                       2  2 10];
                
% Conf volume to adapt length of eigenvectors (i.e. eigenvalues)
sys_para.conf = 0.68;
% Use extended bounds for adapion of Eigenvalues
sys_para.use_extended_bds = true;
% Factor of boundary extension in each direction relative to hard_bds
% Note: Extended Bounds calculated as
%   upper_bds + width*extended_bds
%   lower_bds - width*extended_bds
sys_para.extend_bds = [1 1 1]; 

if isfield(data,'logParamInd')
    sys_para.hard_bds(:,data.logParamInd) = log(sys_para.hard_bds(:,data.logParamInd));
    sys_para.logParamInd = data.logParamInd;
end

% Matlab file name and input parameters of likelihood function
% Note: likelihood function must be structured as (order specific)
%   output - log(likelihood value), struct of other output
%   input - sample, struct of all necessary input parameters (in .para)
sys_para.lik.name = loglike_function;
sys_para.lik.para = data;

%% Settings for the annealing power

% covariance tolerance for weights normalization  
% (takes value 0.1 - 1, 0.1: slow transition, 1: fast transition-default)
sys_para.cv_tol = 1;

% Settings for limits on power increments
    % 0 = no limit, inf = always use limit
sys_para.N_gen_maxP = 0; % 0 = no limit
    % Note: length of maxP = N_gen_maxP or else last element of maxP used
    %   when current generation > length of maxP
sys_para.maxP = 0.2*ones(sys_para.N_gen_maxP,1);
    % 0 = no limit, inf = always use limit
sys_para.N_gen_minP = 0;
    % Note: length of minP = N_gen_minP or else last element of minP used
    %   when current generation > length of minP
sys_para.minP = 0.0001*ones(sys_para.N_gen_minP,1);
    % Min. initial value for optimization (additional to minP above)
    % Note: it also affects opt. tolerance by default (see below)
sys_para.opt_iniP = 1e-8; % Default = 1e-8

% Settings for optimization for scaling annealing power
% sys_para.opt_setup = optimset('Display','iter','TolX',0.00001,...
%     'TolFun',0.00001,'Algorithm','interior-point');
tmp_tol = min(0.00001,sys_para.opt_iniP);
if sys_para.TF_msg
    sys_para.opt_setup = optimset('Display','iter',...
        'TolX',tmp_tol,'TolFun',tmp_tol,'Algorithm','sqp');
else
    sys_para.opt_setup = optimset('Display','off',...
        'TolX',tmp_tol,'TolFun',tmp_tol,'Algorithm','sqp');
end


%% Prior settings 

% Set up for basic prior options, allow independent pdf for each dimension
% (Ref: http://ch.mathworks.com/help/stats/random.html)
% Each cell contains the name and parameters based on 'random' in Matlab
% Store the parameters as vector in the '.para' cells
for i = 1:sys_para.N_dim
    sys_para.pri.name{i} = 'Uniform';
    sys_para.pri.para{i} = sys_para.hard_bds(:,i);
end

% Set up for customized prior random generator and pdf evaluation
% Require:
%   rnd - function that return matrix of samples (N_s x N_dim)
%       (1 output & 1 input - para_custom_rnd)
%   pdf - function that return pdf values (col. vector) of input samples
%       (1 output & 2 input - samples(N_s x N_dim), para_custom_pdf)
%       (*** make sure input order is same as shown above)
%       (*** make sure the output is in log-scale -> -inf when pdf = 0)
% Note: both functions take single para. variable as input, please use
%   struct for 'para_custom_rnd/pdf' if more than one input parameter
%   variable is needed
% Refer to 'tmcmc_master.m' for use of 'name_custom_rnd' and
%   'tmcmc_worker.m' for use of 'name_custom_pdf'
sys_para.pri.TF_custom = false;
sys_para.pri.name_custom_rnd = '';
sys_para.pri.para_custom_rnd = [];
sys_para.pri.name_custom_pdf = '';
sys_para.pri.para_custom_pdf = [];


%% Settings for proposal prior samples and annealing prior

% Set up for proposal samples in the first stage (sub. prior samples)
% Require:
%   prop_theta - provide N_s x N_dim matrix of initial samples
%   prop_val - provide N_s x 1 col. vector of pdf values of ini. samples
% Note: this is only affects the first stage, and is independent to setting
%   of prior. The prior setting is always required even TF_prop = true.
sys_para.TF_prop_pri = false;
sys_para.prop_theta = [];
sys_para.prop_val = [];

% Set up for annealing prior_o
% Note: suggest to use with TF_prop, i.e., TF_anneal = TF_prop always
sys_para.TF_anneal_pri = sys_para.TF_prop_pri;

%% Settings for rng

% Chose seed for random number generation
% Note: available options are positive integer, none or shuffle
seed = 'shuffle';

switch seed

    case 'shuffle'
        rng('shuffle');
        sys_para.seed = seed;
        
    case 'none'
        % do nothing
        sys_para.seed = seed;
        
    otherwise
        if isnumeric(seed) && seed>0
            rng(seed);
            sys_para.seed = seed;
        else
            error('seed must either be a positive integer, shuffle or none')
        end
end

end




