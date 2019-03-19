function run_SimpleSystem_CMA( )
% Note: Script crashes during first run, restart and it works :-)
% Note: Problems with parallel startup

clc; clear

save_filename = 'CMA_simple_system';

data_folder = '../../data/simple_system/';
ode_folder = '../functions/simple_system/simple_system_ode/';

addpath('../functions/tools')
addpath('../functions/CMA')
addpath('../functions/simple_system/')

loglike_func = 'simple_ode_loglike_cma';

load( [ode_folder 'ode_simple.mat' ] );
load( [data_folder 'data_simple_ode.mat' ] );

data.ode = ode;

%% Options
% change warning of ode15s on Tolerance to error for catching
warnId = 'MATLAB:ode15s:IntegrationTolNotMet';
warnstate = warning('error', warnId);

start_parallel( 48, false );

opts.CMA.active = 0;
opts.PopSize = 48;
opts.Resume = 0;
opts.MaxFunEvals = 300000;
opts.LBounds = [-2 -2 0]'; 
opts.UBounds = [2 2 10]';
opts.Noise.on = 0;
opts.LogModulo = 1;
opts.LogPlot = 1;
opts.EvalParallel = 1;
opts.EvalInitialX = 1;
opts.TolX = 1e-8;

opts.SaveFilename      = [ save_filename '.mat'];
opts.LogFilenamePrefix = [ data_folder 'outcmaes_' save_filename  '_' ];

xinit = [0.36 0.19 8.2]';

X = cmaes_parfor( loglike_func,  xinit,[], opts,data);


