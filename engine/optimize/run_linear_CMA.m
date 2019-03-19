function run_linear_CMA( )

  clc; clear

  save_filename = 'CMA_linear';

  data_folder = '../../data/linear/';

  addpath('../functions/tools')
  addpath('../functions/CMA')

  addpath('../functions/linear/')
  
  loglike_func = 'linear_loglike_cma';

  d = load( [data_folder 'data.mat' ] );


  %% Options for CMA optimizer

  start_parallel( 2, false );

  opts.CMA.active = 0;
  opts.PopSize = 10;
  opts.Resume = 0;
  opts.MaxFunEvals = 300000;
  opts.LBounds = [ -5 -5  0 ]'; 
  opts.UBounds = [  5  5 10 ]';
  opts.Noise.on = 0;
  opts.LogModulo = 1;
  opts.LogPlot = 1;
  opts.EvalParallel = 1;
  opts.EvalInitialX = 1;
  opts.TolX = 1e-8;

  opts.SaveFilename      = [ save_filename '.mat'];
  opts.LogFilenamePrefix = [ data_folder 'outcmaes_' save_filename  '_' ];

  xinit = [ 1 1 1 ]';

  X = cmaes_parfor( loglike_func,  xinit,[], opts, d.data );


  fprintf("\n\n Maximum of log-likelihood at:\n" );
  fprintf("  a     = %f \n", X(1) );
  fprintf("  b     = %f \n", X(2) );
  fprintf("  sigma = %f \n\n\n", X(3) );
