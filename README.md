I. SOFTWARE COMPONENTS

You will find the following software packages:

Main algorithm)
- ./engine/functions/BASIS/BASIS_Master.m: control the whole algorithm flow and distribute jobs
- ./engine/functions/BASIS/BASIS_Worker.m: workers that receive individual parallelized job from Master
- ./engine/sample/Input_SimpleSystem.m: prepare system parameters to run BASIS_Master.m (Example)
- ./engine/sample/run_SimpleSystem.m: a dummy file to run the algorithm directly when all setup done (Example)

Supporting materials)
- plot_results_hist.m: plotting the results of BASIS directly
- HANDBOOK.txt: overview of the algorithm structure and recent updates
- README.txt: this file
- AUTHORS: contains all author and contact information
- COPYING: legal terms for using this package

Example) Simple System:

	Y(1) = -th(1)*Y(1);
	Y(2) = th(2)*Y(2); 
 
	y(1) = 50; y(2) = 50;
	f(y|theta,t) = Y(1,t) + Y(2,t)


- run_SimpleSystem.m: main code to run BASIS on the Simple System
- Input_SimpleSystem.m: setup system parameters for the example


II. PREREQUISITES

- MATLAB with optimization and parallel package (this package uses parpool)


III. WRITE YOUR OWN SYSTEM
- copy /engine/functions/simple_system folder
- create folder in /data/
- name your system and rename simple_ode_.. functions
- define your own system in adjoin_ode_simple
- define function of solution of simple system f(y|theta,t) in 
	your_ode_loglike_cma.m
	your_ode_loglike_pMALA.m
	your_ode_loglike_smMALA.m
	your_ode_loglike_SN.m
	your_ode_loglike.m
- run make data, in make data you can specify model parameters, theta and error, respectively sigma
- move new data file from /engine/functions/your_system/your_system_ode/ to /data/your_system/
- copy & rename /engine/sample/Input_SimpleSystem.m
- copy & rename /engine/sample/run_SimpleSystem.m
- adjust Input_YourSystem.m, run_YourSystem.m, especially check function handles and paths, such that it refers to your_ode_..


IV. NOTES

Read HANDBOOK.txt for overview of the algorithm usage
Look into Input_SimpleSystem.m for the important parameters and their purpose

Ask for help, report any problems!: garampat@ethz.ch, wadaniel@ethz.ch

