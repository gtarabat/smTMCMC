
addpath(genpath('../functions/BASIS')) % Path to BASIS Master

% method = 'BASIS'; cov_check = 'NONE';
method = 'smMALA'; cov_check = 'EIG';

param.cov_check = cov_check;

sys_para = Input_linear(method,param);

BASIS = BASIS_Master(sys_para);


