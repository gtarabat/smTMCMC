%% Script to plot final results

%% Clear all
clear all
clf
close all

%% Data
% Assuming dimension > 1 for theta
%==========================================================================
%==========================================================================

file_name = '../../data/simple_system/BASIS_NONE_50.mat';
load(file_name);
plotmatrix_hist(out_master.theta);

load(file_name);
plotmatrix_hist(out_master.theta);


[bestlik,idx] = max(out_master.lik);
besttheta = out_master.theta(idx,:);

fprintf('best likelihhood %f\n',bestlik)
fprintf('found at: ')
fprintf('%7.3f' , besttheta)
fprintf('\n')

return
%% Plot generations

pause
for k = 2:out_master.gen
    load([file_name '_gen_'  sprintf('%02d',k)  '.mat'])
    plotmatrix_hist(out_master.theta);

pause
end

