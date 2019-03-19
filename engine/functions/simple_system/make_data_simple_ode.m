clear; clc

disp('run functions/simple_system/adjoint_ode_simple.m to update ODE system');
addpath(genpath('simple_system_ode'));

ODEsolver = @ ode15s;

theta0 = [ 0.4; 0.2];

%==========================================================================
% Define the ODE system + Adjoint equations
%==========================================================================
ode = adjoint_ode_simple();
%load('ode_simple.mat');
Ny  = ode.Ny;
Nth = ode.Nth;

T = 5; 
tin = [0 T]; 
tt = 0:0.01:T;

Options = odeset('RelTol',1e-8);

%% ========================================================================
sol = ODEsolver( ode.G, tin, ode.g(theta0) , Options, theta0 );

y = deval(sol,tt)';

%% ========================================================================
td = 0:0.5:T;
sigma = 2;
error = sigma*randn(length(td),1);

yd = sum( deval(sol,td) ,1); 
yd = yd' + error;

data.x = td;
data.y = yd;
data.theta = theta0;
data.sigma = sigma;

save('simple_system_ode/data_simple_ode','data')


%% ========================================================================
subplot(2,1,1)
plot(tt, sum(y,2));
hold on;
grid on;
plot(td,yd,'*')

subplot(2,1,2)
plot(tt, y);



