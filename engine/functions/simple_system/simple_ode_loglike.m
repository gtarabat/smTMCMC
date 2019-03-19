function [ res, misc ] = simple_ode_loglike( theta, data )

ode = data.ode;
Nth_ode = ode.Nth;

Nth = length(theta);

time   = data.x;
ydata  = data.y;

% Initial data for adjoint system
Y0 = ode.g(theta);

% Solve adjoint system
Y = simple_ode_ode( ode, time ,  theta, Y0 );

misc.posdef = 1;

if isempty(Y)
    res = -1e8;
    misc.check = 1;
    return;
else
    
    ymodel = sum( Y(:,1:2) ,2 );
    
    dif = ydata - ymodel;
    ss = sum( dif.^2 );
end

N = length(ymodel);

if(Nth>Nth_ode)
    sigma = theta(end);
else
    sigma = data.sigma;
end

sigma2 = sigma.^2;


% log-likelihood
loglike = -0.5*N*( log(2*pi) + 2*log(sigma) ) - 0.5*ss./sigma2;

% retrun results
res = loglike ;

misc.check = 0;
end



