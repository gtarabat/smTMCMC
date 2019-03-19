function [ res, misc ] = simple_ode_loglike_smMALA( theta, data )
ode = data.ode;

Nth = length(theta);

Nth_ode = ode.Nth;
Ny  = ode.Ny;

theta_ode = theta(1:Nth_ode);

time   = data.x;
ydata  = data.y;

% Initial data for adjoint system
y0  = ode.g(theta_ode);
dy0 = ode.dg(theta_ode);
YS0 = [ y0 ; dy0(:) ];

% Solve adjoint system
YS = simple_ode_ode_D1( ode, time , theta_ode, YS0);

if isempty(YS)
    res = -1e8;
    misc.check = 1;
    misc.posdef = false;
    misc.gradient = nan(Nth,1);
    misc.inv_G  = nan(Nth);
    misc.eig.V = nan(Nth);
    misc.eig.D = nan(Nth,1);
    return;
else
    ymodel = sum( YS(:,1:2) ,2 );

    N = length(ymodel);
    
    dif = ydata - ymodel;
    ss = sum( dif.^2 );
   
    term1 = zeros(ode.Nth,1);
    S = zeros(Nth_ode,N);
    for i=1:Nth_ode
        k1 = Ny + Ny*(i-1) + 1 ;
        k2 = k1 + ode.Ny - 1;
        sumS = sum( YS(:, k1:k2  ) ,2);
        term1(i) = sum( dif.*sumS );
        S(i,:) = sumS;
    end
    
end

S2 = S*S';

if(Nth>Nth_ode)
    sigma = theta(end);
else
    sigma = data.sigma;
end

sigma2 = sigma.^2;
sigma3 = sigma.^3;

% log-likelihood
loglike = -0.5*N*(log(2*pi) + 2*log(sigma)) - 0.5*ss./(sigma2);

%  Gradient of log-likelihood
D_loglike = zeros(Nth,1);

D_loglike(1:Nth_ode) =   term1 ./ sigma2;

if(Nth>Nth_ode) % sigma is present
    D_loglike(end)   = -N./sigma + ss./sigma3;
end
    
FIM = zeros(Nth);
FIM(1:Nth_ode,1:Nth_ode) = S2;
if(Nth>Nth_ode) % sigma is present
    FIM(end,end) = 2*N;
end

inv_FIM = inv(FIM)*theta(end)^2;

% retrun results
res = loglike;

if any(~isfinite(inv_FIM(:)))
    misc.check = 2;
    misc.posdef = false;
    misc.gradient = D_loglike;
    misc.inv_G  = nan(Nth);
    misc.eig.V = nan(Nth);
    misc.eig.D = nan(Nth,1);
    fprintf('\n inverse of G contains non finite entries \n');
    return
end
    
[V,D] = eig(inv_FIM);
D = diag(D);
posdef = all( D>0 );

misc.check = 0;
misc.gradient =  D_loglike;
% inv_G ~ cov proposal distribution
misc.inv_G  = inv_FIM;
% eigenvectors and eigenvalues of inv_G
misc.eig.V = V;
misc.eig.D = D; % eigenvalues in vector
misc.posdef = posdef;

end



