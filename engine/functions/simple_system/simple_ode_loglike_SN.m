function [ res, misc ] = simple_ode_loglike_SN( theta, data )

ode = data.ode;

Nth = length(theta);

Nth_ode = ode.Nth;
Ny  = ode.Ny;

time   = data.x;
ydata  = data.y;

theta_ode = theta(1:Nth_ode);

% Initial data for adjoint system
g0 = ode.g(theta_ode);

dg0 = ode.dg(theta_ode);
d2g0 = ode.d2g(theta_ode);
    
YSH0 = [ g0(:) ; dg0(:) ; d2g0(:) ];

% Solve adjoint system
YSH = simple_ode_ode_D2( ode, time ,  theta_ode, YSH0 );

if isempty(YSH)
    res = -1e8;
    misc.check = 1;
    misc.posdef = false;
    misc.gradient = nan(Nth,1);
    misc.inv_G  = nan(Nth);
    misc.eig.V = nan(Nth);
    misc.eig.D = nan(Nth,1);
    return;  
else
    
    ymodel = sum( YSH(:,1:2) ,2 );
   
    N = length(ymodel);
    
    dif = ydata - ymodel;
    ss = sum( dif.^2 );
      
    term1 = zeros(Nth_ode,1);
    for i = 1:Nth_ode
        k1 = ode.Ny + ode.Ny*(i-1) + 1 ;
        k2 = k1 + ode.Ny - 1;
        Dobs = sum( YSH(:, k1:k2  ) ,2);
        term1(i) = sum( dif.*Dobs );
    end
    
    term2 = zeros(Nth_ode,Nth_ode);
    for i = 1 : Nth_ode
        for j = 1 : i
            k1 = (1+Nth_ode)*Ny + (i>1)*( 0.5*(i-1)*i*Ny +  (j-1)*Ny + 1 ) + (i==1);
            k2 = k1 + Ny - 1;

            sumS = sum( YSH(:, k1:k2)   ,2);
            
            term2(i,j) = sum( dif.*sumS );
            term2(j,i) = term2(i,j);
        end
    end
    
    term3 = zeros(Nth_ode,Nth_ode);
    for i = 1 : Nth_ode
        for j = 1 : i
            k1 = ode.Ny*(i) + 1 ;
            k2 = k1 + ode.Ny - 1;
            k3 = ode.Ny*(j) + 1;
            k4 = k3 + ode.Ny - 1;
            
            sumS1 = sum( YSH(:, k1:k2)   ,2);
            sumS2 = sum( YSH(:, k3:k4)   ,2);
            term3(i,j) = - sum( sumS1.*sumS2 );
            term3(j,i) = term3(i,j);
        end
    end
end

if(Nth>Nth_ode)
    sigma = theta(end);
else
    sigma = data.sigma;
end

sigma2 = sigma.^2;
sigma3 = sigma.^3;
sigma4 = sigma.^4;

% log-likelihood
loglike = -0.5*N*( log(2*pi) + 2*log(sigma) ) - 0.5*ss./sigma2;


%%  Gradient of log-likelihood
D_loglike = zeros(Nth,1);
D_loglike(1:Nth_ode) =   term1 ./ sigma2;

if( Nth > Nth_ode ) % sigma is present
    D_loglike(end)     = -N./sigma + ss./sigma3;
end


%% Hessian of the log-likelihood
D2_loglike = zeros(Nth,Nth);
for i = 1 : Nth_ode
    for j = 1 : i
        tmp = ( term2(i,j) + term3(i,j) ) / sigma2;
        D2_loglike(i,j) = tmp;
        D2_loglike(j,i) = tmp;
    end
end

if( Nth > Nth_ode ) % sigma is present
    i = Nth ;
    tmp = -2*term1/sigma3;
    D2_loglike(i,1:Nth_ode) = tmp;
    D2_loglike(1:Nth_ode,i) = tmp';
    D2_loglike(i,i) = N/sigma2 - 3*ss./sigma4;
end

neg_D2_loglike = - D2_loglike;

inv_D2_loglike = inv(neg_D2_loglike);

if any(isfinite(~inv_D2_loglike(:)))
    res = loglike;
    misc.check = 2;
    misc.posdef = false;
    misc.gradient = D_loglike;
    misc.inv_G  = nan(Nth);
    misc.eig.V = nan(Nth);
    misc.eig.D = nan(Nth,1);
    fprintf('\n inverse of G contains non finite entries \n');
    return
end

[V,D] = eig(inv_D2_loglike);
D = diag(D);
posdef = all( D>0 );

% retrun results
res = loglike;

misc.check = 0;
misc.gradient =  D_loglike;
% inv_G ~ cov proposal distribution
misc.inv_G  = inv_D2_loglike;
% eigenvectors and eigenvalues of inv_G
misc.eig.V = V;
misc.eig.D = D; % eigenvalues in vector
misc.posdef = posdef;

end



