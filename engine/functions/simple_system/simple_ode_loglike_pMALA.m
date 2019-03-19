function [ res, misc ] = simple_ode_loglike_pMALA( theta, data )
ode = data.ode;

Nth = length(theta);

Nth_ode = ode.Nth;
Ny  = ode.Ny;

time   = data.x;
ydata  = data.y;

theta_ode = theta(1:Nth_ode);

% Initial data for adjoint system
dg0 = ode.dg(theta_ode);
d2g0 = ode.d2g(theta_ode);

YSH0 = [ ode.g(theta_ode) ; dg0(:) ; d2g0(:) ];

% Solve adjoint system
YSH = simple_ode_ode_D2( ode, time ,  theta_ode, YSH0);

if isempty(YSH)
    res = -1e8;
    misc.check = 1;
    misc.posdef = false;
    misc.gradient = nan(Nth,1);
    misc.inv_G  = nan(Nth);
    misc.eig.V = nan(Nth);
    misc.eig.D = nan(Nth,1);
    misc.dGK = nan(Nth,Nth^2);
    return
end
ymodel = sum( YSH(:,1:2) ,2 );

N = length(ymodel);

dif = ydata - ymodel;
ss = sum( dif.^2 );

term1 = zeros(Nth_ode,1);
S = zeros(Nth_ode,N);
for i=1:Nth_ode
    % checked indices
    k1 = Ny + Ny*(i-1) + 1 ;
    k2 = k1 + Ny - 1;
    
    sumS = sum( YSH(:,k1:k2),2);
    term1(i) = sum( dif.*sumS );
    S(i,:) = sumS;
end

S2 = S*S';

P = zeros(Nth_ode,Nth_ode*N);
for i = 1 : Nth_ode
    for j = 1 : i
        % checked indices
        k1 = (1+Nth_ode)*Ny + (i>1)*( 0.5*(i-1)*i*Ny +  (j-1)*Ny + 1 ) + (i==1);
        k2 = k1 + Ny - 1;
        
        sumH = sum( YSH(:, k1:k2)   ,2);
        P(i,((j-1)*N+1):((j-1)*N+N)) = sumH';
        P(j,((i-1)*N+1):((i-1)*N+N)) = P(i,((j-1)*N+1):((j-1)*N+N));
    end
end

if (Nth>Nth_ode)
    sigma = theta(end);
else
    sigma = data.sigma;
end
sigma2 = sigma.^2;
sigma3 = sigma.^3;

% log-likelihood
loglike = -0.5*N*(log(2*pi) + 2*log(sigma)) - 0.5*ss./(sigma2);

%%  Gradient of log-likelihood
D_loglike = zeros(Nth,1);
D_loglike(1:Nth_ode) =   term1 ./ sigma2;

if( Nth > Nth_ode ) % sigma is present
    D_loglike(end)     = -N./sigma + ss./sigma3;
end

%% FIM
G = zeros(Nth);
G(1:Nth_ode,1:Nth_ode) = S2;
if (Nth > Nth_ode)
    G(end,end) = 2*N;
end
inv_G = inv(G)*theta(end)^2;

%% derivative of FIM
dGK = zeros(Nth,Nth^2);
for k = 1:Nth_ode
    for i = 1 : Nth_ode
        for j = 1 : Nth_ode
            dGK(i,j+(k-1)*Nth) = P(i,((k-1)*N+1):((k-1)*N+N))*S(j,:)' + S(i,:)*P(j,((k-1)*N+1):((k-1)*N+N))';
        end
    end
end

if (Nth > Nth_ode)
    dGK(:,Nth*(Nth-1)+1:end) = -2/sigma*G;
end
dGK = dGK/sigma2;

% retrun results
res = loglike;

if any(~isfinite(inv_G(:)))
    misc.check = 2;
    misc.posdef = false;
    misc.gradient = D_loglike;
    misc.inv_G  = nan(Nth);
    misc.eig.V = nan(Nth);
    misc.eig.D = nan(Nth,1);
    misc.dGK = nan(Nth,Nth^2);
    fprintf('\n inverse of G contains non finite entries \n');
    return
end

[V,D] = eig(inv_G);
D = diag(D);
posdef = all( D>0 );

% gradient
misc.gradient =  D_loglike;
% inv_G ~ cov proposal distribution
misc.inv_G  = inv_G;
% eigenvectors and eigenvalues of inv_G
misc.eig.V = V;
misc.eig.D = D; % eigenvalues in vector
misc.posdef = posdef;

% derivative of G
if any(~isfinite(dGK))
    misc.dGK = nan(Nth,Nth^2);
    misc.check = 3;
else
    misc.dGK = dGK;
    %all good
    misc.check = 0;
end

end



