function [ res, misc ] = linear_loglike_smMALA( theta, data )
%
% log-likelihood for the model,
% 
%   y = a*x + b + error,
%
% where error follows normal distribution with zero mean and std sigma
%
% theta(1) = a
% theta(2) = b
% theta(3) = sigma
%
% 

  Nth = length(theta);

  x  = data.x;
  yd = data.y;
  Nd = data.Nd;
  
  % the model
  Y = theta(1)*x + theta(2);
 
  % the derivatives (row vectors)
  dY1 = x(:)';
  dY2 = ones(1,Nd);
  
  
  
  if isempty( Y ) % something went wrong
    
      res = -1e8;
      misc.check = 1;
      misc.posdef = false;
      misc.gradient = nan(Nth,1);
      misc.inv_G  = nan(Nth);
      misc.eig.V = nan(Nth);
      misc.eig.D = nan(Nth,1);
      return;
      
  else % all good
      
      dif   = yd(:)' - Y(:)';
      sumsq = sum( dif.^2 );

      dsum = zeros( Nth-1,  1 ); 
      S    = zeros( Nth-1, Nd );
      
      dsum(1) = sum( dif.*dY1 );
      dsum(2) = sum( dif.*dY2 );
      
      S(1,:) = dY1;
      S(2,:) = dY2;
  
  end

  
  sigma = theta(end);
  sigma2 = sigma.^2;
  sigma3 = sigma.^3;

  % log-likelihood
  loglike = -0.5*Nd*(log(2*pi) + 2*log(sigma)) - 0.5*sumsq./(sigma2);

  %  Gradient of log-likelihood
  D_loglike          = zeros(Nth,1);
  D_loglike(1:end-1) = dsum ./ sigma2;
  D_loglike(end)     = -Nd./sigma + sumsq./sigma3;
  
  
  % Finsher Information matrix
  S2 = S*S';
  FIM = zeros( Nth );
  FIM(1:end-1,1:end-1) = S2;
  FIM(end,end) = 2*Nd;
  

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



  [V,D]  = eig(inv_FIM);
  D      = diag(D);
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



