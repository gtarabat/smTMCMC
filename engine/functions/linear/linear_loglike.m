function [ res, misc ] = linear_loglike( theta, data )
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



  Y = theta(1)*data.x + theta(2);

  misc.posdef = 1;

  if isempty(Y) % something went wrong
      
    res = -1e8;
      misc.check = 1;
      return;

  else % all good

      dif   = data.y(:) - Y(:) ; % make sure that both have the same dimension
      sumsq = sum( dif.^2 );

  end

  
  Nd = data.Nd;
  
  sigma = theta(end);

  sigma2 = sigma.^2;

  % log-likelihood
  loglike = -0.5*Nd*( log(2*pi) + 2*log(sigma) ) - 0.5*sumsq./sigma2;

  % retrun results
  res = loglike;
  misc.check = 0;

end



