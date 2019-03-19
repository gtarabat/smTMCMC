function [res] = linear_loglike_cma( theta, data )
%
%  return negative of log-likeloihood because CMA minimizes

  [res,~] = linear_loglike( theta, data );
  
  res     = -res;

end
