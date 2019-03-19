function [res] = simple_ode_loglike_cma(theta,data)
[res,~] = simple_ode_loglike(theta,data);
res = -res;
end
