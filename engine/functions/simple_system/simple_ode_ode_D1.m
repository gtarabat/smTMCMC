function yres = simple_ode_ode_D1( ode, time, theta, Y0  )
Y0 = Y0(:)';

ODEsolver = ode.solver;

timein = [time(1) time(end)];

try
    sol = ODEsolver( @rhs_adjoint, timein, Y0, ode.Options, ode, theta(:) );
catch
    yres = [];
    return
end

yres = deval(sol,time)';



end