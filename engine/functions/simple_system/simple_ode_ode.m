function yres = simple_ode_ode( ode, time, theta, Y0  )
Y0 = Y0(:)';

ODEsolver = ode.solver;

timein = [time(1) time(end)];

try
    sol = ODEsolver( ode.G, timein, Y0, ode.Options, theta(:) );    
catch
    yres = [];
    return
end

yres = deval(sol,time)';



end