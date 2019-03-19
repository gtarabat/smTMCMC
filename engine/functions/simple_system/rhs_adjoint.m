function  res = rhs_adjoint( ~, YS, ode, theta)


YS = reshape( YS , ode.Ny, ode.Nth+1);


% ode RHS
Y  = YS(:,1);
YY = ode.G( 0,Y,theta );

% sensitivities RHS
S  = YS(:,2:end);
aa = ode.A ( 0, Y, theta);
bb = ode.B( 0, Y, theta);
     
SS = aa*S + bb;

res = [ YY(:) ; SS(:) ];