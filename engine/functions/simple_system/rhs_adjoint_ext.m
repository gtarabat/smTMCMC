function  res = rhs_adjoint_ext( ~, YSH, ode, theta)


Ny = ode.Ny;

YSH = reshape( YSH , Ny, ode.Neq);


% ode RHS
Y  = YSH(:,1);
YY = ode.G( 0,Y,theta );

% sensitivities RHS
S  = YSH(:,2:ode.Nth+1);
aa = ode.A( 0, Y, theta);
bb = ode.B( 0, Y, theta);
     
SS = aa*S + bb;


% second derivatives RHS
H  = YSH(:,ode.Nth+2:end);
HH = zeros( length(H), 1);

A = ode.A( 0, Y, theta);
C = ode.C(0,Y,theta);

D = cell(ode.Nth,1);
for k=1:ode.Nth
    D{k} = ode.D{k}(0,Y,theta);
end

J = ode.J( 0, Y, theta);

Csl = cell(ode.Nth,1);
for l=1:ode.Nth
    sl = S(:,l);
    Csl{l} = reshape(C*sl,Ny,Ny);
end
    
for k=1:ode.Nth
    for l=1:k

        ind_H = ( 0.5*(k-1)*k + l );
        ind_HH_1 = (k>1)*( 0.5*(k-1)*k*Ny +  (l-1)*Ny + 1 ) + (k==1);
        ind_HH_2 = ind_HH_1 + Ny - 1;
        
        h  = H(:,ind_H);
        sk = S(:,k);
        sl = S(:,l);

        sk = sk(:); sl = sl(:);

        tmp = (sk'*Csl{l})';

        tmp = tmp+D{l}*sk;
        
        tmp = tmp + D{k}*sl;

        tmp = tmp + A*h(:);

        tmp = tmp+J(ind_HH_1:ind_HH_2);

        HH(ind_HH_1:ind_HH_2) = tmp;
        
    end
end

res = [ YY(:) ; SS(:) ; HH(:)];