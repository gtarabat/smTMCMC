% clear
% clear
function ode = adjoint_ode_simple()

Nth = 2;
Ny  = 2;

%==========================================================================
% Nth = 7;
% Ny = 4;
%==========================================================================
th = cell(1, Nth);
for i = 1:Nth
    th{i} = sprintf('th%d',i);
end
th = sym(th, 'real');

Y = cell(1, Ny);
for i = 1:Ny
    Y{i} = sprintf('Y%d',i);
end
Y = sym(Y, 'real');


%==========================================================================
%  System definition
%==========================================================================

save_file = 'simple_system_ode/ode_simple.mat';
 
G = sym(zeros(Ny,1));
g = sym(zeros(Ny,1));

%==========================================================================
%  System 1: linear symmetric 2d
%==========================================================================

 G(1) = -th(1)*Y(1);
 G(2) = th(2)*Y(2); 
 
 g(1) = 50;
 g(2) = 50;

Nh = 0.5*Nth*(Nth+1);
%==========================================================================
%   Gradient of G and Initial data w.r.t. Theta
%==========================================================================
[B, dg] = gradient_theta(G,g,th,Ny,Nth);
%==========================================================================
% Gradient of G w.r.t. Y
%==========================================================================
A = gradient_Y(G,Y,Ny,Nth);
%==========================================================================
% Gradient of A w.r.t. Y
%==========================================================================
C = gradient_Y_Y( G, th, Y, Ny, Nth );
%==========================================================================
% Gradient of G  w.r.t. theta
%==========================================================================
[J, d2g] = gradient_theta_theta( G, g, th, Ny, Nth, Nh );
%==========================================================================
% Gradient of A  w.r.t. theta
%==========================================================================
D = gradient_theta_Y( G, th, Y, Ny, Nth );
%==========================================================================
%   'G' and 'g' from symbolic to function handles
%==========================================================================
[Gc,gc] = sym2func(G,g,Ny,Nth);

%==========================================================================
ode.A  = A;
ode.B = B;
ode.C = C;
ode.D = D;
ode.J = J;
ode.dg = dg;
ode.d2g = d2g;
ode.G  = Gc;
ode.g  = gc;
ode.Ny = Ny;
ode.Nth = Nth;
ode.Nh = Nh;
ode.Neq = 1  +  Nth  +  ode.Nh;  % Y + S + H

ode.Options = odeset('RelTol',1e-8);
ode.solver = @ ode15s;

save( save_file ,'ode');

end


%==========================================================================
%==========================================================================
%==========================================================================
%
function A = gradient_Y(G,Y,Ny,Nth)

A = [];I = [];J = [];
for i=1:Ny
    for j=1:Ny
        
        tmp = diff( G(i) , Y(j) );
        
        if( ~isequal(tmp,0) )
            tmp = char(tmp);
            tmp = subvar(tmp,Ny,Nth);
            
            A = [A  tmp ','];
            I  = [ I num2str(i) ','];
            J  = [ J num2str(j) ','];
        end
        
    end
end
A = str2func( [ '@(t,Y,th) sparse( [' I '],[' J '],[' A '],' num2str(Ny) ',' num2str(Ny) ')'] );
end


%==========================================================================
%==========================================================================
%==========================================================================
%
function [B, dg] = gradient_theta(G,g,th,Ny,Nth)

    B = []; IB = []; JB=[];
    dg = []; Idg=[]; Jdg=[];
    for i=1:Ny
        for j=1:Nth

            tmp1 = diff( G(i) , th(j) );
            tmp2 = diff( g(i), th(j) );

            if( ~isequal(tmp1,0) )
                tmp1 = char(tmp1);
                tmp1 = subvar(tmp1,Ny,Nth);
                B = [B  tmp1 ','];
                IB  = [ IB num2str(i) ','];
                JB  = [ JB num2str(j) ','];
            end
            if( ~isequal(tmp2,0) )
                tmp2 = char(tmp2);
                tmp2 = subvar(tmp2,Ny,Nth);
                dg = [dg  tmp2 ','];
                Idg  = [ Idg num2str(i) ','];
                Jdg  = [ Jdg num2str(j) ','];
            end
        end
    end
    
    B = str2func( [ '@(t,Y,th) sparse( [' IB '],[' JB '],[' B '],' num2str(Ny) ',' num2str(Nth) ')'] );
    dg= str2func( [ '@(Y,th) sparse( [' Idg '],[' Jdg '],[' dg '],' num2str(Ny) ',' num2str(Nth) ')'] );
end


%==========================================================================
%==========================================================================
%==========================================================================
%
function C = gradient_Y_Y( G, th, Y, Ny, Nth )
        
    Ck = [];
    I = [];
    J = [];
    for k=1:Ny
        for i=1:Ny
            for j=1:Ny          
                tmp = diff( diff( G(k) , Y(i) ) , Y(j) );
                
                if( ~isequal(tmp,0) )
                    tmp = char( tmp );
                    tmp = subvar(tmp,Ny,Nth);

                    Ck = [Ck  tmp ','];
                    I  = [ I num2str( (k-1)*Ny+i) ','];
                    J  = [ J num2str(j) ','];
                    
                end
            end
        end
        
    end
    C = str2func( [ '@(t,Y,th) sparse( [' I '],[' J '],[' Ck  '],' num2str(Ny*Ny) ',' num2str(Ny) ')'] );
end


%==========================================================================
%==========================================================================
%==========================================================================
%
function [J, d2g] = gradient_theta_theta( G, g, th, Ny, Nth, Nh )
    
    Jk = [];
    d2gk = [];
    IJ=[]; JJ=[];
    Id=[]; Jd=[];
    for  k=1:Nth
        for l=1:k
            
            for i=1:Ny
            
                tmp1 = diff( diff( G(i) , th(k) ) , th(l) );
                tmp2 = diff( diff( g(i) , th(k) ) , th(l) );
                
                ind = (k>1)*( 0.5*(k-1)*k*Ny +  (l-1)*Ny + 1 ) + (k==1);
                ind = ind-1;
                
                if( ~isequal(tmp1,0) )
                    tmp1 = char( tmp1 );
                    tmp1 = subvar(tmp1,Ny,Nth);
                
                    Jk = [ Jk tmp1 ','];
                    IJ = [ IJ num2str(ind+i) ','];
                    JJ = [ JJ num2str(1) ','];
                end
                if( ~isequal(tmp2,0) )
                    tmp2 = char( tmp2 );
                    tmp2 = subvar(tmp2,Ny,Nth);
                    
                    d2gk = [ d2gk tmp2 ','];
                    Id   = [ Id num2str(ind+i) ','];
                    Jd   = [ Jd num2str(1) ','];
                end
             
            end 
        end
    end

    J = str2func( [ '@(t,Y,th) sparse( [' IJ '],[' JJ '],[' Jk  '],' num2str(Nh*Ny) ',' num2str(1) ')'] );
    d2g = str2func( [ '@(t,Y,th) sparse( [' Id '],[' Jd '],[' d2gk  '],' num2str(Nh*Ny) ',' num2str(1) ')'] );

end


%==========================================================================
%==========================================================================
%==========================================================================
%
function D = gradient_theta_Y( G, th, Y, Ny, Nth )
    
    dGdY = cell(Ny,Ny);
    for i=1:Ny
        for j=1:Ny
            dGdY{i,j} = diff( G(i) , Y(j) );
        end
    end
    
    D = cell(Nth,1);
    for k=1:Nth
        Dk = []; I=[]; J=[];
        for i=1:Ny
            for j=1:Ny

                tmp = diff( dGdY{i,j} , th(k) );
                
                if( ~isequal(tmp,0) )
                    tmp = char( tmp );
                    tmp = subvar(tmp,Ny,Nth);
                    
                    Dk = [Dk tmp ','];
                    I  = [ I num2str(i) ','];
                    J  = [ J num2str(j) ','];
                end
                
            end
        end
   
    D{k} = str2func( [ '@(t,Y,th) sparse( [' I '],[' J '],[' Dk '],' num2str(Ny) ',' num2str(Ny) ')'] );
    end
    
    
end


%==========================================================================
%==========================================================================
%==========================================================================
%
function [Gc,gc] = sym2func(G,g,Ny,Nth)
Gc = [];
gc = [];
for i=1:Ny
        
        tmp1 = char( G(i) );
        tmp2 = char( g(i) );
        
        tmp1 = subvar(tmp1,Ny,Nth);
        tmp2 = subvar(tmp2,Ny,Nth);
        
        Gc = [Gc ';' tmp1];
        gc = [gc ';' tmp2];
              
end
Gc = str2func( [ '@(t,Y,th)  [' Gc  ']' ] );
gc = str2func( [ '@(th)  [' gc  ']' ] );

end


%==========================================================================
%==========================================================================
%==========================================================================
%
%   Substitute substring by a string
function tmp = subvar(tmp,Ny,Nth)

for k=1:Ny
    tmp = strrep( tmp, ['Y' num2str(k)] , ['Y(' num2str(k) ')']);
end
for k=1:Nth
    tmp = strrep( tmp, ['th' num2str(k)] , ['th(' num2str(k) ')']);
end

end