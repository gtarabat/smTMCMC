clc; clear

% number of data
Nd = 40;

x = linspace(1,10,Nd);


% parameters to be inferred 
%
alpha =  2;
beta  = -2;
sigma =  2;


% f(x;alpha,beta) = alpha * x + beta
%
error = normrnd(0,sigma,1,Nd);

y = alpha * x + beta;

yd = y + error;


data.x = x;
data.y = yd;
data.Nd = Nd;

save('data.mat','data')

%% plot data

figure(1); clf
plot(x,y,'LineWidth',3);
hold on
p = plot(x,yd,'o');
p.MarkerSize = 10;
p.MarkerFaceColor = p.Color;

grid on
axis tight

l=legend('exact model','data with noise');
l.Location = 'best';
ax = gca;
ax.FontSize = 15;
