function [h,ax,BigAx,hhist,pax] = plotmatrix_hist(theta, varargin)


sym = '.'; % Default scatter plot symbol.

rows = size(theta,2); cols = rows;
x = theta; y = theta;
dohist = 1;

% Don't plot anything if either x or y is empty
hhist = gobjects(0);
pax   = gobjects(0);
if isempty(rows) || isempty(cols),
    if nargout>0, h = gobjects(0); ax = gobjects(0); BigAx = gobjects(0); end
    return
end

if ~ismatrix(x) || ~ismatrix(y)
    error(message('MATLAB:plotmatrix:InvalidXYMatrices'))
end

if size(x,1)~=size(y,1) || size(x,3)~=size(y,3),
    error(message('MATLAB:plotmatrix:XYSizeMismatch'));
end



% Create/find BigAx and make it invisible
cax = axes;
BigAx = newplot(cax);
fig = ancestor(BigAx,'figure');
hold_state = ishold(BigAx);
set(BigAx,'Visible','off','color','none')

if any(sym=='.'),
    units = get(BigAx,'units');
    set(BigAx,'units','pixels');
    pos = get(BigAx,'Position');
    set(BigAx,'units',units);
    markersize = max(1,min(15,round(15*min(pos(3:4))/max(1,size(x,1))/max(rows,cols))));
else
    markersize = get(0,'DefaultLineMarkerSize');
end



% Create and plot into axes
ax = gobjects(rows,cols);
pos = get(BigAx,'Position');
width = pos(3)/cols;
height = pos(4)/rows;
space = .05; % 2 percent space between axes
pos(1:2) = pos(1:2) + space*[width height];
m = size(y,1);
k = size(y,3);
xlim = zeros([rows cols 2]);
ylim = zeros([rows cols 2]);
BigAxHV = get(BigAx,'HandleVisibility');
BigAxParent = get(BigAx,'Parent');
paxes = findobj(fig,'Type','axes','tag','PlotMatrixScatterAx');
for i=rows:-1:1,
    for j=cols:-1:1,
        axPos = [pos(1)+(j-1)*width pos(2)+(rows-i)*height width*(1-space) height*(1-space)];
        findax = findaxpos(paxes, axPos);
        if isempty(findax),
            ax(i,j) = axes('Position',axPos,'HandleVisibility',BigAxHV,'parent',BigAxParent);
            set(ax(i,j),'visible','on');
        else
            ax(i,j) = findax(1);
        end

        
        tmp1 = reshape(x(:,j),[m k]);
        tmp2 = reshape(y(:,i),[m k]);
        
%         hh(i,j,:) = plot( tmp1, tmp2 ,sym,'parent',ax(i,j))';
        
        if(i<j)
%             hh(i,j,:) = histogram2(tmp1,tmp2, 'DisplayStyle','tile','ShowEmptyBins','on');
%             hh(i,j,:).Normalization = 'pdf';
%             hh(i,j,:).NumBins = [50 50];
%             hh(i,j,:).EdgeColor ='none';
%             colormap parula
            
            if(nargin > 1)
                hh(i,j) = scatter(tmp1,tmp2,[], varargin{1});
            else
                hh(i,j) = plot(tmp1,tmp2,sym);
            end

            
            ax(i,j).XLim = ax(j,i).YLim;
            ax(i,j).YLim = ax(j,i).XLim;
            
        else
            [Z,Xe,Ye]= histcounts2(tmp2,tmp1,20,'Normalization','pdf');
            
            X = Xe(1:end-1) + diff(Xe)/2; 
            Y = Ye(1:end-1) + diff(Ye)/2;
            
            [X,Y] = meshgrid( X, Y );  
            
            Z = smoothn(Z,1);
            

            nn = 100;
            xi = linspace(X(1), X(end), nn);
            yi = linspace(Y(1), Y(end), nn);
            [xi, yi] = meshgrid(xi, yi);
            zi = interp2(X, Y, Z', xi, yi, 'spline');
            
%             [~,hh(i,j,:)] = contour(yi, xi, zi, 30 ,'parent',ax(i,j));
            
            hh(i,j,:) = surf(yi, xi, zi, 'parent',ax(i,j));
            hh(i,j).FaceColor = 'interp';
            hh(i,j).EdgeColor = 'none';
            view(ax(i,j), [0,90]) 


            ax(i,j).XLim = [Y(1) , Y(end)];
            ax(i,j).YLim = [X(1) , X(end)];
        end
        
        disp([i,j])

        set(ax(i,j),'xgrid','off','ygrid','off')

        xlim(i,j,:) = get(ax(i,j),'xlim');
        ylim(i,j,:) = get(ax(i,j),'ylim');        
        
    end
    

end


xlimmin = min(xlim(:,:,1),[],1);
xlimmax = max(xlim(:,:,2),[],1);
ylimmin = min(ylim(:,:,1),[],2);
ylimmax = max(ylim(:,:,2),[],2);



% Try to be smart about axes limits and labels.  Set all the limits of a
% row or column to be the same and inset the tick marks by 10 percent.
inset = .05;
for i=1:rows,
    set(ax(i,1),'ylim',[ylimmin(i,1) ylimmax(i,1)])
    dy = diff(get(ax(i,1),'ylim'))*inset;
    set(ax(i,:),'ylim',[ylimmin(i,1)-dy ylimmax(i,1)+dy])
end
dx = zeros(1,cols);
for j=1:cols,
    set(ax(1,j),'xlim',[xlimmin(1,j) xlimmax(1,j)])
    dx(j) = diff(get(ax(1,j),'xlim'))*inset;
    set(ax(:,j),'xlim',[xlimmin(1,j)-dx(j) xlimmax(1,j)+dx(j)])
end

set(ax(1:rows-1,:),'xticklabel','')
set(ax(:,2:cols),'yticklabel','')
set(BigAx,'XTick',get(ax(rows,1),'xtick'),'YTick',get(ax(rows,1),'ytick'), ...
    'userdata',ax,'tag','PlotMatrixBigAx')
set(ax,'tag','PlotMatrixScatterAx');

if dohist, % Put a histogram on the diagonal for plotmatrix(y) case
    paxes = findobj(fig,'Type','axes','tag','PlotMatrixHistAx');
    pax = gobjects(1, rows);
    for i=rows:-1:1,
        axPos = get(ax(i,i),'Position');
        findax = findaxpos(paxes, axPos);
        if isempty(findax),
            histax = axes('Position',axPos,'HandleVisibility',BigAxHV,'parent',BigAxParent);
            set(histax,'visible','on');
        else
            histax = findax(1);
        end
        
        hhist(i) = histogram(histax,y(:,i,:),50);
        hhist(i).Normalization = 'pdf';
        hhist(i).FaceColor = [51,255,51]/255;

        
%         hold on
%         xx = hhist(i).BinEdges + hhist(i).BinWidth(1)/2; xx=xx(1:end-1);
%         plot(histax,xx',smoothn(hhist(i).Values,0.1))
        

        
        set(histax,'xtick',[],'ytick',[],'xgrid','off','ygrid','off');
        set(histax,'xlim',[xlimmin(1,i)-dx(i) xlimmax(1,i)+dx(i)])
        set(histax,'tag','PlotMatrixHistAx');
        pax(i) = histax;  % ax handles for histograms
    end
end

% Make BigAx the CurrentAxes
set(fig,'CurrentAx',BigAx)
if ~hold_state,
    set(fig,'NextPlot','replacechildren')
end

% Also set Title and X/YLabel visibility to on and strings to empty
set([get(BigAx,'Title'); get(BigAx,'XLabel'); get(BigAx,'YLabel')], ...
    'String','','Visible','on')

if nargout~=0,
    h = hh;
end



function findax = findaxpos(ax, axpos)
tol = eps;
findax = [];
for i = 1:length(ax)
    axipos = get(ax(i),'Position');
    diffpos = axipos - axpos;
    if (max(max(abs(diffpos))) < tol)
        findax = ax(i);
        break;
    end
end
