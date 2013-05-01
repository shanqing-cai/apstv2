function t_int_sum(chg_IUInt,chg_IYInt,chg_IU2Int,chg_IY2Int,chg_IU3Int,chg_IY3Int,colors,varargin)
%% Config
lw=2;
fontSize=12;
alpha=0.025;

windowWidth=500;
windowHeight=260;

if ~isempty(fsic(varargin,'fontSize'))
    fontSize=varargin{fsic(varargin,'fontSize')+1};
end
if ~isempty(fsic(varargin,'windowWidth'))
    windowWidth=varargin{fsic(varargin,'windowWidth')+1};
end
if ~isempty(fsic(varargin,'windowHeight'))
    windowHeight=varargin{fsic(varargin,'windowHeight')+1};
end

%%
tab_int_chg.accel=[chg_IUInt(:,1),chg_IYInt(:,1),chg_IU2Int(:,1),chg_IY2Int(:,1),chg_IU3Int(:,1),chg_IY3Int(:,1)];
mean_int_chg.accel=1e3*mean(tab_int_chg.accel);
ste_int_chg.accel=1e3*ste(tab_int_chg.accel);
             
tab_int_chg.decel=[chg_IUInt(:,2),chg_IYInt(:,2),chg_IU2Int(:,2),chg_IY2Int(:,2),chg_IU3Int(:,2),chg_IY3Int(:,2)];
mean_int_chg.decel=1e3*mean(tab_int_chg.decel);
ste_int_chg.decel=1e3*ste(tab_int_chg.decel);

pvals_pairedT=nan(1,numel(mean_int_chg.accel));
pvals_accel=nan(1,numel(mean_int_chg.accel));
pvals_decel=nan(1,numel(mean_int_chg.decel));
for i1=1:numel(mean_int_chg.accel)
    [h,t_p]=ttest(tab_int_chg.accel(:,i1),tab_int_chg.decel(:,i1));    
    pvals_pairedT(i1)=t_p;
    [h,t_p]=ttest(tab_int_chg.accel(:,i1));
    pvals_accel(i1)=t_p;
    [h,t_p]=ttest(tab_int_chg.decel(:,i1));
    pvals_decel(i1)=t_p;
end
            
if isempty(fsic(varargin, 'noFigure'))
    figure('Position',[200,200,windowWidth,windowHeight]);
end
set(gca,'FontSize',fontSize);
errorbar(1:numel(mean_int_chg.accel),mean_int_chg.accel,ste_int_chg.accel,'color',colors.accel,'LineWidth',lw); hold on;
errorbar(1:numel(mean_int_chg.decel),mean_int_chg.decel,ste_int_chg.decel,'color',colors.decel,'LineWidth',lw);
plot(1:numel(mean_int_chg.accel),mean_int_chg.accel,'o','LineWidth',lw,'color',colors.accel);
plot(1:numel(mean_int_chg.decel),mean_int_chg.decel,'o','LineWidth',lw,'color',colors.decel);

if ~isempty(fsic(varargin, 'YLim'))
    y_lim = varargin{fsic(varargin, 'YLim') + 1};
    set(gca, 'YLim', y_lim);
end

xs=get(gca,'XLim'); ys=get(gca,'YLim');
for i1=1:numel(pvals_pairedT)
    if pvals_pairedT(i1)<alpha
        plot(i1,ys(2)-0.05*range(ys),'*','color','k','MarkerSize',11);
    end
    if pvals_accel(i1)<alpha
        plot(i1, mean_int_chg.accel(i1), 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', colors.accel);
    else
        plot(i1, mean_int_chg.accel(i1), 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w');
    end
    if pvals_decel(i1)<alpha
        plot(i1, mean_int_chg.decel(i1), 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', colors.decel);
    else
        plot(i1, mean_int_chg.decel(i1), 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w');
    end
end
plot(xs(1)+0.05*range(xs),ys(1)+0.075*range(ys),'*','color','k','MarkerSize',11);
text(xs(1)+0.075*range(xs),ys(1)+0.075*range(ys),sprintf('p<%.3f',alpha),'FontSize',fontSize-2);

if isempty(fsic(varargin, 'noLabel'))
    xlabel('Interval name');
end
ylabel('Time interval change (ms)', 'FontSize', fontSize);

set(gca,'XLim',[0,numel(mean_int_chg.accel)+1]);
plot(xs,[0,0],'k-','LineWidth',1);
set(gca,'FontSize',fontSize);

if isempty(fsic(varargin, 'noLabel'))
    if isempty(fsic(varargin, 'upDown'))
        legend({'Accel - noPert','Decel - noPert'},'FontSize',fontSize,'Location','Southeast')
    else
        legend({'Down - noPert','Up - noPert'},'FontSize',fontSize,'Location','Southeast')
    end
end

title('Change in time interval (ms) (mean \pm SEM)','FontSize',fontSize);

set(gca,'XTick',[]);
for i1=1:6
    plot([i1,i1],ys(1)+[0,0.02*range(ys)],'k-'); hold on;
end
ys=get(gca,'YLim');

if isempty(fsic(varargin, 'noLabel'))
    text(1-0.2,ys(1)-0.05*range(ys),'[i]-[u]_1','FontSize',fontSize);
    text(2-0.2,ys(1)-0.05*range(ys),'[i]-[j]_1','FontSize',fontSize);
    text(3-0.22,ys(1)-0.05*range(ys),'[i]-[u]_2','FontSize',fontSize);
    text(4-0.22,ys(1)-0.05*range(ys),'[i]-[j]_2','FontSize',fontSize);
    text(5-0.22,ys(1)-0.05*range(ys),'[i]-[u]_3','FontSize',fontSize);
    text(6-0.22,ys(1)-0.05*range(ys),'[i]-[j]_3','FontSize',fontSize);
end
return
