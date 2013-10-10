function test_APSTV2_makeFig_BLSC_fig
%%
% repDir='E:\APSTV2\TS_20100617_1\main\rep9'; 
repDir='E:\DATA\APSTV2\APSTV2_PS01\main\rep3';

mvaWinWidth=21;
fontSize=15;
namStyle='accelDecel';
% traceColor=[139/256,69/256,19/256];
traceColor=[1,1,1];
% traceColor=[0,0,0];
% traceColor_pert='c';
traceColor_pert='m';
traceWidth=3;
axesColor='k';
pointerColor='k';
pointerWidth=3;
bgColor='none';
timeBarColor=[0.9,0.9,0.9];
timeBarColor=[0, 0, 0];
textColor= [1, 1, 1];
focusIntervalLineClr = [1, 0, 1];

colors.none=[1, 1, 1];
colors.accel=[1,0.5,0];
colors.decel=[0,0.5,0];
colors.down='b';
colors.up='r';

% timeLim=[0.5,1.2];
timeLim_overview=[0.52,1.9];
timeLim=[0.72,1.0];

windowWidth_overview=1000;
windowHeight_overview=260;
windowWidth=300;
windowHeight=260;

lw=2;

scaling=2;

windowWidth_schem=(windowWidth*5-windowWidth_overview)/2;
windowHeight_schem=windowHeight_overview;

%% Find the accel perturbation
d1=dir(fullfile(repDir,'trial-*-1.mat'));
for i1=1:length(d1)
    load(fullfile(repDir,d1(i1).name));  % gives data
    if (data.params.bShift==1 && data.params.bTimeWarp==1 && data.params.bDecelWarp==0)
        dataFN_accel=fullfile(repDir,d1(i1).name);
        fprintf('Same-rep accel trial found: %s\n',dataFN_accel);
    elseif (data.params.bShift==1 && data.params.bTimeWarp==1 && data.params.bDecelWarp==1)
        dataFN_decel=fullfile(repDir,d1(i1).name);
        fprintf('Same-rep accel trial found: %s\n',dataFN_decel);
    elseif (data.params.bShift==1 && data.params.bTimeWarp==0 && ~isempty(data.params.pertField))
        if isequal(data.params.pertField.pertDirection,'up')
            dataFN_up=fullfile(repDir,d1(i1).name);
            fprintf('Same-rep up trial found: %s\n',dataFN_up);
        elseif isequal(data.params.pertField.pertDirection,'down')
            dataFN_down=fullfile(repDir,d1(i1).name);
            fprintf('Same-rep up trial found: %s\n',dataFN_down);
        end
    end
   
end

%%
for i1 = 1 : 5
    load(dataFN_decel);   % gives data

    p=data.params;
    p_decel=p;

    if i1==1
        p.bShift=0;
    elseif i1==2
        data_accel=load(dataFN_accel);
        p=data_accel.data.params;
        p_accel=p;
        
        % Manipulate the lambda - F1 map (p.lambdaF1Map)
%         p.lambdaF1Map=p.lambdaF1Map+100;
%         p.mapF1Max=p.mapF1Max+100;
%         p.twF1Max=p.twF1Max+100;
        d_maxF1=70;
        [maxF1,idxMaxF1]=max(p.lambdaF1Map);
        deltaF1=zeros(size(p.lambdaF1Map));
        deltaF1(1:idxMaxF1)=linspace(0,d_maxF1,idxMaxF1);
        deltaF1(idxMaxF1:end)=linspace(d_maxF1,0,length(deltaF1)-idxMaxF1+1);
        p.lambdaF1Map=p.lambdaF1Map+deltaF1;                
    elseif i1==3
        pause(0);
    elseif i1==4
        data_down=load(dataFN_down);
        p=data_down.data.params;
    elseif i1==5
        data_up=load(dataFN_up);
        p=data_up.data.params;
    end
   
    MexIO('init',p);
    MexIO('reset');    
    
    sigIn=resample(data.signalIn,48000,data.params.sr);
    sigInCell=makecell(sigIn,64);
    procTimes=nan(1,length(sigInCell));
    for n = 1 : length(sigInCell)
        tic;
        TransShiftMex(5,sigInCell{n});
        procTimes(n)=toc;
    end

    data=MexIO('getData');
    
    if i1==1
        repeats=2;
    else
        repeats=1;
    end
    
    for k1=1:repeats
        if i1==1 && k1==1
            figure('Position',[100,100,windowWidth_overview*scaling,windowHeight_overview*scaling]);
        else
            figure('Position',[100,100,windowWidth*scaling,windowHeight*scaling*0.75]);
        end
        set(gca,'FontSize',fontSize*scaling);

        sig=data.signalOut;    
        fs=data.params.sr;
%         [s,f,t]=spectrogram(sig, 64, 32, 1024, fs);
        [s,f,t]=spectrogram(sig, 128, 112, 1024, fs); % --- Use narrow-band spectogram to enhance frequency resolution --- %

        imagesc(t,f/1e3,10*log10(abs(s))); hold on;
        colormap gray;
        cmap=colormap;
        ncmap = size(cmap, 1);
%         cmap = [linspace(0, 0, ncmap)', ...
%                 (linspace(0, 1, ncmap)') .^ 4, ...
%                 linspace(0, 0, ncmap)'];
        cmap = [linspace(0, 0, ncmap)', ...
                (0.5 - 0.5 * cos(linspace(0, 1, ncmap)' * pi)) .^ 24, ...
                linspace(0, 0, ncmap)'];
            
%         cmap = 1 - repmat((0.5 - 0.5 * cos(linspace(0, 1, ncmap)' * pi)) .^ 5, 1, 3);
        Nshift = 2;
        cmap = [cmap(Nshift + 1: end, :); repmat(cmap(end, :), Nshift, 1)];
%         cmap = 1 - repmat(((0.5 - 0.5 * cos(linspace(0, 1, ncmap)' * pi)) .^ 5), 1, 3);        
        colormap(cmap);

        axis xy;
        hold on;

        frameDur=data.params.frameLen/data.params.sr;
        f1v=data.fmts(:,1)/1e3;
        f2v=data.fmts(:,2)/1e3;
        f1v=mva_nz(f1v,mvaWinWidth,'Hamming');
        f2v=mva_nz(f2v,mvaWinWidth,'Hamming');
        sentStat=data.sentStat;
        iouOnset=min(find(sentStat==2));
        youOnset=min(find(sentStat==4));
        taxis1=0:(frameDur):(frameDur*(length(f1v)-1));
        
        if i1==1
            f2v_noPert=f2v;
        end

        lineStyle='--';
        plot(taxis1,f1v,lineStyle,'LineWidth',traceWidth,'Color',traceColor);
        hold on;
        plot(taxis1,f2v,lineStyle,'LineWidth',traceWidth,'Color',traceColor);
        set(gca,'LineWidth',2,'XColor',axesColor,'YColor',axesColor);

        if p.bShift==1
            f1s=data.sfmts(:,1)/1e3;
            f2s=data.sfmts(:,2)/1e3;
            f1s=mva_nz(f1s,mvaWinWidth,'Hamming');
            f2s=mva_nz(f2s,mvaWinWidth,'Hamming');        
            
            if i1==4
                f2s_down=f2s;
            elseif i1==5
                f2s_up=f2s;
            end

            plot(taxis1(f2s>0),f2s(f2s>0),'--','LineWidth',traceWidth,'Color',traceColor_pert);
            plot(taxis1(f1s>0),f1s(f1s>0),'--','LineWidth',traceWidth,'Color',traceColor_pert);        
        end

        set(gca,'YLim',[0,2.5]);
        lim1=min(find(f1v>0));
        lim2=max(find(f1v>0));
        if i1==1 & k1==1
            set(gca,'XLim',timeLim_overview,'XTick',[]);
        else
            set(gca,'XLim',timeLim,'XTick',[]);
        end
        xlabel('Time');
        ylabel('Frequency (kHz)');

        
        if i1==1 && k1==1
            plot([0.6,0.7],[0.15,0.15],'-','LineWidth',traceWidth,'Color',timeBarColor);
            text(0.6,0.28,'100 ms','FontSize',(fontSize-3)*scaling*1.2,'FontWeight','Bold','Color',timeBarColor);
        else
            plot([0.75,0.85],[0.15,0.15],'-','LineWidth',traceWidth,'Color',timeBarColor);
            text(0.765,0.28,'100 ms','FontSize',(fontSize-3)*scaling*1.2,'FontWeight','Bold','Color',timeBarColor);
        end
        ys=get(gca,'YLim');
    
        if i1==1 || k1==1
            % Figure out [i] and [j]1 timing.
            df1v=diff(f1v); df2v=diff(f2v);
            df1v=mva(df1v,mvaWinWidth,'Hamming');
            df2v=mva(df2v,mvaWinWidth,'Hamming');
            idx_peaks=find(df2v(1:end-1)>0 & df2v(2:end)<0)+1;
            idx_troughs=find(df2v(1:end-1)<0 & df2v(2:end)>0)+1;
            
            t_a=taxis1(find(f2v>0,1));
            t_i=taxis1(idx_peaks(1));
            t_j1=taxis1(idx_peaks(2));
            t_j2=taxis1(idx_peaks(3));
            t_j3=taxis1(idx_peaks(4));
            t_u1=taxis1(idx_troughs(1));
            t_u2=taxis1(idx_troughs(2));
            t_u3=taxis1(idx_troughs(3));
            
            text(t_a,f2v(taxis1==t_a),'[a]','FontSize',fontSize*scaling*1.0,'Color',textColor);
            text(t_i,f2v(taxis1==t_i),'[i]','FontSize',fontSize*scaling*1.0,'Color',textColor);
            text(t_j1,f2v(taxis1==t_j1),'[j]_1','FontSize',fontSize*scaling*1.0,'Color',textColor);
            text(t_j2,f2v(taxis1==t_j2),'[j]_2','FontSize',fontSize*scaling*1.0,'Color',textColor);
            text(t_j3,f2v(taxis1==t_j3),'[j]_3','FontSize',fontSize*scaling*1.0,'Color',textColor);
            text(t_u1,f2v(taxis1==t_u1),'[u]_1','FontSize',fontSize*scaling*1.0,'Color',textColor);
            text(t_u2,f2v(taxis1==t_u2),'[u]_2','FontSize',fontSize*scaling*1.0,'Color',textColor);
            text(t_u3,f2v(taxis1==t_u3),'[u]_3','FontSize',fontSize*scaling*1.0,'Color',textColor);
%             plot(repmat(t_1,1,2),ys,'b-','LineWidth',2);
%             plot(repmat(t_2,1,2),ys,'b-','LineWidth',2);
            
            idx_downCross=find(f2v(1:end-1)>p.F2UB/1e3 & f2v(2:end)<p.F2UB/1e3);
            idx_upCross=find(f2v(1:end-1)<p.F2UB/1e3 & f2v(2:end)>p.F2UB/1e3);
            t_1=taxis1(idx_downCross(1));
            t_2=taxis1(idx_upCross(2));
            
            if i1==1 && k1==1
                plot(repmat(t_1,1,2),ys, '-', 'Color', focusIntervalLineClr, 'LineWidth',traceWidth);
                plot(repmat(t_2,1,2),ys, '-', 'Color', focusIntervalLineClr, 'LineWidth',traceWidth);
            end
        end
    end
    
%     if i1==1
%         title('A. No perturbation (noPert)','FontSize',fontSize+2,'FontWeight','Bold','Color',axesColor);
%     elseif i1==2
%         title('B. Accel perturbation','FontSize',fontSize+2,'FontWeight','Bold','Color',axesColor);
%     elseif i1==3
%         title('C. Decel perturbation','FontSize',fontSize+2,'FontWeight','Bold','Color',axesColor);
%     elseif i1==4
%         title('D. Down perturbation','FontSize',fontSize+2,'FontWeight','Bold','Color',axesColor);
%     elseif i1==5
%         title('E. Up perturbation','FontSize',fontSize+2,'FontWeight','Bold','Color',axesColor);        
%     end   
    
    if i1>1 && k1==1   
        legend_left=timeLim(1)+0.01;
        legend_bottom=2.00;
        legend_width=0.175;
        legend_height=0.50;
        rectangle('Position',[legend_left,legend_bottom,legend_width,legend_height],'EdgeColor','k','FaceColor',[0.25,0.25,0.25])
        plot([legend_left+0.01,legend_left+0.05],[2.34,2.34],'--','LineWidth',traceWidth,'Color',traceColor);
        plot([legend_left+0.01,legend_left+0.05],[2.12,2.12],'--','LineWidth',traceWidth,'Color',traceColor_pert);
        text(legend_left+0.055,2.35,'Mic. input','Color',traceColor,'FontSize',(fontSize-5)*scaling);
        text(legend_left+0.055,2.14,'Perturbed AF','Color',traceColor_pert,'FontSize',(fontSize-4.5)*scaling);
    end
    

end

%% Schematic drawing 1: mapping of F2 in Down and Up
figure('Position',[100,80,windowWidth_schem*scaling,windowHeight_schem*scaling*0.8]);
set(gca,'FontSize',(fontSize-1)*scaling);

f2s_nz=f2s_down(f2s_down>0);
f2s_down_min=min(f2s_nz);
idx_ctr=find(f2s_down==f2s_down_min);

f2s_down(f2s_down==0 & f2v_noPert>0)=f2v_noPert(f2s_down==0 & f2v_noPert>0);
f2s_up(f2s_up==0 & f2v_noPert>0)=f2v_noPert(f2s_up==0 & f2v_noPert>0);

plot(taxis1,f2s_down,'Color',colors.down,'LineWidth',lw);
hold on;
plot(taxis1,f2s_up,'Color',colors.up,'LineWidth',lw);
plot(taxis1,f2v_noPert,'Color',colors.none,'LineWidth',lw);

set(gca,'XLim',[timeLim(1)+0.125*range(timeLim),timeLim(2)-0.125*range(timeLim)],'XTick',[]);
set(gca,'YLim',[0.9,p.F2UB*1.05/1e3],'YTick',[]);
xlabel('Time','FontSize',fontSize*scaling);
ylabel('Frequency','FontSize',fontSize*scaling);

arrow([taxis1(idx_ctr),f2v_noPert(idx_ctr)],[taxis1(idx_ctr),f2s_down(idx_ctr)],'width',2,'Length',19,...
    'FaceColor',colors.down,'EdgeColor',colors.down);
text(taxis1(idx_ctr),mean([f2v_noPert(idx_ctr),f2s_down(idx_ctr)]),'\DeltaF_2',...
    'FontSize',fontSize*scaling*0.8,'FontName','Times New Roman','FontAngle','italic','Color',colors.down);
arrow([taxis1(idx_ctr),f2v_noPert(idx_ctr)],[taxis1(idx_ctr),f2s_up(idx_ctr)],'width',2,'Length',19,...
    'FaceColor',colors.up,'EdgeColor',colors.up);
text(taxis1(idx_ctr),mean([f2v_noPert(idx_ctr),f2s_up(idx_ctr)]),'\DeltaF_2',...
    'FontSize',fontSize*scaling*0.8,'FontName','Times New Roman','FontAngle','italic','Color',colors.up);

legend({'Down','Up','Baseline'},'Location','Southwest','FontSize',fontSize*scaling*0.6);

xs=get(gca,'XLim');
plot(xs,repmat(p.F2UB/1e3,1,2),'k--','LineWidth',lw);
text(xs(1)+0.05*range(xs),p.F2UB/1e3*1.05,'F_2^m^a^x','FontSize',fontSize*scaling*0.8,...
    'FontName','Times New Roman','FontAngle','italic');
% set(gca,'XLim',[timeLim(1)+0.05*range(timeLim),timeLim(2)-0.05*range(timeLim)]);
box off;

%% Schematic drawing 2: time-warping in Accel and Decel 
figure('Position',[100,80,windowWidth_schem*scaling,windowHeight_schem*scaling]);
set(gca,'FontSize',(fontSize-1)*scaling);
x=linspace(0,1,numel(p_accel.twLambdaLUT));
% plot(x,x.*p_accel.twLambdaLUT,'Color',colors.accel);
plot(x,p_accel.twLambdaLUT,'Color',colors.accel,'LineWidth',lw);
hold on;
plot(x,p_decel.twLambdaLUT,'Color',colors.decel,'LineWidth',lw);
plot([x(1),x(end)],[x(1),x(end)],'--','Color',[0.5,0.5,0.5],'LineWidth',lw);
axis equal;
set(gca,'XLim',[x(1),x(end)],'YLim',[x(1),x(end)]);
set(gca,'XTick',[0:0.2:1],'YTick',[0:0.2:1]);
xlabel('Orignal normalized time','FontSize',fontSize*scaling*0.85);
ylabel('Shifted normalized time','FontSize',fontSize*scaling*0.85);
legend({'Accel','Decel'},'Location','Southeast','FontSize',fontSize*scaling*0.75);
return