function [vowel_acu, endPerts] = analyze_STUT_HA(subjID, runNums, varargin)
%%
hostName = getHostName;

if isequal(hostName, 'ubuntu')
    percDataDir=fullfile('/media/DATA/STUT_DATA',subjID,'APSTV2_STUT_EHPERC/DATA');
    percStimDir=fullfile('/media/DATA/STUT/DATA',subjID,'APSTV2_STUT_EHPERC/STIM');
else
    percDataDir=fullfile('E:\STUT_DATA\',subjID,'APSTV2_STUT_EHPERC\DATA');
    percStimDir=fullfile('E:\STUT_DATA\',subjID,'APSTV2_STUT_EHPERC\STIM');
end
stgNum=6;

colors={'r',[0,0.5,0],'b','k','m',[0.5,0,0],[0,0.5,0.5],[0.5,0,0.5],[0.5,0.5,0],[0.5,0.5,0.5],[0.25,0.25,0.25]};

%%
vowel_acu = NaN;

bPlot = 1;
if ~isempty(fsic(varargin, 'noPlot'))
    bPlot = 0;
end

if bPlot, figure; end
endPerts=[];
runCnt=1;
legendText={};
for k1=1:numel(runNums)
    runNum=runNums(k1);
    
    runDataDir=fullfile(percDataDir,[subjID,'_updown',num2str(runNum)]);
    if ~isdir(runDataDir)
        fprintf('Error: data for subject %s VID run %d not found.\n',subjID,runNum);
        return
    end

    d1=dir(fullfile(runDataDir,['session_info_',subjID,'*.mat']));

    load(fullfile(runDataDir,d1(1).name));
    
    for i1=1:numel(sess_info.stage(6).pert)
        stairPerts=sess_info.stage(6).pert{i1};
        if bPlot
            plot(1:numel(stairPerts),stairPerts,'o-','Color',colors{runCnt});
        end
        legendText{end+1}=sprintf('Staircase #%d',runCnt);
        runCnt=runCnt+1;
        if bPlot, hold on; end
        endPerts(end+1)=sess_info.stage(6).pert{i1}(end);
    end

end


if bPlot
    xlabel('Trial # in a run');
    ylabel('Pert (re. full pert)');
    legend(legendText,'Location','Northeast');
end

%%
meanEndPert = nanmean(endPerts);
medianEndPert = nanmedian(endPerts);
lastEndPert = endPerts(end);
firstEndPert = endPerts(1);
bestEndPert = nanmin(endPerts);

medianEndPert_last4 = nanmedian(endPerts(end - 3 : end));
meanEndPert_last4 = nanmean(endPerts(end - 3 : end));

if bPlot
    figure;
    plot(1:numel(endPerts),endPerts,'bo-');
    set(gca,'XLim',[0,numel(endPerts)+1]);
    hold on;
    
    plot([0,numel(endPerts)+1],repmat(meanEndPert,1,2),'ko--');
    xlabel('Staircase #');
    ylabel('Pert at end of staircase');

    legend({'Single-stairacase end perts','Mean end pert'});
end

%%
% vowel_acu = medianEndPert;
% vowel_acu = medianEndPert_last4;
vowel_acu = meanEndPert_last4;
% vowel_acu = mean(endPerts(4:5));
return