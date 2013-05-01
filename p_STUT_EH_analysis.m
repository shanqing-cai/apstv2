function p_STUT_EH_analysis(baseField, varargin)
%% CONFIG
subjIDs.PFS = get_STUT_EH_subjIDs('PFS');
subjIDs.PWS = get_STUT_EH_subjIDs('PWS');

% subjIDs.PFS = subjIDs.PFS(1 : end - 1);

FRAME_DUR=16/12000;

colors.down= [0, 0, 1];
colors.up= [1, 0, 0];
colors.down_neg = 0.5 + 0.5 * colors.down;
colors.up_neg = 0.5 + 0.5 * colors.up;

colors.PFS = 'k';
colors.PWS = [160, 32, 240] / 255;

N_SAMP = 225;   % Default: 165 * 1.333 = 220 (ms)
if ~isempty(fsic(varargin, 'N_SAMP'))
    N_SAMP = varargin{fsic(varargin, 'N_SAMP') + 1};
    fprintf('INFO: N_SAMP set to %d\n', N_SAMP);
end

% QUANTIFY_TIMES = {0.15, 0.185, 0.22, [0.18, 0.21]};  % s
QUANTIFY_TIME_ADAPT = [0.03, 0.10];  % s
QUANTIFY_TIMES = {[0.30, 0.30]};  % s. Old value: {0.15, 0.185, [0.22, 0.22]}

xlsSysDir = 'e:\speechres\stut\xlsSysEH';
fontSize = 12;

FDA_dataSet_FN = 'e:\speechres\stut\mcode\p_BBAP_FDA_analysis_ds.mat';

FDR = 0.05;

MAX_LAG = 3;

downFld = 'down';
upFld = 'up';

% baseField = 'none';
if ~exist('baseField')
    error('baseField not specified.')
end

F1_SD_QUANT_TIME = [0.298, 0.30];    % Unit of time: s

%% Settings related to the mixed ANOVA
% mAOV_tPoints = {[0.00, 0.00], [0.05, 0.05], [0.10, 0.10], [0.15, 0.15], ...
%                 [0.20, 0.20], [0.25, 0.25], [0.30, 0.30]};
% mAOV_tPoints = {[0.00, 0.00], [0.025, 0.025], [0.05, 0.05], [0.075, 0.075], ...
%                 [0.10, 0.10], [0.125, 0.125], [0.15, 0.15], [0.175, 0.175], ...
%                 [0.20, 0.20], [0.225, 0.225], [0.25, 0.25], [0.275, 0.275], [0.30, 0.30]};
% mAOV_tPoints = {[0.00, 0.020], [0.020, 0.040], [0.040, 0.060], [0.060, 0.080], ...
%                 [0.080, 0.100], [0.100, 0.120], [0.120, 0.140], [0.140, 0.160], ...
%                 [0.160, 0.180], [0.180, 0.200], [0.200, 0.220], [0.220, 0.240], ...
%                 [0.240, 0.260], [0.260, 0.280], [0.280, 0.300]};
% mAOV_tPoints = {[0.000, 0.015], [0.015, 0.030], [0.030, 0.045], [0.045, 0.060], ...
%                 [0.060, 0.075], [0.075, 0.090], [0.090, 0.105], [0.105, 0.120], ...
%                 [0.120, 0.135], [0.135, 0.150], [0.150, 0.165], [0.165, 0.180], ...
%                 [0.180, 0.195], [0.195, 0.210], [0.210, 0.225], [0.225, 0.240], ...
%                 [0.240, 0.255], [0.255, 0.270], [0.270, 0.285], [0.285, 0.300]};
% mAOV_tPoints = {[0.00, 0.025], [0.025, 0.05], [0.05, 0.075], [0.075, 0.1], ...
%                 [0.10, 0.125], [0.125, 0.15], [0.15, 0.175], [0.175, 0.2], ...
%                 [0.20, 0.225], [0.225, 0.25], [0.25, 0.275], [0.275, 0.3]};
mAOV_tPoints = {[0.000, 0.000], [0.030, 0.030], [0.060, 0.060], [0.090, 0.090], ...
                [0.120, 0.120], [0.150, 0.150], [0.180, 0.180], [0.210, 0.210], ...
                [0.240, 0.240], [0.270, 0.270], [0.300, 0.300]};

if ~isempty(fsic(varargin, '1'))
    mAOV_xlsFN = fullfile(cdds, 'stut\xlsSysEH', 'STUT_EH_mAOV.xls');
    mAOV_downUp_xlsFN= fullfile(cdds, 'stut\xlsSysEH', 'STUT_EH_mAOV_downUp.xls');
else
    mAOV_xlsFN = fullfile(cdds, 'stut\xlsSysEH', 'STUT_EH_mAOV_fullDS.xls');
    mAOV_downUp_xlsFN = fullfile(cdds, 'stut\xlsSysEH', 'STUT_EH_mAOV_downUp_fullDS.xls');
end

spacingInt = [0.00, 0.050]; % Unit: sec
            
%% Obtain the hearing acuity data

if ~isempty(fsic(varargin, '1'))
    downFld = [downFld, '1'];
    downFld_neg = [downFld, 'neg'];
    upFld = [upFld, '1'];
    upFld_neg = [upFld, 'neg'];
end

groups = fields(subjIDs);
for i1 = 1 : numel(groups)
    grp = groups{i1};
    vowelAcuity.(grp) = nan(1, numel(subjIDs.(grp)));
    
    for i2 = 1 : numel(subjIDs.(grp))
        t_subjID = subjIDs.(grp){i2};
        if isequal(t_subjID(end-1:end), '_1')
            t_subjID = t_subjID(1:end-2);
        end
        vowelAcuity.(grp)(i2) = analyze_STUT_HA(t_subjID,[1,2],'noPlot');
    end
end

%%
hostName=getHostName;
if isequal(hostName,'smcg-w510') || isequal(hostName,'smcgw510') || isequal(hostName,'smcg_w510')
    dacacheDir='e:/speechres/apstv2/mcode/dacache_STUT';
%     rawDataDir='e:/STUT_DATA';
elseif isequal(hostName,'glossa')
    dacacheDir='z:/speechres/apstv2/mcode/dacache_STUT';
%     rawDataDir='f:/STUT_DATA';
elseif isequal(hostName, 'ubuntu')
    dacacheDir='/media/DATA/speechres/apstv2/mcode/dacache_STUT';   
else
    dacacheDir='z:/speechres/apstv2/mcode/dacache_STUT';
%     rawDataDir='z:/STUT_DATA/';
end

upDownContra.PWS = nan(numel(subjIDs.PWS), numel(QUANTIFY_TIMES));
upDownContra.PFS = nan(numel(subjIDs.PFS), numel(QUANTIFY_TIMES));
sdF1.PWS.none = nan(numel(subjIDs.PWS), numel(QUANTIFY_TIMES));
sdF1.PWS.up = nan(numel(subjIDs.PWS), numel(QUANTIFY_TIMES));
sdF1.PWS.down = nan(numel(subjIDs.PWS), numel(QUANTIFY_TIMES));
sdF1.PFS.none = nan(numel(subjIDs.PFS), numel(QUANTIFY_TIMES));
sdF1.PFS.up = nan(numel(subjIDs.PFS), numel(QUANTIFY_TIMES));
sdF1.PFS.down = nan(numel(subjIDs.PFS), numel(QUANTIFY_TIMES));
sdRatio.PWS = nan(numel(subjIDs.PWS), numel(QUANTIFY_TIMES));
sdRatio.PFS = nan(numel(subjIDs.PFS), numel(QUANTIFY_TIMES));

bpRMS.downUp.PWS = nan(numel(subjIDs.PWS), 1);
bpRMS.downUp.PFS = nan(numel(subjIDs.PFS), 1);

aftUpDownContra.PWS = nan(numel(subjIDs.PWS), 1);
aftUpDownContra.PFS = nan(numel(subjIDs.PFS), 1);

% a_vowelOnset_manCorr.all.PWS = [];
% a_vowelOnset_manCorr.all.PFS = [];
% a_vowelOnset_manCorr.pert.PWS = [];
% a_vowelOnset_manCorr.pert.PFS = [];

SSI4.total.PWS = nan(numel(subjIDs.PWS), 1);
SSI4.freq.PWS = nan(numel(subjIDs.PWS), 1);
SSI4.dur.PWS = nan(numel(subjIDs.PWS), 1);
SSI4.concom.PWS = nan(numel(subjIDs.PWS), 1);

figure('Position',[100,100,720,320]);
set(gca, 'FontSize', fontSize);
a_avgChg_traj_F1_ML.PWS.down = nan(N_SAMP,0);
a_avgChg_traj_F1_ML.PWS.up = nan(N_SAMP,0);
a_avgChg_traj_F1_ML.PFS.down = nan(N_SAMP,0);
a_avgChg_traj_F1_ML.PFS.up = nan(N_SAMP,0);

a_avgChg_traj_F1_ns_ML.PWS.down = nan(N_SAMP,0);    % Not smoothed
a_avgChg_traj_F1_ns_ML.PWS.up = nan(N_SAMP,0);
a_avgChg_traj_F1_ns_ML.PFS.down = nan(N_SAMP,0);
a_avgChg_traj_F1_ns_ML.PFS.up = nan(N_SAMP,0);

if isequal(baseField, 'none0')
    a_avgF1Traj_none0neg_ML.PFS.none0 = nan(N_SAMP, 0);
    a_avgF1Traj_none0neg_ML.PFS.down = nan(N_SAMP, 0);
    a_avgF1Traj_none0neg_ML.PFS.up = nan(N_SAMP, 0);
    a_avgF1Traj_none0neg_ML.PWS.none0 = nan(N_SAMP, 0);
    a_avgF1Traj_none0neg_ML.PWS.down = nan(N_SAMP, 0);
    a_avgF1Traj_none0neg_ML.PWS.up = nan(N_SAMP, 0);
elseif isequal(baseField, 'none1')
    a_avgF1Traj_none1neg_ML.PFS.none1 = nan(N_SAMP, 0);
    a_avgF1Traj_none1neg_ML.PFS.down = nan(N_SAMP, 0);
    a_avgF1Traj_none1neg_ML.PFS.up = nan(N_SAMP, 0);
    a_avgF1Traj_none1neg_ML.PWS.none1 = nan(N_SAMP, 0);
    a_avgF1Traj_none1neg_ML.PWS.down = nan(N_SAMP, 0);
    a_avgF1Traj_none1neg_ML.PWS.up = nan(N_SAMP, 0);
end

if ~isempty(fsic(varargin, '1'))
    a_avgChg_traj_F1_ML.PWS.down_neg = nan(N_SAMP,0);
    a_avgChg_traj_F1_ML.PWS.up_neg = nan(N_SAMP,0);
    a_avgChg_traj_F1_ML.PFS.down_neg = nan(N_SAMP,0);
    a_avgChg_traj_F1_ML.PFS.up_neg = nan(N_SAMP,0);
end

a_avgChg_traj_F1_ML.PWS.downUp = nan(N_SAMP, 0);
a_avgChg_traj_F1_ML.PFS.downUp = nan(N_SAMP, 0);

a_avgChg_traj_F1_ns_ML.PWS.downUp = nan(N_SAMP, 0);
a_avgChg_traj_F1_ns_ML.PFS.downUp = nan(N_SAMP, 0);

a_sd_F1_ML.PWS.none = nan(N_SAMP,0);
a_sd_F1_ML.PWS.none0 = nan(N_SAMP,0);
a_sd_F1_ML.PWS.none1 = nan(N_SAMP,0);
a_sd_F1_ML.PWS.down = nan(N_SAMP,0);
a_sd_F1_ML.PWS.up = nan(N_SAMP,0);
a_sd_F1_ML.PFS.none = nan(N_SAMP,0);
a_sd_F1_ML.PFS.none0 = nan(N_SAMP,0);
a_sd_F1_ML.PFS.none1 = nan(N_SAMP,0);
a_sd_F1_ML.PFS.down = nan(N_SAMP,0);
a_sd_F1_ML.PFS.up = nan(N_SAMP,0);

a_sdRatio_F1_ML.PWS.down = nan(N_SAMP, 0);
a_sdRatio_F1_ML.PWS.up = nan(N_SAMP, 0);
a_sdRatio_F1_ML.PFS.down = nan(N_SAMP, 0);
a_sdRatio_F1_ML.PFS.up = nan(N_SAMP, 0);

adapt_avgChg_F1_ML.PWS.aftDown = nan(N_SAMP, 0);
adapt_avgChg_F1_ML.PWS.aftUp = nan(N_SAMP, 0);
adapt_avgChg_F1_ML.PFS.aftDown = nan(N_SAMP, 0);
adapt_avgChg_F1_ML.PFS.aftUp = nan(N_SAMP, 0);

a_F1_downPert.PWS = [];
a_F1_downPert.PFS = [];
a_F1_upPert.PWS = [];
a_F1_upPert.PFS = [];

a_lat.PWS = [];
a_lat.PFS = [];

a_compenSlope.PWS = [];
a_compenSlope.PFS = [];

a_p_quadVline.PWS = [];
a_p_quadVline.PFS = [];
a_R2_lineSpline.PWS = [];
a_R2_lineSpline.PFS = [];

lagN_F1Chg.PWS.down = nan(0, MAX_LAG);
lagN_F1Chg.PWS.up = nan(0, MAX_LAG);
lagN_F1Chg.PFS.down = nan(0, MAX_LAG);
lagN_F1Chg.PFS.up = nan(0, MAX_LAG);
earlyFlw_downUp.PWS.down = [];
earlyFlw_downUp.PWS.up = [];
earlyFlw_downUp.PFS.down = [];
earlyFlw_downUp.PFS.up = [];

ratioKeptTrials.PWS.down = [];  ratioKeptTrials.PWS.up = [];
ratioKeptTrials.PFS.down = [];  ratioKeptTrials.PFS.up = [];

ratioLenDiscard.PWS.none = []; ratioLenDiscard.PWS.down = []; ratioLenDiscard.PWS.up = [];
ratioLenDiscard.PFS.none = []; ratioLenDiscard.PFS.down = []; ratioLenDiscard.PFS.up = [];

a_CD_thresh = [];

a_nDscd_prodErr.PWS = [];   a_nDscd_prodErr.PFS = [];
a_nTotTrials.PWS = [];      a_nTotTrials.PFS = [];
a_nDscd_all.PWS = [];       a_nDscd_all.PFS = [];

a_avgDur.PFS = [];  a_avgDur.PWS = [];  % Columns: [none, down, up]
a_avgLv.PFS = [];   a_avgLv.PWS = [];   % Columns: [none, down, up]

spacing_chgF1_row1 = {'SNUM', 'GRP', 'DOWNL', 'DOWNS', 'UPL', 'UPS'};
spacing_chgF1_tab = nan(0, 6); % A row: [SNUM, GRP, DOWN_LARGE_SPAC, DOWN_SMALL_SPAC, UP_LARGE_SPACE, DOWN_SMALL_SPAC]


propComp.down.PWS = [];   propComp.down.PFS = [];
propComp.up.PWS = [];   propComp.up.PFS = [];
propFoll.down.PWS = [];   propFoll.down.PFS = [];
propFoll.up.PWS = [];   propFoll.up.PFS = [];

groups=fields(subjIDs);
for i0=1:numel(groups)
    subplot('Position',[0.14+(i0-1)*0.42,0.12,0.42,0.8]);
    set(gca, 'FontSize', fontSize);
    grp=groups{i0};
    
    for i1=1:numel(subjIDs.(grp))
        load(fullfile(dacacheDir,[subjIDs.(grp){i1},'_EH.mat']));   % gives pdata
        
        [nDscd_prodErr, nTotTrials, nDscd_all] = get_EH_discard_stats(pdata);
        a_nDscd_prodErr.(grp)(end + 1) = nDscd_prodErr;
        a_nTotTrials.(grp)(end + 1) = nTotTrials;
        a_nDscd_all.(grp)(end + 1) = nDscd_all;
        
        ratioLenDiscard.(grp).none(end + 1) = pdata.stage2.lenDiscardRatio.none;
        ratioLenDiscard.(grp).down(end + 1) = pdata.stage2.lenDiscardRatio.down;
        ratioLenDiscard.(grp).up(end + 1) = pdata.stage2.lenDiscardRatio.up;
        
        lagN_F1Chg.(grp).down = [lagN_F1Chg.(grp).down; ...
            pdata.stage2.avg_lagF1Win.down(2 : MAX_LAG + 1) - pdata.stage2.avg_lagF1Win.down(1)];
        lagN_F1Chg.(grp).up = [lagN_F1Chg.(grp).up; ...
            pdata.stage2.avg_lagF1Win.up(2 : MAX_LAG + 1) - pdata.stage2.avg_lagF1Win.up(1)];
        
        earlyFlw_downUp.(grp).down(end + 1) =  pdata.stage2.earlyF1Chg.down;
        earlyFlw_downUp.(grp).up(end + 1) =  pdata.stage2.earlyF1Chg.up;
        
        t_chg = (pdata.stage2.avg_traj_F1_ML.(downFld)(:,1)-pdata.stage2.avg_traj_F1_ML.(baseField)(:,1))./...
                 pdata.stage2.avg_traj_F1_ML.(baseField)(:,1);
        t_chg_ns = (pdata.stage2.avg_traj_F1_ns_ML.(downFld)(:,1) - pdata.stage2.avg_traj_F1_ns_ML.(baseField)(:,1))./...
                    pdata.stage2.avg_traj_F1_ns_ML.(baseField)(:,1);
             
             
        a_avgChg_traj_F1_ML.(grp).down = [a_avgChg_traj_F1_ML.(grp).down, t_chg(1:N_SAMP,:)];
        if ~isequal(downFld, 'down')
            t_chg_neg = (pdata.stage2.avg_traj_F1_ML.(downFld_neg)(:,1)-pdata.stage2.avg_traj_F1_ML.(baseField)(:,1))./...
                         pdata.stage2.avg_traj_F1_ML.(baseField)(:,1);
            a_avgChg_traj_F1_ML.(grp).down_neg = [a_avgChg_traj_F1_ML.(grp).down_neg, t_chg_neg];
        end
        
        a_avgChg_traj_F1_ns_ML.(grp).down = [a_avgChg_traj_F1_ns_ML.(grp).down, t_chg_ns(1:N_SAMP,:)];        
        
        ratioKeptTrials.(grp).down(end + 1) = pdata.stage2.avg_traj_F1_ML.(downFld)(1, 3) / pdata.stage2.avg_traj_F1_ML.down(1, 3);
        ratioKeptTrials.(grp).up(end + 1) = pdata.stage2.avg_traj_F1_ML.(upFld)(1, 3) / pdata.stage2.avg_traj_F1_ML.up(1, 3);
        
        t_chg = (pdata.stage2.avg_traj_F1_ML.(upFld)(:,1)-pdata.stage2.avg_traj_F1_ML.(baseField)(:,1))./...
                 pdata.stage2.avg_traj_F1_ML.(baseField)(:,1);
        t_chg_ns = (pdata.stage2.avg_traj_F1_ns_ML.(upFld)(:,1)-pdata.stage2.avg_traj_F1_ns_ML.(baseField)(:,1))./...
                    pdata.stage2.avg_traj_F1_ns_ML.(baseField)(:,1);
             
        a_avgChg_traj_F1_ML.(grp).up=[a_avgChg_traj_F1_ML.(grp).up, t_chg(1:N_SAMP,:)];
        if ~isequal(upFld, 'up')
            t_chg_neg = (pdata.stage2.avg_traj_F1_ML.(upFld_neg)(:,1)-pdata.stage2.avg_traj_F1_ML.(baseField)(:,1))./...
                         pdata.stage2.avg_traj_F1_ML.(baseField)(:,1);
            a_avgChg_traj_F1_ML.(grp).up_neg = [a_avgChg_traj_F1_ML.(grp).up_neg, t_chg_neg];
        end
        
        a_avgChg_traj_F1_ns_ML.(grp).up = [a_avgChg_traj_F1_ns_ML.(grp).up, t_chg_ns(1:N_SAMP,:)];
        
        a_avgChg_traj_F1_ML.(grp).downUp = [a_avgChg_traj_F1_ML.(grp).downUp, a_avgChg_traj_F1_ML.(grp).down(:, end) - a_avgChg_traj_F1_ML.(grp).up(:, end)];
        a_avgChg_traj_F1_ns_ML.(grp).downUp = [a_avgChg_traj_F1_ns_ML.(grp).downUp, a_avgChg_traj_F1_ns_ML.(grp).down(:, end) - a_avgChg_traj_F1_ns_ML.(grp).up(:, end)];
        
        a_sd_F1_ML.(grp).(baseField) = [a_sd_F1_ML.(grp).(baseField), pdata.stage2.avg_traj_F1_ML.(baseField)(:,2)];
        a_sd_F1_ML.(grp).down = [a_sd_F1_ML.(grp).down, pdata.stage2.avg_traj_F1_ML.(downFld)(:,2)];
        a_sd_F1_ML.(grp).up = [a_sd_F1_ML.(grp).up, pdata.stage2.avg_traj_F1_ML.(upFld)(:,2)];
        
        a_sdRatio_F1_ML.(grp).down = [a_sdRatio_F1_ML.(grp).down, ...
            pdata.stage2.avg_traj_F1_ML.(downFld)(:,2) ./ pdata.stage2.avg_traj_F1_ML.(baseField)(:,2)];
        a_sdRatio_F1_ML.(grp).up = [a_sdRatio_F1_ML.(grp).up, ...
            pdata.stage2.avg_traj_F1_ML.(upFld)(:,2) ./ pdata.stage2.avg_traj_F1_ML.(baseField)(:,2)];
        
        len = size(pdata.stage2.avg_traj_F1_ML.(downFld), 1);
        adapt_chg = (pdata.stage2.adapt_avg_traj_F1.aftDown(1:len,1)-pdata.stage2.adapt_avg_traj_F1.aftNone(1:len,1))./...
            pdata.stage2.adapt_avg_traj_F1.aftNone(1:len,1);
        adapt_avgChg_F1_ML.(grp).aftDown = [adapt_avgChg_F1_ML.(grp).aftDown, adapt_chg(1:N_SAMP,:)];
        
        adapt_chg = (pdata.stage2.adapt_avg_traj_F1.aftUp(1:len,1)-pdata.stage2.adapt_avg_traj_F1.aftNone(1:len,1))./...
            pdata.stage2.adapt_avg_traj_F1.aftNone(1:len,1);
        adapt_avgChg_F1_ML.(grp).aftUp = [adapt_avgChg_F1_ML.(grp).aftUp, adapt_chg(1:N_SAMP,:)];
        
        taxis0 = 1e3 * (0:FRAME_DUR:FRAME_DUR*(size(a_avgChg_traj_F1_ML.(grp).down,1)-1));
        plot(taxis0,a_avgChg_traj_F1_ML.(grp).down(:,end),'color',colors.down);
        hold on;
        plot(taxis0,a_avgChg_traj_F1_ML.(grp).up(:,end),'color',colors.up);        
        
        if isequal(baseField, 'none0')
            a_avgF1Traj_none0neg_ML.(grp).none0 = [a_avgF1Traj_none0neg_ML.(grp).none0, ...
                pdata.stage2.avg_traj_F1_ML.none0(:, 1)];
            a_avgF1Traj_none0neg_ML.(grp).down = [a_avgF1Traj_none0neg_ML.(grp).down, ...
                pdata.stage2.avg_traj_F1_ML.none0neg_down(:, 1)];
            a_avgF1Traj_none0neg_ML.(grp).up = [a_avgF1Traj_none0neg_ML.(grp).up, ...
                pdata.stage2.avg_traj_F1_ML.none0neg_up(:, 1)];
        elseif isequal(baseField, 'none1')
            a_avgF1Traj_none1neg_ML.(grp).none1 = [a_avgF1Traj_none1neg_ML.(grp).none1, ...
                pdata.stage2.avg_traj_F1_ML.none1(:, 1)];
            a_avgF1Traj_none1neg_ML.(grp).down = [a_avgF1Traj_none1neg_ML.(grp).down, ...
                pdata.stage2.avg_traj_F1_ML.none1neg_down(:, 1)];
            a_avgF1Traj_none1neg_ML.(grp).up = [a_avgF1Traj_none1neg_ML.(grp).up, ...
                pdata.stage2.avg_traj_F1_ML.none1neg_up(:, 1)];
        end
        
        if isequal(grp, 'PWS')
            SSI4.total.PWS(i1) = get_PWS_SSI4(strrep(subjIDs.(grp){i1}, '_1', ''));
            SSI4.freq.PWS(i1) = get_PWS_SSI4(strrep(subjIDs.(grp){i1}, '_1', ''), 'freq');
            SSI4.dur.PWS(i1) = get_PWS_SSI4(strrep(subjIDs.(grp){i1}, '_1', ''), 'dur');
            SSI4.concom.PWS(i1) = get_PWS_SSI4(strrep(subjIDs.(grp){i1}, '_1', ''), 'concom');
        end
        
        a_avgDur.(grp) = [a_avgDur.(grp); ...
                          [pdata.stage2.avg_vowelDur.none, pdata.stage2.avg_vowelDur.down, pdata.stage2.avg_vowelDur.up]];
        a_avgLv.(grp) = [a_avgLv.(grp); ...
                         [pdata.stage2.avg_vowelLv.none, pdata.stage2.avg_vowelLv.down, pdata.stage2.avg_vowelLv.up]];             
        
        for k1 = 1 : numel(QUANTIFY_TIMES)
            if numel(QUANTIFY_TIMES{k1}) == 1
                idx = round(QUANTIFY_TIMES{k1} / FRAME_DUR);
                upDownContra.(grp)(i1,k1) = a_avgChg_traj_F1_ML.(grp).down(idx, end) - a_avgChg_traj_F1_ML.(grp).up(idx, end);
                
                sdF1.(grp).(baseField)(i1,k1) = pdata.stage2.avg_traj_F1_ML.(baseField)(idx, 2);
                sdF1.(grp).up(i1,k1) = pdata.stage2.avg_traj_F1_ML.(upFld)(idx, 2);
                sdF1.(grp).down(i1,k1) = pdata.stage2.avg_traj_F1_ML.(downFld)(idx, 2);
                
                sdRatio.(grp)(i1,k1) = rms([sdF1.(grp).up(i1, k1), sdF1.(grp).down(i1, k1)]) / sdF1.(grp).(baseField)(i1, k1);
            elseif numel(QUANTIFY_TIMES{k1}) == 2
                idx1 = round(QUANTIFY_TIMES{k1}(1) / FRAME_DUR);
                idx2 = round(QUANTIFY_TIMES{k1}(2) / FRAME_DUR);
                upDownContra.(grp)(i1,k1) = mean(a_avgChg_traj_F1_ML.(grp).down(idx1 : idx2, end) - a_avgChg_traj_F1_ML.(grp).up(idx1 : idx2, end));
                
                sdF1.(grp).(baseField)(i1,k1) = rms(pdata.stage2.avg_traj_F1_ML.(baseField)(idx1 : idx2, 2));
                sdF1.(grp).up(i1,k1) = rms(pdata.stage2.avg_traj_F1_ML.(upFld)(idx1 : idx2, 2));
                sdF1.(grp).down(i1,k1) = rms(pdata.stage2.avg_traj_F1_ML.(downFld)(idx1 : idx2, 2));
                
                sdRatio.(grp)(i1,k1) = rms([sdF1.(grp).up(i1,k1), sdF1.(grp).down(i1,k1)]) / sdF1.(grp).(baseField)(i1, k1);
            end
        end
        
        idx1 = round(QUANTIFY_TIME_ADAPT(1) / FRAME_DUR);
        idx2 = round(QUANTIFY_TIME_ADAPT(2) / FRAME_DUR);
        aftUpDownContra.(grp)(i1) = mean(adapt_avgChg_F1_ML.(grp).aftDown(idx1 : idx2, end) - adapt_avgChg_F1_ML.(grp).aftUp(idx1 : idx2, end));
        
%         a_vowelOnset_manCorr.all.(grp) = [a_vowelOnset_manCorr.all.(grp), ...
%             [pdata.stage2.a_vowelOnset_manCorr.none, pdata.stage2.a_vowelOnset_manCorr.up, pdata.stage2.a_vowelOnset_manCorr.down]];
%         a_vowelOnset_manCorr.pert.(grp) = [a_vowelOnset_manCorr.all.(grp), ...
%             [pdata.stage2.a_vowelOnset_manCorr.up, pdata.stage2.a_vowelOnset_manCorr.down]];
        
        a_F1_downPert.(grp)(end + 1) = mean(pdata.stage2.avg_traj_F1.(downFld)(:, 1)) * pdata.subject.pertRatio;
        a_F1_upPert.(grp)(end + 1) = mean(pdata.stage2.avg_traj_F1.(upFld)(:, 1)) * pdata.subject.pertRatio;
        
        if ~isempty(fsic(varargin, '1'))
%             latFld = 'lat_subsetMode';
%             latFld = 'lat_noSD_subsetMode';
            latFld = 'lat_fit_subsetMode';
            if ~isempty(fsic(varargin, 'linTurnPoint'))
                latFld = 'lat_fit_linTurnPoint_subsetMode';
            end
            
            compenSlopeFld = 'compenSlope_subsetMode';
            SSreg_lineFld = 'SSreg_lineSpline_subsetMode';
            SSreg_quadFld = 'SSreg_quadSpline_subsetMode';
            df_lineFld = 'df_lineSpline_subsetMode';
            df_quadFld = 'df_quadSpline_subsetMode';
            R2_lineFld = 'R2_lineSpline_subsetMode';
        else
%             latFld = 'lat_nonSubsetMode';
%             latFld = 'lat_noSD_nonSubsetMode';
            latFld = 'lat_fit_nonSubsetMode';
            if ~isempty(fsic(varargin, 'linTurnPoint'))
                latFld = 'lat_fit_linTurnPoint_nonSubsetMode';
            end
            
            compenSlopeFld = 'compenSlope_nonSubsetMode';
            SSreg_lineFld = 'SSreg_lineSpline_nonSubsetMode';
            SSreg_quadFld = 'SSreg_quadSpline_nonSubsetMode';
            df_lineFld = 'df_lineSpline_nonSubsetMode';
            df_quadFld = 'df_quadSpline_nonSubsetMode';
            R2_lineFld = 'R2_lineSpline_nonSubsetMode';
        end
        
        a_lat.(grp)(end + 1) = pdata.stage2.(latFld);
        a_compenSlope.(grp)(end + 1) = pdata.stage2.(compenSlopeFld);
        
        a_CD_thresh(end + 1) = pdata.stage2.CD_TRAJ_LATENCY_THRESH;
        
        if ~isempty(pdata.stage2.(compenSlopeFld)) && ~isnan(pdata.stage2.(compenSlopeFld))
            SS1 = pdata.stage2.(SSreg_lineFld);
            SS2 = pdata.stage2.(SSreg_quadFld);
            df1 = pdata.stage2.(df_lineFld);
            df2 = pdata.stage2.(df_quadFld);
            
            F = (SS1 - SS2) / SS2 / ((df2 - df1) / df2);
            a_p_quadVline.(grp)(end + 1) = 1 - fcdf(abs(F), df1, df2);  
            
            a_R2_lineSpline.(grp)(end + 1) = pdata.stage2.(R2_lineFld);
        else
            a_p_quadVline.(grp)(end + 1) = NaN;
            
            a_R2_lineSpline.(grp)(end + 1) = NaN;
        end
        
        % Analyze the f1_trend
        meanF1Trend.(baseField) = nanmean(pdata.stage2.f1_trend.(baseField));
        stdF1Trend.(baseField) = nanstd(pdata.stage2.f1_trend.(baseField));
        f1Trend_UB = meanF1Trend.(baseField) + 1. * stdF1Trend.(baseField);
        f1Trend_LB = meanF1Trend.(baseField) - 1. * stdF1Trend.(baseField);
        
        propComp.down.(grp)(end + 1) = calc_prop(pdata.stage2.f1_trend.(downFld), f1Trend_UB, 'gt');
        propFoll.down.(grp)(end + 1) = calc_prop(pdata.stage2.f1_trend.(downFld), f1Trend_LB, 'lt');
        propComp.up.(grp)(end + 1) = calc_prop(pdata.stage2.f1_trend.(upFld), f1Trend_LB, 'lt');
        propFoll.up.(grp)(end + 1) = calc_prop(pdata.stage2.f1_trend.(upFld), f1Trend_UB, 'gt');
        
        % Noisiness in the [20, 40]-Hz frequency band
        bpRMS.downUp.(grp)(i1) = rms([pdata.stage2.bpRMS.([downFld, '_none']), ...
                                      pdata.stage2.bpRMS.([upFld, '_none'])]);
        
    end
    
    plot([taxis0(1),taxis0(end)],[0,0],'-','Color',[0.5,0.5,0.5]);
    
    if i0==2
        set(gca,'YTickLabel',{});
    else
        ylabel('Fraction F1 change from noPert');
    end
    
    avg_avgChg_F1_ML.(grp).down=mean(a_avgChg_traj_F1_ML.(grp).down');
    avg_avgChg_F1_ML.(grp).up=mean(a_avgChg_traj_F1_ML.(grp).up');
    sem_avgChg_F1_ML.(grp).down=std(a_avgChg_traj_F1_ML.(grp).down')/sqrt(size(a_avgChg_traj_F1_ML.(grp).down,2));
    sem_avgChg_F1_ML.(grp).up=std(a_avgChg_traj_F1_ML.(grp).up')/sqrt(size(a_avgChg_traj_F1_ML.(grp).up,2));
    
    if ~isempty(fsic(varargin, '1'))
        avg_avgChg_F1_ML.(grp).down_neg=mean(a_avgChg_traj_F1_ML.(grp).down_neg');
        avg_avgChg_F1_ML.(grp).up_neg=mean(a_avgChg_traj_F1_ML.(grp).up_neg');
        sem_avgChg_F1_ML.(grp).down_neg=std(a_avgChg_traj_F1_ML.(grp).down_neg')/sqrt(size(a_avgChg_traj_F1_ML.(grp).down_neg,2));
        sem_avgChg_F1_ML.(grp).up_neg=std(a_avgChg_traj_F1_ML.(grp).up_neg')/sqrt(size(a_avgChg_traj_F1_ML.(grp).up_neg,2));
        
        [foo, idx_0] = min(abs(taxis0 - 1e3 * spacingInt(1)));
        [foo, idx_1] = min(abs(taxis0 - 1e3 * spacingInt(2)));
        column_snum = (1 : numel(subjIDs.(grp)))';
        column_grp = i0 * ones(numel(subjIDs.(grp)), 1);
        column_down_large = mean(a_avgChg_traj_F1_ML.(grp).down(idx_0: idx_1, :))';
        column_up_large = mean(a_avgChg_traj_F1_ML.(grp).up(idx_0: idx_1, :))';
        column_down_small = mean(a_avgChg_traj_F1_ML.(grp).down_neg(idx_0: idx_1, :))';
        column_up_small = mean(a_avgChg_traj_F1_ML.(grp).up_neg(idx_0: idx_1, :))';
        spacing_chgF1_tab = [spacing_chgF1_tab; [column_snum, column_grp, column_down_large, column_down_small, column_up_large, column_up_small]];
    end
    
    % -- Construct the xls table for SYSTAT ANOVA analysis --
    
    
    avg_sd_F1_ML.(grp).none = mean(a_sd_F1_ML.(grp).none');
    avg_sd_F1_ML.(grp).down = mean(a_sd_F1_ML.(grp).down');
    avg_sd_F1_ML.(grp).up = mean(a_sd_F1_ML.(grp).up');
    sem_sd_F1_ML.(grp).none = std(a_sd_F1_ML.(grp).none') / sqrt(size(a_sd_F1_ML.(grp).none, 2));
    sem_sd_F1_ML.(grp).down = std(a_sd_F1_ML.(grp).down') / sqrt(size(a_sd_F1_ML.(grp).down, 2));
    sem_sd_F1_ML.(grp).up = std(a_sd_F1_ML.(grp).up')  / sqrt(size(a_sd_F1_ML.(grp).up, 2));
    
    avg_sdRatio_F1_ML.(grp).down = mean(a_sdRatio_F1_ML.(grp).down');
    avg_sdRatio_F1_ML.(grp).up = mean(a_sdRatio_F1_ML.(grp).up');
    sem_sdRatio_F1_ML.(grp).down = std(a_sdRatio_F1_ML.(grp).down') / sqrt(size(a_sdRatio_F1_ML.(grp).down, 2));
    sem_sdRatio_F1_ML.(grp).up = std(a_sdRatio_F1_ML.(grp).up') / sqrt(size(a_sdRatio_F1_ML.(grp).up, 2));
    
    avg_avgChg_F1_ML.(grp).downUp=mean(a_avgChg_traj_F1_ML.(grp).downUp');    
    sem_avgChg_F1_ML.(grp).downUp=std(a_avgChg_traj_F1_ML.(grp).downUp')/sqrt(size(a_avgChg_traj_F1_ML.(grp).downUp,2));
    
    if ~isempty(fsic(varargin, '1'))
        set(gca,'XLim',[taxis0(1),taxis0(end)],'YLim',[-0.175,0.175]);
    else
        set(gca,'XLim',[taxis0(1),taxis0(end)],'YLim',[-0.125,0.125]);
    end
    
    if i0 == 1
        xs = get(gca, 'XLim'); ys = get(gca, 'YLim');
        ezlegend([xs(1) + 0.05 * range(xs), ys(2) - 0.3 * range(ys), 0.5 * range(xs), 0.2 * range(ys)], ...
                  0.45, {'Down', 'Up'}, ...
                  {colors.down, colors.up}, ...
                  repmat(fontSize, 1, 4), {'-', '-', '-', '-'}, {colors.down, colors.up}, ...
                  [1, 1], [0, 0]);
    end
    
    xlabel('Time (ms)');
    title(grp);
end

%% Correlation between bpRMS and upDownContra (noisiness of the compensatory responses and the
%% magnitude of the compensatory responses)
figure;
for i0 = 1 : numel(groups)
    grp = groups{i0};
    plot(bpRMS.downUp.(grp), upDownContra.(grp), ...
         'o', 'Color', colors.(grp));
     
    [k, r2, p] = lincorr(bpRMS.downUp.(grp), upDownContra.(grp));
    
    fprintf('%s: r2 = %f, p = %f\n', grp, r2, p);
    hold on;
end

[h, p] = ttest2(bpRMS.downUp.PFS, bpRMS.downUp.PWS);
fprintf('Between-group t-test: p = %.4f\n', p);

[foo, idx] = sort(bpRMS.downUp.PWS, 'descend');
fprintf('Five PWS with the noisiest F1 compensation in the [20, 40]-Hz band: \n');
for i0 = 1 : 5
    fprintf('\t%s\n', subjIDs.PWS{idx(i0)});
end

%% Mixed ANOVA on level and duration
grps = fields(a_avgLv);
for i1 = 1 : numel(grps)
    grp = grps{i1};
    fprintf('Group: %s\n', grp);
    fprintf('\tLevels (mean+/-1SE): %.2f+/-%.2f s (none); %.2f+/-%.2f s (down); %.2f+/-%.2f s (up)\n', ...
            mean(a_avgLv.(grp)(:, 1)), ste(a_avgLv.(grp)(:, 1)), ...
            mean(a_avgLv.(grp)(:, 2)), ste(a_avgLv.(grp)(:, 2)), ...
            mean(a_avgLv.(grp)(:, 3)), ste(a_avgLv.(grp)(:, 3)));
    
    fprintf('\tDurations (mean+/-1SE): %.4f+/-%.4f s (none); %.4f+/-%.4f s (down); %.4f+/-%.4f s (up)\n', ...
            mean(a_avgDur.(grp)(:, 1)), ste(a_avgDur.(grp)(:, 1)), ...
            mean(a_avgDur.(grp)(:, 2)), ste(a_avgDur.(grp)(:, 2)), ...
            mean(a_avgDur.(grp)(:, 3)), ste(a_avgDur.(grp)(:, 3)));
end
    

do_RMAOV_1G1W(a_avgLv, 0.05);
do_RMAOV_1G1W(a_avgDur, 0.05);

%% Analyze the cross-trial adaptation effects in the noPert trials
figure('Position',[100,100,720,320]);
if isequal(baseField, 'none0') || isequal(baseField, 'none1')
    if isequal(baseField, 'none1')
        a_avgF1Traj_none0neg_ML = a_avgF1Traj_none1neg_ML;
        a_avgF1Traj_none0neg_ML.PFS.none0 = a_avgF1Traj_none0neg_ML.PFS.none1;
        a_avgF1Traj_none0neg_ML.PWS.none0 = a_avgF1Traj_none0neg_ML.PWS.none1;
    end
    
    groups = fields(subjIDs);
    for i0 = 1 : numel(groups)
        subplot('Position',[0.14+(i0-1)*0.42,0.12,0.42,0.8]);
        set(gca, 'FontSize', 11);
        grp = groups{i0};
        avg_chg.(grp).down = mean(a_avgF1Traj_none0neg_ML.(grp).down - a_avgF1Traj_none0neg_ML.(grp).none0, 2);
        avg_chg.(grp).up = mean(a_avgF1Traj_none0neg_ML.(grp).up - a_avgF1Traj_none0neg_ML.(grp).none0, 2);
        ste_chg.(grp).down = std(a_avgF1Traj_none0neg_ML.(grp).down - a_avgF1Traj_none0neg_ML.(grp).none0, [], 2) / ...
                          sqrt(size(a_avgF1Traj_none0neg_ML.(grp).down, 2));
        ste_chg.(grp).up = std(a_avgF1Traj_none0neg_ML.(grp).up - a_avgF1Traj_none0neg_ML.(grp).none0, [], 2) / ...
                          sqrt(size(a_avgF1Traj_none0neg_ML.(grp).up, 2));
                      
        taxis0 = 1e3 * (0:FRAME_DUR:FRAME_DUR*(size(a_avgF1Traj_none0neg_ML.(grp).none0, 1)-1));

        
        hold on;
        sdirs = {'down', 'up'};
        for i1 = 1 : numel(sdirs)
            sdir = sdirs{i1};
            plot(taxis0, avg_chg.(grp).(sdir), '-', 'Color', colors.(sdir), 'LineWidth', 2);
            plot(taxis0, avg_chg.(grp).(sdir) - ste_chg.(grp).(sdir), '--', 'Color', colors.(sdir));
            plot(taxis0, avg_chg.(grp).(sdir) + ste_chg.(grp).(sdir), '--', 'Color', colors.(sdir));
        end
        plot([taxis0(1), taxis0(end)], [0, 0], '-', 'Color', [0.5, 0.5, 0.5]);
        
        set(gca, 'XLim', [taxis0(1), taxis0(end)]);
        set(gca, 'YLim', [-12, 12]);
        
        if i0 == 1
            title('PFS');
            ylabel('F1 change from clean noPert (Hz) (mean\pm1 SEM)');
            xs = get(gca, 'XLim'); ys = get(gca, 'YLim');
            
        else
            title('PWS');
            set(gca, 'YTickLabel', []);
            ezlegend([xs(1) + 0.025 * range(xs), ys(1) + 0.03 * range(ys), 0.95 * range(xs), 0.15 * range(ys)], ...
                  0.45, {'contaminated by Down', 'contaminated by Up'}, ...
                  {colors.down, colors.up}, ...
                  repmat(fontSize - 2, 1, 4), {'-', '-', '-', '-'}, {colors.down, colors.up}, ...
                  [2, 2], [0, 0]);
        end
        
        xlabel('Time (ms)');
        
    end
end

%% Proportion of compensating and following responses
for i0 = 1 : numel(groups)
    grp = groups{i0};
    propComp.both.(grp) = mean([propComp.down.(grp); propComp.up.(grp)]);
    propFoll.both.(grp) = mean([propFoll.down.(grp); propFoll.up.(grp)]);
    propNull.both.(grp) = ones(size(propComp.both.(grp))) - propComp.both.(grp) - propFoll.both.(grp);
    
    fprintf('%s: propComp (both) = %.4f +/- %.4f; propNull (both) = %.4f +/- %.4f; propFoll = %.4f +/- %.4f\n', ...
            grp, ...
            mean(propComp.both.(grp)), ste(propComp.both.(grp)), ...
            mean(propNull.both.(grp)), ste(propNull.both.(grp)), ...
            mean(propFoll.both.(grp)), ste(propFoll.both.(grp)));
end

[h, p]  = ttest2(propComp.both.PWS, propComp.both.PFS);

%% Spacing F1chg figure
if ~isempty(fsic(varargin, '1'))
    chgF1_ssDown_PFS = spacing_chgF1_tab(spacing_chgF1_tab(:, 2) == 1, 4);
    chgF1_ssDown_PWS = spacing_chgF1_tab(spacing_chgF1_tab(:, 2) == 2, 4);
    chgF1_ssUp_PFS = spacing_chgF1_tab(spacing_chgF1_tab(:, 2) == 1, 6);
    chgF1_ssUp_PWS = spacing_chgF1_tab(spacing_chgF1_tab(:, 2) == 2, 6);

    figure('Position', [50, 50, 400, 300]);
    set(gca, 'FontSize', 12);
    hold on;
    plot([0.75, 2.25], [0, 0], '-', 'Color', [0.5, 0.5, 0.5]);
    errorbar([1, 2], [mean(chgF1_ssDown_PFS), mean(chgF1_ssUp_PFS)], ...
                     [ste(chgF1_ssDown_PFS), ste(chgF1_ssUp_PFS)], 'Color', 'k');
    plot([1, 2], [mean(chgF1_ssDown_PFS), mean(chgF1_ssUp_PFS)], 'o', 'Color', [0, 0, 0]);
    errorbar([1, 2], [mean(chgF1_ssDown_PWS), mean(chgF1_ssUp_PWS)], ...
                     [ste(chgF1_ssDown_PWS), ste(chgF1_ssUp_PWS)], 'Color', [160, 32, 240] / 255);
    plot([1, 2], [mean(chgF1_ssDown_PWS), mean(chgF1_ssUp_PWS)], 's', 'Color', [160, 32, 240] / 255);
    set(gca, 'XLim', [0.75, 2.25], 'XTick', [1, 2], 'XTickLabel', {'Down', 'Up'});
    ylabel('Fraction F1 change from noPert in small-spacing trials (Mean\pm1 SEM)');
    xs = get(gca, 'XLim'); ys = get(gca, 'YLim');
    ezlegend([xs(1) + 0.05 * range(xs), ys(2) - 0.3 * range(ys), 0.45 * range(xs), 0.2 * range(ys)], ...
                      0.5, {'PFS', 'PWS'}, ...
                      {[0, 0, 0], [160, 32, 240] / 255}, ...
                      repmat(fontSize, 1, 4), {'o-', 's-'}, {[0, 0, 0], [160, 32, 240] / 255}, ...
                      [1, 1], [0, 0]);
end

%% Write to xls file for spacing analysis
if ~isempty(fsic(varargin, '1'))
    spacing_xls_fn = fullfile(cdds, 'stut', 'xlsSysEH', 'spacing_chgF1.xls');
    if isfile(spacing_xls_fn)
        delete(spacing_xls_fn);
    end

    spacing_chgF1_cell = [spacing_chgF1_row1; num2cell(spacing_chgF1_tab)];
    xlswrite(spacing_xls_fn, spacing_chgF1_cell, 1);
    fprintf('Wrote spacing chgF1 data to xls file: %s\n', spacing_xls_fn);
end


%% Extract the data for mixed ANOVA use
flds = {'down', 'up'};
maov_tab = struct;
maov_tab_downUp = struct;
ntp = numel(mAOV_tPoints);
rmaov32_tab = nan(0, 5);

if taxis0(end) > 10
    taxis0 = taxis0 / 1e3;
end

% 1st row for the xls file:
row1 = {'SNUM', 'GRP'};
for i1 = 1 : numel(flds)
    for i2 = 1 : ntp
        row1{end + 1} = sprintf('S%dT%d', i1, i2);
    end
end
xls_tab = [];

% 1st row for the xls_downUp file:
row1_downUp = {'SNUM', 'GRP'};
for i1 = 1 : ntp
    row1_downUp{end + 1} = sprintf('T%d', i1); 
end
xls_tab_downUp = [];

for i0 = 1 : numel(groups)
    grp = groups{i0};
    maov_tab.(grp) = struct;
    
    maov_tab.(grp).down = nan(numel(subjIDs.(grp)), ntp);
    maov_tab.(grp).up = nan(numel(subjIDs.(grp)), ntp);
        
    for i2 = 1 : numel(subjIDs.(grp))
        xls_row = [i2, i0];
        xls_row_downUp = [i2, i0];
        
        for i1 = 1 : numel(flds)
            for i3 = 1 : ntp
                fld = flds{i1};
                [foo, idx0] = min(abs(taxis0 - mAOV_tPoints{i3}(1)));
                [foo, idx1] = min(abs(taxis0 - mAOV_tPoints{i3}(2)));
                
                if isempty(fsic(varargin, 'noSmoothAOV'))
                    avg_chg = a_avgChg_traj_F1_ML.(grp).(fld);
                else
                    avg_chg = a_avgChg_traj_F1_ns_ML.(grp).(fld);
                end
                
                if ~isempty(fsic(varargin, 'tPointEnd'))
                    maov_tab.(grp).(fld)(i2, i3) = avg_chg(idx1, i2);
                    rmaov32_tab = [rmaov32_tab; [avg_chg(idx1, i2), i0, i1, i3, i2]];
                    xls_row(end + 1) = avg_chg(idx1, i2);
                else
                    maov_tab.(grp).(fld)(i2, i3) = nanmean(avg_chg(idx0 : idx1, i2));
                    rmaov32_tab = [rmaov32_tab; [nanmean(avg_chg(idx0 : idx1, i2)), i0, i1, i3, i2]];
                    xls_row(end + 1) = nanmean(avg_chg(idx0 : idx1, i2));
                end
            end
        end
        xls_row_downUp = [xls_row_downUp, xls_row(3 : 3 + ntp - 1) - xls_row(3 + ntp : end)];
        
        xls_tab = [xls_tab; xls_row];
        xls_tab_downUp = [xls_tab_downUp; xls_row_downUp];
    end
end

% Apply RMAOV32 
%    (2 repeated measures: fld and time)
%    (1 between-subjects: grp)
% rmaov32_res = RMAOV32(rmaov32_tab, 0.05); % Caveat: requires equal number of subjects in both groups. Have to use SYSTAT for unbalanced group sizes.

%
xls_cell = [row1; num2cell(xls_tab)];
if isfile(mAOV_xlsFN)
    delete(mAOV_xlsFN);
end
xlswrite(mAOV_xlsFN, xls_cell, 1);
fprintf('Wrote mixed ANOVA (1G2W) data table to %s\n', mAOV_xlsFN);

xls_cell_downUp = [row1_downUp; num2cell(xls_tab_downUp)];
if isfile(mAOV_downUp_xlsFN)
    delete(mAOV_downUp_xlsFN);
end
xlswrite(mAOV_downUp_xlsFN, xls_cell_downUp, 1);
fprintf('Wrote mixed ANOVA (1G1W) down-up data table to %s\n', mAOV_downUp_xlsFN);

%% Coarse time-scale composite response curves
figure('Position', [100, 100, 480, 320]);
set(gca, 'FontSize', fontSize);
hold on;
plot([0, numel(mAOV_tPoints) + 1], [0, 0], '-', 'Color', [0.5, 0.5, 0.5], ...
     'LineWidth', 1);
errorbar(1 : numel(mAOV_tPoints), mean(maov_tab.PFS.down - maov_tab.PFS.up), ...
         ste(maov_tab.PFS.down - maov_tab.PFS.up), ...
         'Color', colors.PFS, 'LineWidth', 1.5);
errorbar(1 : numel(mAOV_tPoints), mean(maov_tab.PWS.down - maov_tab.PWS.up), ...
         ste(maov_tab.PWS.down - maov_tab.PWS.up), ...
         'Color', colors.PWS, 'LineWidth', 1.5);
set(gca, 'XLim', [0, numel(mAOV_tPoints) + 1]);
set(gca, 'XTick', 1 : 12);
xs = get(gca, 'XLim'); ys = get(gca, 'YLim');

for i1 = 1 : numel(mAOV_tPoints)
    [h, p] = ttest2(maov_tab.PFS.down(:, i1) - maov_tab.PFS.up(:, i1), ...
                    maov_tab.PWS.down(:, i1) - maov_tab.PWS.up(:, i1));
    if p < 0.05
        plot(i1, ys(2) - 0.05 * range(ys), 'k*', 'MarkerSize', 11);
    end
end

if ~isempty(fsic(varargin, 'tPointEnd'))
    xlabel('Time bin #');
else
    xlabel('Time point #');
end
ylabel('F1 difference between Down and Up (fract.)');

ezlegend([xs(1) + 0.65 * range(xs), ys(1) + 0.03 * range(ys), 0.33 * range(xs), 0.14 * range(ys)], ...
          0.5, {'PFS', 'PWS'}, ...
          {[0, 0, 0], [160, 32, 240] / 255}, ...
          repmat(fontSize, 1, 4), {'-', '-'}, {[0, 0, 0], [160, 32, 240] / 255}, ...
          [1.5, 1.5], [0, 0]);
     
%% Print stats about discard:
for i1 = 1 : numel(groups)
    grp = groups{i1};
    fprintf('%s: %d of %d (%.4f%%) trials discarded,\n', grp, ...
            sum(a_nDscd_all.(grp)), sum(a_nTotTrials.(grp)), ...
            1e2 * sum(a_nDscd_all.(grp)) / sum(a_nTotTrials.(grp)));
    fprintf('\t%d (%.4f%%) because of dysfluencies or speech error.\n', ...
            sum(a_nDscd_prodErr.(grp)), ...
            1e2 * sum(a_nDscd_prodErr.(grp)) / sum(a_nTotTrials.(grp)));
end

%% Print the ratio of discarded trials due to minimum length
fprintf('\n');
for i1 = 1 : numel(groups)
    grp = groups{i1};
    fprintf('minLen discarded: %s: none = %.2f +/- %.2f%% (SD); down = %.2f +/- %.2f%% (SD); up = %.2f +/- %.2f%% (SD)\n', ...
        grp, mean(ratioLenDiscard.(grp).none) * 1e2, std(ratioLenDiscard.(grp).none) * 1e2, ...
             mean(ratioLenDiscard.(grp).down) * 1e2, std(ratioLenDiscard.(grp).down) * 1e2, ...
             mean(ratioLenDiscard.(grp).up) * 1e2, std(ratioLenDiscard.(grp).up) * 1e2);
    totLenDiscardRatio.(grp) = (ratioLenDiscard.(grp).none * 80 + ratioLenDiscard.(grp).down * 20 + ratioLenDiscard.(grp).up * 20) / 120;
    fprintf('\tTotal = %.2f +/- %.2f%% (SD)\n', mean(totLenDiscardRatio.(grp)) * 1e2, std(totLenDiscardRatio.(grp)) * 1e2);
end
fprintf('\n');

[h, p, ci, t_stats] = ttest2(totLenDiscardRatio.PWS, totLenDiscardRatio.PFS);
fprintf('\nBetween-group comparison of totLenDiscardRatio:\n\tt(%d) = %f, p = %f (t-test)\n', t_stats.df, t_stats.tstat, p);

[p_rs, h_rs] = ranksum(totLenDiscardRatio.PWS, totLenDiscardRatio.PFS);
fprintf('\nBetween-group comparison of totLenDiscardRatio:\n\tp = %f (ranksum)\n', p_rs);

%% Print the ratio of kept trials under downFld and upFld (lag criterion)
for i1 = 1 : numel(groups)
    grp = groups{i1};
    fprintf('Lag criterion: %s: downFld = %.2f +/- %.2f%% (SD); upFld = %.2f +/- %.2f%% (SD)\n', ...
        grp, mean(ratioKeptTrials.(grp).down) * 1e2, std(ratioKeptTrials.(grp).down) * 1e2, ...
             mean(ratioKeptTrials.(grp).up) * 1e2, std(ratioKeptTrials.(grp).up) * 1e2);
end


%% Stats: P-values for within and between group comparisons
p_t_onePert = struct;
p_t_twoPert = struct;
FDR_p_thresh_onePert = struct;
FDR_P_thresh_twoPert = struct;
for i1 = 1 : numel(groups)
    grp = groups{i1};
    p_t_onePert.(grp) = struct;
    FDR_p_thresh_onePert.(grp) = struct;
    
    perts = {'down', 'up'};
    
    for i2 = 1 : numel(perts)
        pert = perts{i2};
        p_t_onePert.(grp).(pert) = nan(N_SAMP, 1);        
        
        mat = a_avgChg_traj_F1_ML.(grp).(pert);
        for i3 = 1 : N_SAMP
            [h, p_t_onePert.(grp).(pert)(i3)] = ttest(a_avgChg_traj_F1_ML.(grp).(pert)(i3, :));
        end
        FDR_p_thresh_onePert.(grp).(pert) = fdr(p_t_onePert.(grp).(pert), FDR);
    end
    
    p_t_twoPert.(grp) = nan(N_SAMP, 1);
    for i3 = 1 : N_SAMP
        [h, p_t_twoPert.(grp)(i3)] = ttest(a_avgChg_traj_F1_ML.(grp).(perts{1})(i3, :), a_avgChg_traj_F1_ML.(grp).(perts{2})(i3, :));
    end
    
    FDR_p_thresh_twoPert.(grp) = fdr(p_t_twoPert.(grp), FDR);
end

p_t_twoGrp = nan(N_SAMP, 1);
p_onePert_twoGrp.down = nan(N_SAMP, 1);
p_onePert_twoGrp.up = nan(N_SAMP, 1);

for i1 = 1 : N_SAMP
    [h, p_t_twoGrp(i1)] = ttest2(a_avgChg_traj_F1_ML.PWS.downUp(i1, :), a_avgChg_traj_F1_ML.PFS.downUp(i1, :));
    
%     [h, p_onePert_twoGrp.down(i1)] = ttest2(a_avgChg_traj_F1_ML.PWS.down(i1, :), a_avgChg_traj_F1_ML.PFS.down(i1, :));
%     [h, p_onePert_twoGrp.up(i1)] = ttest2(a_avgChg_traj_F1_ML.PWS.up(i1, :), a_avgChg_traj_F1_ML.PFS.up(i1, :));
%     [p_t_twoGrp(i1), h] = ranksum(a_avgChg_traj_F1_ML.PWS.downUp(i1, :), a_avgChg_traj_F1_ML.PFS.downUp(i1, :));

    [p_onePert_twoGrp.down(i1), h] = ranksum(a_avgChg_traj_F1_ML.PWS.down(i1, :), a_avgChg_traj_F1_ML.PFS.down(i1, :));
    [p_onePert_twoGrp.up(i1), h] = ranksum(a_avgChg_traj_F1_ML.PWS.up(i1, :), a_avgChg_traj_F1_ML.PFS.up(i1, :));
end
FDR_p_thresh_twoGrp = fdr(p_t_twoGrp, FDR);

%% Latency comparison
CD_thresh = unique(a_CD_thresh);
if length(CD_thresh) ~= 1
    error('CD threshold unequal among the subects.');
end

figure('Position', [100, 100, 380, 300]);
set(gca, 'FontSize', fontSize);
boxplot(1e3 * [a_lat.PFS(~isnan(a_lat.PFS)), a_lat.PWS(~isnan(a_lat.PWS))], ...
    [1 * ones(size(a_lat.PFS(~isnan(a_lat.PFS)))), 2 * ones(size(a_lat.PWS(~isnan(a_lat.PWS))))]);
hold on;
for i1 = 1 : numel(groups)
    grp = groups{i1};
    
    t_lat = a_lat.(grp); 
    t_lat = t_lat(find(~isnan(t_lat)));    
    plot(i1 + 0.25, t_lat * 1e3, 'o', 'Color', colors.(grp));
    hold on;
    
    plot(i1 + 0.35, 1e3 * nanmean(a_lat.(grp)), 'o', 'Color', colors.(grp));
    plot(repmat(i1 + 0.35, 1, 2), 1e3 * (nanmean(a_lat.(grp)) + [-1, 1] * nanste(a_lat.(grp))), ...
         '-', 'Color', colors.(grp));
end

set(gca, 'XLim', [0, 3], 'XTick', [1, 2], 'XTickLabel', groups);



[h_t, p_t, ci, stats] = ttest2(a_lat.PFS, a_lat.PWS);
[p_rs, h_rs] = ranksum(a_lat.PFS(~isnan(a_lat.PFS)), a_lat.PWS(~isnan(a_lat.PWS)));
cohend = abs(cohen_d(a_lat.PFS(~isnan(a_lat.PFS)), a_lat.PWS(~isnan(a_lat.PWS))));
xs = get(gca, 'XLim'); ys = get(gca, 'YLim');
text(xs(1) + 0.5 * range(xs), ys(2) - 0.055 * range(ys), sprintf('t-test: t=%.4f; p=%.4f', stats.tstat, p_t));
text(xs(1) + 0.5 * range(xs), ys(2) - 0.11 * range(ys), sprintf('Cohen''s d=%.4f', cohend));
text(xs(1) + 0.5 * range(xs), ys(2) - 0.165 * range(ys), sprintf('Ranksum: p=%.4f', p_rs));
ylabel('Response latency (ms)');

figure('Position', [100, 50, 300, 300]);
set(gca, 'FontSize', fontSize + 2);
for i1 = 1 : 2
    if i1 == 1; grp = 'PFS'; 
    else; grp = 'PWS';
    end
    bar(i1, 1e3 * nanmean(a_lat.(grp)), 'EdgeColor', colors.(grp), 'FaceColor', 'none', 'LineWidth', 2);
    hold on;
    plot([i1, i1], 1e3 * nanmean(a_lat.(grp)) + [-1, 1] * 1e3 * nanste(a_lat.(grp)), '-', 'Color', colors.(grp), 'LineWidth', 2);
end
set(gca, 'XTick', [1, 2], 'XTickLabel', {'PFS', 'PWS'});
ylabel('Response latency (mean \pm 1 SEM)');
xs = get(gca, 'XLim'); ys = get(gca, 'YLim');
text(xs(1) + 0.05 * range(xs), ys(2) - 0.055 * range(ys), sprintf('t-test: p=%.4f', p_t), 'FontSize', fontSize + 1);
box off;

%% cdtraj fitting performances


%% Compensation slope comparison
for i0 = 1 : 2
    if i0 == 2
        a_compenSlope.PFS(isnan(a_compenSlope.PFS)) = 0;
        a_compenSlope.PWS(isnan(a_compenSlope.PWS)) = 0;
    end
    
    [h_t_compenSlope, p_t_compenSlope, ci, stats] = ttest2(a_compenSlope.PFS, a_compenSlope.PWS);
    [p_rs_compenSlope, h_rs_compenSlope] = ranksum(a_lat.PFS(~isnan(a_compenSlope.PFS)), a_lat.PWS(~isnan(a_compenSlope.PWS)));

    if i0 == 1
        figure('Position', [100, 50, 350, 350], 'Name', 'Fit slope comparison');
    else
        figure('Position', [100, 50, 350, 350], 'Name', 'Fit slope comparison (nans zeroed)');
    end
    
    set(gca, 'FontSize', fontSize + 2);
    for i1 = 1 : 2
        if i1 == 1; grp = 'PFS'; 
        else; grp = 'PWS';
        end
        bar(i1, 1e3 * nanmean(a_compenSlope.(grp)), 'EdgeColor', colors.(grp), 'FaceColor', 'none', 'LineWidth', 2);
        hold on;
        plot([i1, i1], 1e3 * nanmean(a_compenSlope.(grp)) + [-1, 1] * 1e3 * nanste(a_compenSlope.(grp)), '-', 'Color', colors.(grp), 'LineWidth', 2);
    end
    set(gca, 'XTick', [1, 2], 'XTickLabel', {'PFS', 'PWS'});
    ylabel('Response slope (mean \pm 1 SEM)');
    xs = get(gca, 'XLim'); ys = get(gca, 'YLim');
    text(xs(1) + 0.05 * range(xs), ys(2) - 0.055 * range(ys), sprintf('t-test: p=%.4f', p_t_compenSlope), 'FontSize', fontSize + 1);
    box off;
end
%% Print the number of subjects without latency: 
fprintf('\n');
for i1 = 1 : numel(groups)
    grp = groups{i1};
    n_noLat = numel(find(isnan(a_lat.(grp))));
    fprintf('%s: %d of %d subjects are without latency under CD threshold %f\n', ...
            grp, n_noLat, numel(a_lat.(grp)), CD_thresh);
end
fprintf('\n');

%% Show the absolute magnitude of the down and up perturbations. 
figure;
boxplot([a_F1_downPert.PFS, a_F1_downPert.PWS], [1 * ones(size(a_F1_downPert.PFS)), 2 * ones(size(a_F1_downPert.PWS))]);
set(gca, 'XLim', [0, 3], 'XTick', [1, 2], 'XTickLabel', {'PFS', 'PWS'});
ylabel('Absolute magnitude of down perturbation (Hz)');
[h_t, p_t] = ttest2(a_F1_downPert.PFS, a_F1_downPert.PWS);
cohend = cohen_d(a_F1_downPert.PFS, a_F1_downPert.PWS);
xs = get(gca, 'XLim'); ys = get(gca, 'YLim');
text(xs(1) + 0.05 * range(xs), ys(2) - 0.05 * range(ys),  sprintf('t-test: p = %.4f', p_t));
text(xs(1) + 0.05 * range(xs), ys(2) - 0.10 * range(ys),  sprintf('t-test: Cohen''s d = %.4f', cohend));

fprintf('Down pert: PFS: %.2f +/- %.2f; PWS: %.2f +/- %.2f\n', ...
    mean(a_F1_downPert.PFS), ste(a_F1_downPert.PFS), mean(a_F1_downPert.PWS), ste(a_F1_downPert.PWS));

figure;
boxplot([a_F1_upPert.PFS, a_F1_upPert.PWS], [1 * ones(size(a_F1_upPert.PFS)), 2 * ones(size(a_F1_upPert.PWS))]);
set(gca, 'XLim', [0, 3], 'XTick', [1, 2], 'XTickLabel', {'PFS', 'PWS'});
ylabel('Absolute magnitude of up perturbation (Hz)');
[h_t, p_t] = ttest2(a_F1_upPert.PFS, a_F1_upPert.PWS);
cohend = cohen_d(a_F1_upPert.PFS, a_F1_upPert.PWS);
xs = get(gca, 'XLim'); ys = get(gca, 'YLim');
text(xs(1) + 0.05 * range(xs), ys(2) - 0.05 * range(ys),  sprintf('t-test: p = %.4f', p_t));
text(xs(1) + 0.05 * range(xs), ys(2) - 0.10 * range(ys),  sprintf('t-test: Cohen''s d = %.4f', cohend));

fprintf('Up pert: PFS: %.2f +/- %.2f; PWS: %.2f +/- %.2f\n', ...
    mean(a_F1_upPert.PFS), ste(a_F1_upPert.PFS), mean(a_F1_upPert.PWS), ste(a_F1_upPert.PWS));

%% Adaptation plot
figure('Position', [100,100,720,320], 'Name', 'F1 Adaptation Effects');
set(gca, 'FontSize', fontSize);
for i1 = 1 : numel(groups)
    grp = groups{i1};
    subplot('Position',[0.125+(i1-1)*0.42,0.12,0.42,0.825]);
    
%     for i2 = 1 : numel(adapt_avgChg_F1_ML.(grp))        
    taxis0=0:FRAME_DUR:FRAME_DUR*(size(adapt_avgChg_F1_ML.(grp).aftDown,1)-1);
    plot(taxis0, adapt_avgChg_F1_ML.(grp).aftDown, 'Color', colors.down);
    hold on;
    plot(taxis0, adapt_avgChg_F1_ML.(grp).aftUp, 'Color', colors.up);        
%     end
    plot([taxis0(1),taxis0(end)],[0,0],'-','Color',[0.5,0.5,0.5]);
    
    avg_adapt_avgChg_F1_ML.(grp).aftDown=mean(adapt_avgChg_F1_ML.(grp).aftDown');
    avg_adapt_avgChg_F1_ML.(grp).aftUp=mean(adapt_avgChg_F1_ML.(grp).aftUp');
    sem_adapt_avgChg_F1_ML.(grp).aftDown=std(adapt_avgChg_F1_ML.(grp).aftDown')/sqrt(size(adapt_avgChg_F1_ML.(grp).aftDown,2));
    sem_adapt_avgChg_F1_ML.(grp).aftUp=std(adapt_avgChg_F1_ML.(grp).aftUp')/sqrt(size(adapt_avgChg_F1_ML.(grp).aftUp,2));

    title(grp);
    
    set(gca,'XLim',[taxis0(1),taxis0(end)],'YLim',[-0.14,0.14]);
end

%% Subset mode ('1'): showing the 1 and 1neg comparison
if ~isempty(fsic(varargin, '1'))
    p_t_1neg = struct;    
    FDR_p_thresh_1neg = struct;    
    for i1 = 1 : numel(groups)
        grp = groups{i1};
        p_t_1neg.(grp) = struct;
        FDR_p_thresh_1neg.(grp) = struct;

        perts = {'down', 'up'};

        for i2 = 1 : numel(perts)
            pert = perts{i2};
            pert_neg = [pert, '_neg'];
            p_t_1neg.(grp).(pert) = nan(N_SAMP, 1);        
            
            for i3 = 1 : N_SAMP
                [h, p_t_1neg.(grp).(pert)(i3)] = ttest(a_avgChg_traj_F1_ML.(grp).(pert)(i3, :), a_avgChg_traj_F1_ML.(grp).(pert_neg)(i3, :));
            end
            FDR_p_thresh_1neg.(grp).(pert) = fdr(p_t_onePert.(grp).(pert), FDR);
        end
        
        
    end        

    figure('Position', [100,100,720,700], 'Name', 'Subset mode: ''1'' vs. ''1neg''');
    
    for i1 = 1 : numel(groups)
        grp = groups{i1};
    
        taxis0=0:FRAME_DUR:FRAME_DUR*(length(avg_avgChg_F1_ML.(grp).down)-1);
        for i2 = 1 : 2
            subplot('Position',[0.125+(i1-1)*0.42,0.06 + (i2-1)*.44,0.42,0.44]);
            set(gca, 'FontSize', fontSize);
            
            if i2 == 1 fld = 'down'; else fld = 'up'; end
            fld_neg = [fld, '_neg'];
            plot(taxis0, avg_avgChg_F1_ML.(grp).(fld), '-', 'LineWidth', 1.5, 'Color', colors.(fld));
            hold on;
            plot(taxis0, avg_avgChg_F1_ML.(grp).(fld) - sem_avgChg_F1_ML.(grp).(fld), '--', 'Color', colors.(fld));
            plot(taxis0, avg_avgChg_F1_ML.(grp).(fld) + sem_avgChg_F1_ML.(grp).(fld), '--', 'Color', colors.(fld));
            plot(taxis0, avg_avgChg_F1_ML.(grp).(fld_neg), 'Color', colors.(fld_neg));
            plot(taxis0, avg_avgChg_F1_ML.(grp).(fld_neg) - sem_avgChg_F1_ML.(grp).(fld), '--', 'Color', colors.(fld_neg));
            plot(taxis0, avg_avgChg_F1_ML.(grp).(fld_neg) + sem_avgChg_F1_ML.(grp).(fld), '--', 'Color', colors.(fld_neg));
            
            plot([taxis0(1),taxis0(end)],[0,0],'-','Color',[0.5,0.5,0.5]);
            
            set(gca,'XLim',[taxis0(1),taxis0(end)],'YLim',[-0.1,0.1]);
            if i1 == 2
                set(gca, 'YTickLabel', []);
            else
                ylabel('F1 change from noPert baseline (fraction)');
            end
            if i2 == 2
                set(gca, 'XTickLabel', []);   
            else
                xlabel('Time (ms)');
            end
            xs = get(gca, 'XLim'); ys = get(gca, 'YLim');
            
            draw_sgn_bar(taxis0, p_t_1neg.(grp).(fld), 0.05, FDR_p_thresh_1neg.(grp).(pert), ...
                 ys(1) + 0.20 * range(ys), 0.02 * range(ys), 'k', [0.75, 0.75, 0.75], [0, 0, 0]);
            text(xs(1) + 0.05 * range(xs), ys(2) - 0.06 * range(ys), [grp, ': ', fld], 'Color', 'k', 'FontSize', fontSize);
           
            if i1 == 1 && i2 == 2
                ezlegend([xs(1) + 0.05 * range(xs), ys(2) - 0.4 * range(ys), 0.77 * range(xs), 0.3 * range(ys)], ...
                     0.45, {'Up: large spacing', 'Up: small spacing', 'Down: large spacing', 'Down: small spacing'}, ...
                     {colors.up, colors.up_neg, colors.down, colors.down_neg}, ...
                     repmat(fontSize, 1, 4), {'-', '-', '-', '-'}, {colors.up, colors.up_neg, colors.down, colors.down_neg}, ...
                     [1, 1, 1, 1], [1, 1, 1, 1]);
            end
        end
        
       
%         title(grp);

        
    end
    
    % Write to xls file for ANOVA (1G2W in SYSTAT)
    % 1G: PWS, 
    % 2W: 1) SDIR (), 2) SPACING ({SMALL, LARGE})
   
    
end

%%
figure('Position',[100,100,800,400]);
for i0=1:numel(groups)
    grp=groups{i0};
    subplot('Position',[0.125+(i0-1)*0.42,0.12,0.42,0.8]);
    set(gca, 'FontSize', fontSize);
    
    taxis0 = 1e3 * (0:FRAME_DUR:FRAME_DUR*(size(avg_avgChg_F1_ML.(grp).down,2)-1));
    plot(taxis0,avg_avgChg_F1_ML.(grp).down,'-','Color',colors.down,'LineWidth',2); hold on;
    plot(taxis0,avg_avgChg_F1_ML.(grp).down-sem_avgChg_F1_ML.(grp).down,'--','Color',colors.down);    
    plot(taxis0,avg_avgChg_F1_ML.(grp).down+sem_avgChg_F1_ML.(grp).down,'--','Color',colors.down);    
    plot(taxis0,avg_avgChg_F1_ML.(grp).up,'-','Color',colors.up);
    plot(taxis0,avg_avgChg_F1_ML.(grp).up-sem_avgChg_F1_ML.(grp).up,'--','Color',colors.up);
    plot(taxis0,avg_avgChg_F1_ML.(grp).up+sem_avgChg_F1_ML.(grp).up,'--','Color',colors.up);
    
    plot([taxis0(1),taxis0(end)],[0,0],'-','color',[0.5,0.5,0.5]);
    
    if i0==2
        set(gca,'YTickLabel',{});
    else
        ylabel('Fraction F1 change from noPert');
    end
    
    set(gca,'XLim',[taxis0(1),taxis0(end)],'YLim',[-0.06,0.06]);
    xlabel('Time (ms)');
    title(grp);
    
    % -- Draw the significance bars --
    xs = get(gca, 'XLim'); ys = get(gca, 'YLim');
    draw_sgn_bar(taxis0, p_t_onePert.(grp).down, 0.05, FDR_p_thresh_onePert.(grp).down, ...
                 ys(1) + 0.20 * range(ys), 0.02 * range(ys), 'k', 0.75 + 0.25 * colors.down, colors.down);
    draw_sgn_bar(taxis0, p_t_onePert.(grp).up, 0.05, FDR_p_thresh_onePert.(grp).up, ...
                 ys(1) + 0.12 * range(ys), 0.02 * range(ys), 'k', 0.75 + 0.25 * colors.up, colors.up);
    draw_sgn_bar(taxis0, p_t_twoPert.(grp), 0.05, FDR_p_thresh_twoPert.(grp), ...
                 ys(1) + 0.04 * range(ys), 0.02 * range(ys), 'k', [0.75, 0.75, 0.75], 'k');
             
    if i0 == 1        
        ezlegend([xs(1) + 0.05 * range(xs), ys(2) - 0.25 * range(ys), 0.76 * range(xs), 0.225 * range(ys)], ...
                 0.45, {'Down (mean\pm1 SEM)', 'Up (mean\pm 1 SEM)'}, {colors.down, colors.up}, ...
                 [fontSize - 1, fontSize - 1], {'-', '-'}, {colors.down, colors.up}, ...
                 [1, 1], [1, 1]);
             
        rectangle('Position', [xs(1) + 0.05 * range(xs), ys(1) + 0.32 * range(ys), 0.30 * range(xs), 0.02 * range(ys)], ...
                  'EdgeColor', 'k', 'FaceColor', 'w');
        text(xs(1) + 0.36 * range(xs), ys(1) + 0.335 * range(ys), 'n.s.');
        
        rectangle('Position', [xs(1) + 0.05 * range(xs), ys(1) + 0.29 * range(ys), 0.10 * range(xs), 0.02 * range(ys)], ...
                  'EdgeColor', 'k', 'FaceColor', 0.5 + 0.5 * colors.down);
        rectangle('Position', [xs(1) + 0.15 * range(xs), ys(1) + 0.29 * range(ys), 0.10 * range(xs), 0.02 * range(ys)], ...
                  'EdgeColor', 'k', 'FaceColor', 0.5 + 0.5 * colors.up);
        rectangle('Position', [xs(1) + 0.25 * range(xs), ys(1) + 0.29 * range(ys), 0.10 * range(xs), 0.02 * range(ys)], ...
                  'EdgeColor', 'k', 'FaceColor', 0.5 + 0.5 * [0, 0, 0]);
        text(xs(1) + 0.36 * range(xs), ys(1) + 0.305 * range(ys), 'p<0.05 (unc.)');
              
        rectangle('Position', [xs(1) + 0.05 * range(xs), ys(1) + 0.26 * range(ys), 0.10 * range(xs), 0.02 * range(ys)], ...
                  'EdgeColor', 'k', 'FaceColor', colors.down);
        rectangle('Position', [xs(1) + 0.15 * range(xs), ys(1) + 0.26 * range(ys), 0.10 * range(xs), 0.02 * range(ys)], ...
                  'EdgeColor', 'k', 'FaceColor', colors.up);
        rectangle('Position', [xs(1) + 0.25 * range(xs), ys(1) + 0.26 * range(ys), 0.10 * range(xs), 0.02 * range(ys)], ...
                  'EdgeColor', 'k', 'FaceColor', [0, 0, 0]);
        text(xs(1) + 0.36 * range(xs), ys(1) + 0.275 * range(ys), 'q=0.05 (FDR)');
    end
end



% -- Same figure plot: down and up -- 
figure('Position', [100, 100, 600, 400], 'Color','w');
set(gca, 'FontSize', fontSize);
for i0=1:numel(groups)
    grp=groups{i0};    
    
    taxis0=0:FRAME_DUR:FRAME_DUR*(size(avg_avgChg_F1_ML.(grp).down,2)-1);
    plot_sd_t(taxis0 * 1e3, avg_avgChg_F1_ML.(grp).down, sem_avgChg_F1_ML.(grp).down, colors.(grp), 'patch', 'FaceAlpha', 0.3);
    plot(taxis0 * 1e3, avg_avgChg_F1_ML.(grp).down, '-', 'Color', colors.(grp), 'LineWidth', 2);
    hold on;
    
    plot_sd_t(taxis0 * 1e3, avg_avgChg_F1_ML.(grp).up, sem_avgChg_F1_ML.(grp).up, colors.(grp), 'patch', 'FaceAlpha', 0.3);
    plot(taxis0 * 1e3, avg_avgChg_F1_ML.(grp).up, '-', 'Color', colors.(grp), 'LineWidth', 2);   
end
set(gca,'XLim',[taxis0(1),taxis0(end)] * 1e3,'YLim',[-0.06,0.06]);
xlabel('Time (ms)');    
ylabel('F1 change from baseline (normalized)');
% set(gca, 'YLim', [-0.04, 0.04]);
xs = get(gca, 'XLim'); ys = get(gca, 'YLim'); 
plot(xs, [0, 0], 'Color', [0.5, 0.5, 0.5]);
line(xs, repmat(ys(1) + 0.002 * range(ys), 1, 2), 'Color', 'k');
line(repmat(xs(1) + 0.002 * range(xs), 1, 2), ys, 'Color', 'k');

plot(xs(1) + [0.7, 0.8] * range(xs), repmat(ys(1) + 0.2 * range(ys), 1, 2), '-', 'Color', colors.PFS, 'LineWidth', 2);
text(xs(1) + 0.81 * range(xs), ys(1) + 0.2 * range(ys), 'PFS', 'Color', colors.PFS, 'FontSize', 12);
plot(xs(1) + [0.7, 0.8] * range(xs), repmat(ys(1) + 0.1 * range(ys), 1, 2), '-', 'Color', colors.PWS, 'LineWidth', 2);
text(xs(1) + 0.81 * range(xs), ys(1) + 0.1 * range(ys), 'PWS', 'Color', colors.PWS, 'FontSize', 12);

% -- Same figure plot: down vs. up -- 
figure('Position', [100, 50, 1000, 750], 'Color','w');
set(gca, 'FontSize', fontSize + 14);
for i0=1:numel(groups)
    grp=groups{i0};    
    
    taxis0=0:FRAME_DUR:FRAME_DUR*(size(avg_avgChg_F1_ML.(grp).down,2)-1);
    plot_sd_t(1e3 * taxis0, avg_avgChg_F1_ML.(grp).downUp, sem_avgChg_F1_ML.(grp).downUp, colors.(grp), 'patch', 'FaceAlpha', 0.3);
    plot(1e3 * taxis0, avg_avgChg_F1_ML.(grp).downUp, '-', 'Color', colors.(grp), 'LineWidth', 2);
    
    chgP = changePointTest(avg_avgChg_F1_ML.(grp).downUp);
    hold on;
end
% set(gca,'XLim',[taxis0(1),taxis0(end)],'YLim',[-0.06,0.06]);
xlabel('Time (ms)');    
ylabel('F1 difference between Down and Up (fract.)');
xs = get(gca, 'XLim'); ys = get(gca, 'YLim'); 
plot(xs, [0, 0], 'Color', [0.5, 0.5, 0.5]);
line(xs, repmat(ys(1) + 0.002 * range(ys), 1, 2), 'Color', 'k');
line(repmat(xs(1) + 0.002 * range(xs), 1, 2), ys, 'Color', 'k');

xs = get(gca,'XLim'); ys = get(gca, 'YLim');
draw_sgn_bar(1e3 * taxis0, p_t_twoGrp, 0.05, FDR_p_thresh_twoGrp, ...
             ys(1) + 0.09 * range(ys), 0.02 * range(ys), 'k', [0.75, 0.75, 0.75], 'k');
         
ezlegend([xs(1) + 0.05 * range(xs), ys(2) - 0.3 * range(ys), 0.75 * range(xs), 0.25 * range(ys)], ...
         0.45, {'PWS (mean \pm 1 SEM)', 'PFS (mean\pm 1 SEM)'}, {colors.PWS, colors.PFS}, ...
         [fontSize + 12, fontSize + 12], {'-', '-'}, {colors.PWS, colors.PFS}, ...
         [1, 1], [2, 2]);         
     

%% SD of F1 (variability of production and compensation)
figure('Position', [50, 250, 1350, 450]);
perts = {'none', 'down', 'up'};
for i0 = 1 : 3
    pert = perts{i0};
%     if isequal(pert, 'none')
%         pert = baseField;
%     end
    
    subplot('Position', [0.075 + (0.275 + 0.035) * (i0 - 1), 0.125, 0.275, 0.8]);
    set(gca, 'FontSize', fontSize);
    for i1 = 1 : numel(groups)
        grp = groups{i1};
        plot_sd_t(taxis0, avg_sd_F1_ML.(grp).(pert), sem_sd_F1_ML.(grp).(pert), ...
                  colors.(grp), 'patch', 'FaceAlpha', 0.3);
        hold on;
        plot(taxis0, avg_sd_F1_ML.(grp).(pert), '-', 'Color', colors.(grp));
    end
    
    set(gca, 'XLim', [taxis0(1), taxis0(end)]);
    set(gca, 'YLim', [20, 90]);
    
    p_grpComp = nan(N_SAMP, 1);
    for i1 = 1 : N_SAMP
        [h, p_grpComp(i1)] = ttest2(a_sd_F1_ML.PWS.(pert)(i1, :), a_sd_F1_ML.PFS.(pert)(i1, :));
    end
    FDR_p_thresh_grpComp = fdr(p_grpComp, FDR);
    
    xs = get(gca, 'XLim'); ys = get(gca, 'YLim');
    draw_sgn_bar(taxis0, p_grpComp, 0.05, FDR_p_thresh_grpComp, ...
             ys(1) + 0.09 * range(ys), 0.02 * range(ys), 'k', [0.75, 0.75, 0.75], 'k');
        
    title(strrep(pert, 'none', 'noPert'), 'FontSize', fontSize, 'FontWeight', 'Bold');
         
    grid on;
    if i0 == 2
        xlabel('Time (s)');
    end
    if i0 == 1
        ylabel('Within-subject SD of F1 (Hz) (Mean \pm 1 SEM)');
    end
    
    draw_xy_axes;
    
    if i0 == 1        
        ezlegend([xs(1) + 0.05 * range(xs), ys(2) - 0.3 * range(ys), 0.77 * range(xs), 0.25 * range(ys)], ...
                 0.45, {'PWS (mean \pm 1 SEM)', 'PFS (mean\pm 1 SEM)'}, {colors.PWS, colors.PFS}, ...
                 [fontSize - 1, fontSize - 1], {'-', '-'}, {colors.PWS, colors.PFS}, ...
                 [1, 1], [2, 2]);
    end
end

figure('Position', [50, 250, 900, 420]);
perts = {'down', 'up'};
for i0 = 1 : 2
    pert = perts{i0};
    subplot('Position', [0.075 + (0.425 + 0.035) * (i0 - 1), 0.125, 0.425, 0.8]);
    
    for i1 = 1 : numel(groups)
        grp = groups{i1};
        plot_sd_t(taxis0, avg_sdRatio_F1_ML.(grp).(pert), sem_sdRatio_F1_ML.(grp).(pert), ...
                  colors.(grp), 'patch', 'FaceAlpha', 0.3);
        hold on;
        plot(taxis0, avg_sdRatio_F1_ML.(grp).(pert), '-', 'Color', colors.(grp));
    end
    
    set(gca, 'XLim', [taxis0(1), taxis0(end)]);
    set(gca, 'YLim', [0.75, 1.25]);
    title([pert, ' / noPert'], 'FontSize', fontSize, 'FontWeight', 'Bold');
    
    p_grpComp = nan(N_SAMP, 1);
    for i1 = 1 : N_SAMP
        [h, p_grpComp(i1)] = ttest2(a_sdRatio_F1_ML.PWS.(pert)(i1, :), a_sdRatio_F1_ML.PFS.(pert)(i1, :));
    end
    FDR_p_thresh_grpComp = fdr(p_grpComp, FDR);
    xs = get(gca, 'XLim'); ys = get(gca, 'YLim');
    draw_sgn_bar(taxis0, p_grpComp, 0.05, FDR_p_thresh_grpComp, ...
             ys(1) + 0.09 * range(ys), 0.02 * range(ys), 'k', [0.75, 0.75, 0.75], 'k');
    
    xs = get(gca, 'XLim'); ys = get(gca, 'YLim');
    draw_sgn_bar(taxis0, p_grpComp, 0.05, FDR_p_thresh_grpComp, ...
             ys(1) + 0.09 * range(ys), 0.02 * range(ys), 'k', [0.75, 0.75, 0.75], 'k');
    
    plot(xs, [1, 1], '-', 'Color', [0.5, 0.5, 0.5]);     
    
    grid off;
    xlabel('Time (s)');
    if i0 == 1
        ylabel('Ratio SD of F1 between with and without perturbation (Hz) (Mean \pm 1 SEM)');
    end
    draw_xy_axes;
    
    if i0 == 1        
        ezlegend([xs(1) + 0.05 * range(xs), ys(2) - 0.225 * range(ys), 0.77 * range(xs), 0.2 * range(ys)], ...
                 0.45, {'PWS (mean \pm 1 SEM)', 'PFS (mean\pm 1 SEM)'}, {colors.PWS, colors.PFS}, ...
                 [fontSize - 1, fontSize - 1], {'-', '-'}, {colors.PWS, colors.PFS}, ...
                 [1, 1], [2, 2]);
    end
end

%% SD of F1: at specific time points
flds = {'none', 'down', 'up'};
tAxis = 0 : FRAME_DUR : FRAME_DUR * (size(a_sdRatio_F1_ML.PWS.down, 1) - 1);
for i1 = 1 : numel(grps)
    grp = grps{i1};
    a_F1_SD.(grp) = nan(numel(subjIDs.(grp)), 3);
    
    for i2 = 1 : numel(flds)
        fld = flds{i2};
        
        for i3 = 1 : numel(subjIDs.(grp))
            t_SDs = a_sd_F1_ML.(grp).(fld)(tAxis >= F1_SD_QUANT_TIME(1) & tAxis <= F1_SD_QUANT_TIME(2), i3);
        
            a_F1_SD.(grp)(i3, i2) = mean(t_SDs);
        end
    end
end

if ~isempty(fsic(varargin, 'Fig5XJitter'))
    xJitter = varargin{fsic(varargin, 'Fig5XJitter') + 1};
else
    xJitter = 0.;
end

figure('Position', [200, 200, 400, 300]);
set(gca, 'FontSize', 11);
for i1 = 1 : numel(grps)
    grp = grps{i1};
    if isequal(grp, 'PFS')
        symbl = 'o-';
    else
        symbl = 's-';
    end
    errorbar((1 : 3) + xJitter * (-1) ^ i1, mean(a_F1_SD.(grp)), ste(a_F1_SD.(grp)), ...
             symbl, 'Color', colors.(grp), 'LineWidth', 1.5);
    hold on;
end
set(gca, 'XLim', [0.5, 3.5]);
set(gca, 'XTick', [1, 2, 3], 'XTickLabel', {'noPert', 'down', 'up'});
ylabel('Within-subject SD of F1 (Hz) (Mean \pm 1 SEM)');

legend(grps);

do_RMAOV_1G1W(a_F1_SD, 0.05);

%% Subject-by-subject latencies
figure('Position', [100, 100, 800, 400], 'Color','w');
sph = nan(1, numel(groups));
yLims = nan(0, 2);
for i0 = 1 : numel(groups)
    grp = groups{i0};
    sph(i0) = subplot('Position',[0.125+(i0-1)*0.42,0.12,0.42,0.825]);
    
    for i1 = 1 : size(a_avgChg_traj_F1_ML.(grp).downUp, 2)
        t_sig = a_avgChg_traj_F1_ML.(grp).downUp(:, i1); 
        
        taxis0=0:FRAME_DUR:FRAME_DUR*(length(t_sig)-1);
        plot(taxis0, t_sig, '-', 'Color', colors.(grp));
        hold on;
        set(gca, 'XLim', [taxis0(1), taxis0(end)])
        
    end
    yLims = [yLims; get(gca, 'YLim')];
    xlabel('Time (s)');
    if i0 == 1
        ylabel('Composite response (normalized)');
    end
end

yLim = [min(yLims(:, 1)), max(yLims(:, 2))];
for i1 = 1 : numel(groups)
    set(gcf, 'CurrentAxes', sph(i1));
    set(gca, 'YLim', yLim);
end

%% Adaptation plot: mean 
figure('Position',[100,100,800,400], 'Name', 'F1 Adaptation Effects: Mean');
set(gca, 'FontSize', fontSize);
for i0=1:numel(groups)
    grp=groups{i0};
    subplot('Position',[0.125+(i0-1)*0.42,0.12,0.42,0.825]);
    
    taxis0=0:FRAME_DUR:FRAME_DUR*(size(avg_avgChg_F1_ML.(grp).down,2)-1);
    plot(taxis0,avg_adapt_avgChg_F1_ML.(grp).aftDown,'-','Color',colors.down,'LineWidth',2); hold on;
    plot(taxis0,avg_adapt_avgChg_F1_ML.(grp).aftDown-sem_adapt_avgChg_F1_ML.(grp).aftDown,'--','Color',colors.down);    
    plot(taxis0,avg_adapt_avgChg_F1_ML.(grp).aftDown+sem_adapt_avgChg_F1_ML.(grp).aftDown,'--','Color',colors.down);    
    plot(taxis0,avg_adapt_avgChg_F1_ML.(grp).aftUp,'-','Color',colors.up);
    plot(taxis0,avg_adapt_avgChg_F1_ML.(grp).aftUp-sem_adapt_avgChg_F1_ML.(grp).aftUp,'--','Color',colors.up);
    plot(taxis0,avg_adapt_avgChg_F1_ML.(grp).aftUp+sem_adapt_avgChg_F1_ML.(grp).aftUp,'--','Color',colors.up);
    
    plot([taxis0(1),taxis0(end)],[0,0],'-','color',[0.5,0.5,0.5]);
    
    if i0==2
        set(gca,'YTickLabel',{});
    else
        ylabel('Fraction F2 change from aftNone');
    end
    
    set(gca,'XLim',[taxis0(1),taxis0(end)],'YLim',[-0.06,0.06]);
    xlabel('Time (s)');
    title(grp);
end

%%
for i1 = 1 : numel(QUANTIFY_TIMES)    
    if numel(QUANTIFY_TIMES{i1}) == 1
        figure('Name',sprintf('QUANTIFY_TIME = %.3f s', QUANTIFY_TIMES{i1}), ...
               'Position', [200, 200, 400, 360]);
        set(gca, 'FontSize', fontSize);
    elseif numel(QUANTIFY_TIMES{i1}) == 2
        figure('Name',sprintf('QUANTIFY_TIME = [%.3f, %.3f] s', QUANTIFY_TIMES{i1}(1), QUANTIFY_TIMES{i1}(2)), ...
               'Position', [200, 200, 400, 360]);
        set(gca, 'FontSize', fontSize);
    end
    
    for i2 = 1 : numel(groups)
        grp = groups{i2};
        
        boxplot([upDownContra.PFS(:, i1); upDownContra.PWS(:, i1)], ...
                [1 * ones(size(upDownContra.PFS(:, i1))); 2 * ones(size(upDownContra.PWS(:, i1)))], ...
                'Colors', [0, 0, 0]);
        hold on;
        for i3 = 1 :  numel(subjIDs.(grp))
            plot(i2 - 0.25, upDownContra.(grp)(i3, i1), 'o', 'Color', colors.(grp));
%             text(i2+0.35, upDownContra.(grp)(i3, i1), strrep(strrep(subjIDs.(grp){i3}, '_1', ''), '_', '\_'), ...
%                 'Color', colors.(grp));
            hold on;
        end
        
        plot(i2 + 0.25, mean(upDownContra.(grp)), 'o', 'Color', colors.(grp), ...
             'LineWidth', 1);
        plot(repmat(i2 + 0.25, 1, 2), mean(upDownContra.(grp)) + ...
             [-1, 1] * ste(upDownContra.(grp)), '-', 'Color', colors.(grp), ...
             'LineWidth', 1);
    end
    
    [h, p] = ftest(upDownContra.PFS(:, i1), upDownContra.PWS(:, i1), 0.05);
    [h_t, p_t] = ttest2(upDownContra.PFS(:, i1), upDownContra.PWS(:, i1));
    [p_rs, h_rs] = ranksum(upDownContra.PFS(:, i1), upDownContra.PWS(:, i1));
    
%     xlabel('Vowel F1 difference limen (pct. perturbation)')
    ylabel('Ratio of compensation');
    set(gca, 'XTick', [1, 2], 'XTickLabel', {'PFS', 'PWS'});
    set(gca,'XLim',[0.5, 3 - 0.5]);
    
    xs = get(gca, 'XLim'); 
    ys = get(gca, 'YLim');
    text(xs(1) + 0.05 * range(xs), ys(2) - 0.075 * range(ys), sprintf('F-test: p = %.5f (one-tailed)', p), 'FontSize', fontSize);
    text(xs(1) + 0.05 * range(xs), ys(2) - 0.15 * range(ys), sprintf('t-test: p = %.5f (two-tailed)', p_t), 'FontSize', fontSize);
    text(xs(1) + 0.05 * range(xs), ys(2) - 0.225 * range(ys), sprintf('Wilcoxon Rank-sum test: p = %.5f (two-tailed)', p_rs), 'FontSize', fontSize);
end

%% Rationale for preserving the two PWS (PWS_M05 and PWS_F05) with mild hearing loss: JND comparison
figure('Position', [100, 100, 450, 320]);
plot(1, [vowelAcuity.PFS, vowelAcuity.PWS], 'bo'); 
hold on;
plot(2, vowelAcuity.PFS, 'o', 'Color', colors.PFS);
plot(3, vowelAcuity.PWS, 'o', 'Color', colors.PWS);

idx_S_1 = fsic(subjIDs.PWS, 'PWS_M05');
idx_S_2 = fsic(subjIDs.PWS, 'PWS_F05_1');
plot(4, vowelAcuity.PWS([idx_S_1, idx_S_2]), 'o', 'Color', colors.PWS);
set(gca, 'XLim', [0, 5]);
text(4 + 0.1, vowelAcuity.PWS(idx_S_1), 'PWS\_M05', 'Color', colors.PWS);
text(4 + 0.1, vowelAcuity.PWS(idx_S_2), 'PWS\_F05', 'Color', colors.PWS);

set(gca, 'XTick', [1 : 4], 'XTickLabel', {'All Ss', 'PFS', 'PWS', 'PWS with HL'});
ylabel('Vowel F1 difference limen (pct. of perturbation)');
set(gca, 'YLim', [0, 0.6]);
ys = get(gca, 'YLim');
set(gca, 'YLim', [ys(1) - 0.1 * range(ys), ys(2)]);
set(gca, 'YLim', [0, 0.6]);

%% Lag-N analysis and early-following response
figure('Position', [50, 100, 800, 380], 'Name', 'Lag-N mean F1 changes (Hz)');
for i0 = 1 : 2
    subplot(1, 2, i0);
    set(gca, 'FontSize', fontSize);
    if i0 == 1  pert = 'down';
    else        pert = 'up';
    end

    for i1 = 1 : numel(groups)
        grp = groups{i1};
        errorbar(mean(lagN_F1Chg.(grp).(pert)), ste(lagN_F1Chg.(grp).(pert)), 'Color', colors.(grp));
        set(gca, 'YLim', [-15, 15]);
        
        hold on;
    end
    
    ys = get(gca, 'YLim');
    nLag = size(lagN_F1Chg.PFS.down, 2);
    p_t_twoGrp = nan(1, size(lagN_F1Chg.PFS.down, 2));
    for i1 = 1 : nLag
        [h, p_t_twoGrp(i1)] = ttest2(lagN_F1Chg.PFS.(pert)(:, i1), lagN_F1Chg.PWS.(pert)(:, i1));
        
        if p_t_twoGrp(i1) < 0.05
            plot(i1, ys(2) - 0.05 * range(ys), 'k*');
        end
    end
        
    set(gca, 'XTick', 1 : MAX_LAG, 'XLim', [0, MAX_LAG + 1]);
    xlabel('Lag N');
    ylabel('F1 Change (Hz)');
    
    
    xs = get(gca, 'XLim');
    plot(xs, [0, 0], '-', 'Color', [0.5, 0.5, 0.5]  );
    
    title(['After ', pert]);
    
    if i0 == 1
        legend(groups, 'Location', 'Southwest');
    end
end

lagN_F1Chg_down.PWS = lagN_F1Chg.PWS.down;
lagN_F1Chg_up.PWS = lagN_F1Chg.PWS.up;
lagN_F1Chg_down.PFS = lagN_F1Chg.PFS.down;
lagN_F1Chg_up.PFS = lagN_F1Chg.PFS.up;

do_RMAOV_1G1W(lagN_F1Chg_down, 0.05);
do_RMAOV_1G1W(lagN_F1Chg_up, 0.05);

%% Comparison of the vowel acuity from the two groups
vowelAcuity_1.PFS = vowelAcuity.PFS(~isnan(vowelAcuity.PFS));
vowelAcuity_1.PWS = vowelAcuity.PWS(~isnan(vowelAcuity.PWS));
metaPlot_2grp(vowelAcuity_1, 'PFS', 'PWS', 'F2 perturbation: Up - Down', 'JND of vowel F1 (Frac. pert.)', colors, 'zeroLine');

figure('Position', [100, 50, 350, 350]);
set(gca, 'FontSize', fontSize + 2);
for i1 = 1 : 2
    if i1 == 1; grp = 'PFS'; 
    else; grp = 'PWS';
    end
    bar(i1, nanmean(vowelAcuity.(grp)), 'EdgeColor', colors.(grp), 'FaceColor', 'none', 'LineWidth', 2);
    hold on;
    plot([i1, i1], nanmean(vowelAcuity.(grp)) + [-1, 1] * nanste(vowelAcuity.(grp)), '-', 'Color', colors.(grp), 'LineWidth', 2);
end
set(gca, 'XTick', [1, 2], 'XTickLabel', {'PFS', 'PWS'});
ylabel('JND of F1 (Fraction of pert.) (mean \pm 1 SEM)');
xs = get(gca, 'XLim'); ys = get(gca, 'YLim');
[h_t, p_t] = ttest2(vowelAcuity.PFS, vowelAcuity.PWS);
text(xs(1) + 0.05 * range(xs), ys(2) - 0.055 * range(ys), sprintf('t-test: p=%.4f', p_t), 'FontSize', fontSize + 1);
box off;


%% Relation between vowel acuity and compensation to perturbation
figure('name', 'Vowel acuity - compensation relation', 'Position', [200, 200, 500, 360]);
set(gca, 'FontSize', fontSize + 1);
for i1 = 1 : numel(QUANTIFY_TIMES)    
    subplot(1, 1, i1);
    title(sprintf('QUANTIFY_TIME = %.3f s', QUANTIFY_TIMES{i1}));
    
    t_x = [];
    t_y = [];
    
    for i2 = 1 : numel(groups)
        grp = groups{i2};
        
        for i3 = 1 :  numel(subjIDs.(grp))
            if isequal(grp, 'PWS')
                symbl = 's';
            else
                symbl = 'o';
            end
            
            plot(vowelAcuity.(grp)(i3), upDownContra.(grp)(i3,i1), symbl, 'Color', colors.(grp));
            hold on;
            
            t_x(end + 1) = vowelAcuity.(grp)(i3);
            t_y(end + 1) = upDownContra.(grp)(i3,i1);
%             plot(vowelAcuity.PWS, upDownContra.PWS(:,i1), 'ro');
        end
        
    end
    xlabel('JND to vowel F1 change (fraction of perturbation)');
    ylabel('Ratio of compensation');
    
    [k,r2,p]=lincorr(t_x,t_y);
    [rho_spear,t_spear,p_spear] = spear(t_x', t_y');
    xs=get(gca,'XLim'); ys = get(gca,'YLim');
    plot(xs, xs *k(2) + k(1), '--', 'Color', [0.5, 0.5, 0.5]);
    text(xs(1) + 0.025 * range(xs), ys(2) - 0.06 * range(ys), sprintf('Linear correlation: R^2 = %.4f, p = %.4f', r2, p), 'FontSize', fontSize);
    text(xs(1) + 0.025 * range(xs), ys(2) - 0.12 * range(ys), sprintf('Spearman: p = %.4f',p_spear), 'FontSize', fontSize);
    
    %%%%
    xs = get(gca, 'XLim'); ys = get(gca, 'YLim');
    ezlegend([xs(1) + 0.65 * range(xs), ys(1) + 0.025 * range(ys), 0.25 * range(xs), 0.15 * range(ys)], ...
              0.5, {'PWS', 'PFS'}, ...
              {colors.PWS, colors.PFS}, ...
              repmat(fontSize, 1, 2), {'o', 'o'}, {colors.PWS, colors.PFS}, ...
              [1, 1], [0, 0]);
end

%% Relation between vowel acuity and adaptation to perturbation
figure('Name',sprintf('QUANTIFY_TIME_ADAPT = %.3f - %.3f s', QUANTIFY_TIME_ADAPT(1), QUANTIFY_TIME_ADAPT(2)));
set(gca, 'FontSize', fontSize);
t_x = []; 
t_y = [];
for i2 = 1 : numel(groups)
    grp = groups{i2};

    for i3 = 1 :  numel(subjIDs.(grp))
        plot(vowelAcuity.(grp)(i3) * 1e2, aftUpDownContra.(grp)(i3), 'o', 'Color', colors.(grp));
        t_x(end + 1) = vowelAcuity.(grp)(i3);
        t_y(end + 1) = aftUpDownContra.(grp)(i3);
        
        
        hold on;
    end
end
idx = find(~isnan(t_x) & ~isnan(t_y));
t_x = t_x(idx);
t_y = t_y(idx);
[k,r2,p]=lincorr(t_x,t_y);
[rho_spear,t_spear,p_spear] = spear(t_x', t_y');
xs=get(gca,'XLim'); ys = get(gca,'YLim');
plot(xs, xs *k(2) + k(1), '--', 'Color', [0.5, 0.5, 0.5]);
text(xs(1) + 0.025 * range(xs), ys(2) - 0.06 * range(ys), sprintf('Lincorr: p = %.4f',p));
text(xs(1) + 0.025 * range(xs), ys(2) - 0.12 * range(ys), sprintf('Spearman: p = %.4f',p_spear));
set(gca,'XLim',xs,'YLim',ys);
xlabel('Vowel F1 difference limen (pct. perturbation)');
ylabel('Ratio of compensation');

%% Correlation between SSI4 and compensation to perturbation
SSI4_fields = {'total', 'freq', 'dur', 'concom'};
for i1 = 1 : numel(SSI4_fields)
    field = SSI4_fields{i1};
    figure('Name', sprintf('Correlation between SSI4 and compensation in the PWS group'), ...
           'Position', [150, 150, 400, 300]);
    set(gca, 'FontSize', fontSize);
    plot(SSI4.(field).PWS, upDownContra.PWS(:, end), 'o', 'Color', colors.(grp));
    xlabel(sprintf('SSI4 - %s', field));
    ylabel('Ratio of compensation');
    idx = find(~isnan(SSI4.(field).PWS));
    [k,r2,p]=lincorr(SSI4.(field).PWS(idx), upDownContra.PWS(idx, end));
    xs = get(gca, 'XLim'); ys = get(gca, 'YLim'); 
    text(xs(1) + 0.05 * range(xs), ys(2) - 0.05 * range(ys), ...
        sprintf('Lin. corr.: r^2 = %f; p = %f', r2, p), 'FontSize', fontSize);
    [rho, t, p_s] = spear(SSI4.(field).PWS(idx), upDownContra.PWS(idx, end));
    text(xs(1) + 0.05 * range(xs), ys(2) - 0.10 * range(ys), ...
        sprintf('Spearman: rho = %f, p = %f', rho, p_s), 'FontSize', fontSize);
    hold on;
    plot(xs, k(1) + k(2) * xs, '--', 'Color', colors.PWS);
end

%% sdRatio: for testing the hypothesis of more unstable inverse model in PWS
figure;
set(gca, 'FontSize', fontSize);
errorbar(1 : 3, [mean(sdF1.PWS.(baseField)(:, end)), mean(sdF1.PWS.down(:, end)), mean(sdF1.PWS.up(:, end))], ...
         [ste(sdF1.PWS.(baseField)(:, 1)), ste(sdF1.PWS.down(:, 1)), ste(sdF1.PWS.up(:, 1))], 'o-', 'Color', colors.PWS);
hold on;
errorbar(1 : 3, [mean(sdF1.PFS.(baseField)(:, end)), mean(sdF1.PFS.down(:, end)), mean(sdF1.PFS.up(:, end))], ...
         [ste(sdF1.PFS.(baseField)(:, 1)), ste(sdF1.PFS.down(:, 1)), ste(sdF1.PFS.up(:, 1))], 'o-', 'Color', colors.PFS);
set(gca, 'XTick', [0 : 4], 'XTickLabel', {baseField, 'down', 'up'});
set(gca, 'XLim', [0, 4]);
ylabel(sprintf('sdF1 between %.2f and %.2f s re. perturb. onset', QUANTIFY_TIMES{end}(1), QUANTIFY_TIMES{end}(2)));

figure;
set(gca, 'FontSize', fontSize);
errorbar(1 : numel(QUANTIFY_TIMES), mean(sdRatio.PWS), ste(sdRatio.PWS), 'o-', 'Color', colors.PWS);
hold on;
errorbar(1 : numel(QUANTIFY_TIMES), mean(sdRatio.PFS), ste(sdRatio.PFS), 'o-', 'Color', colors.PFS);
xlabel('QUANTIFY TIME');
ylabel('sdRatio');
set(gca, 'XTick', [1 : numel(QUANTIFY_TIMES)]);

%% Write the sd data to an xls file for later use in SYSTAT
sdF1_cell.PWS = [sdF1.PWS.(baseField)(:, end), sdF1.PWS.down(:, end), sdF1.PWS.up(:, end)];
sdF1_cell.PFS = [sdF1.PFS.(baseField)(:, end), sdF1.PFS.down(:, end), sdF1.PFS.up(:, end)];
a_sdF1_cell = [num2cell([0 * ones(numel(subjIDs.PWS), 1), sdF1_cell.PWS]); ...
               num2cell([1 * ones(numel(subjIDs.PFS), 1), sdF1_cell.PFS])];
a_sdF1_cell = [{'GRP', 'None', 'Down', 'Up'}; a_sdF1_cell];

t_xls_fn = fullfile(xlsSysDir, 'sdF1');
status=xlswrite(t_xls_fn, a_sdF1_cell);
if status==1
    fprintf('sdF1 data successfully written to %s.\n',t_xls_fn);    
else
    fprintf('Writing sdF1 data to %s was unsuccessful.\n',t_xls_fn);   
end

%% Stats related to vowelOnset correction
% for i1 = 1 : numel(groups)
%     grp = groups{i1};
%     fprintf('%f%% of all trials from %s contain vowelOnset manual correction.\n', ...
%         numel(find(a_vowelOnset_manCorr.all.(grp)) ~= 0) / numel(a_vowelOnset_manCorr.all.(grp)) * 100, grp);
%     fprintf('\tMedian correction: %f s\n', median(a_vowelOnset_manCorr.all.(grp)));
%     fprintf('\tMean +/- SD correction: %f +/- %f s\n', mean(a_vowelOnset_manCorr.all.(grp)), std(a_vowelOnset_manCorr.all.(grp)));
% end
% 
% fprintf('\n');
% for i1 = 1 : numel(groups)
%     grp = groups{i1};
%     fprintf('%f%% of pert trials from %s contain vowelOnset manual correction.\n', ...
%         numel(find(a_vowelOnset_manCorr.pert.(grp)) ~= 0) / numel(a_vowelOnset_manCorr.pert.(grp)) * 100, grp);
%     fprintf('\tMedian correction: %f s\n', median(a_vowelOnset_manCorr.pert.(grp)));
%     fprintf('\tMean +/- SD correction: %f +/- %f s\n', mean(a_vowelOnset_manCorr.pert.(grp)), std(a_vowelOnset_manCorr.pert.(grp)));
% end



%% Stats
for i1 = 1 : numel(QUANTIFY_TIMES)
    [h_omni, p_omni] = ttest([upDownContra.PFS(:, i1); upDownContra.PWS(:, i1)]);
    if numel(QUANTIFY_TIMES{i1}) == 1
        fprintf('Change at %f s following perturbation onset:\n\tp_omni = %f\n', ...
            QUANTIFY_TIMES{i1}, p_omni);
    else
        fprintf('Change within [%f, %f] s following perturbation onset:\n\tp_omni = %f\n', ...
            QUANTIFY_TIMES{i1}(1), QUANTIFY_TIMES{i1}(2), p_omni);
    end

    [h_PFS, p_PFS] = ttest(upDownContra.PFS(:, i1));
    fprintf('\tp_PFS = %f\n', p_PFS);

    [h_PWS, p_PWS] = ttest(upDownContra.PWS(:, i1));
    fprintf('\tp_PWS = %f\n', p_PWS);

    [h_t, p_t, ci, t_stats] = ttest2(upDownContra.PFS(:, i1),  upDownContra.PWS(:, i1));
    fprintf('\tBetween-group t-test: t(%d) = %f, p = %f\n', t_stats.df, t_stats.tstat, p_t);
    
    cohend = cohen_d(upDownContra.PFS(:, i1),  upDownContra.PWS(:, i1));
    fprintf('\tBetween-group Cohen''s d = %f\n', cohend);
    
    [p_rs, h_rs] = ranksum(upDownContra.PFS(:, i1),  upDownContra.PWS(:, i1));
    fprintf('\tBetwen-group rank-sum test: p = %f\n\n', p_rs);

end

%% Correlation with FDA-SVI
load(FDA_dataSet_FN);
a_isPWS = [];
a_subjIDs = {};
a_FDA_SVI = [];
a_sdF1 = [];

for i1 = 1 : numel(groups)
    grp = groups{i1};
    
    for i2 = 1 : numel(subjIDs.(grp))
        subjID = subjIDs.(grp){i2};
        
        if isequal(subjID(end - 1 : end), '_1')
            subjID = subjID(1 : end - 2);
        end
        if ~isempty(fsic(grpData.subjIDs.(grp), subjID))
            idx_fda = fsic(grpData.subjIDs.(grp), subjID); 
            
            a_isPWS(end + 1) = isequal(grp, 'PWS');
            a_subjIDs{end + 1} = subjID;
            a_FDA_SVI(end + 1) = mean(grpData.fdaSV.(grp)(idx_fda, :));
            a_sdF1(end + 1) = sdF1.(grp).none(i2, end); 
        end
    end
end

figure;
set(gca, 'FontSize', fontSize)
plot(a_FDA_SVI(a_isPWS == 1), a_sdF1(a_isPWS == 1), 'o', 'Color', colors.PWS);
hold on;
plot(a_FDA_SVI(a_isPWS == 0), a_sdF1(a_isPWS == 0), 'o', 'Color', colors.PFS);
xlabel('FDA-SVI');
ylabel('sdF1'); 

%% Correlation between variability under noPert and compensation
figure;
plot(sdF1.PWS.none(:, 1), upDownContra.PWS(:, end), 'o', 'Color', colors.PWS);
hold on;
plot(sdF1.PFS.none(:, 1), upDownContra.PFS(:, end), 'o', 'Color', colors.PFS);
xlabel('SD of F1 under noPert (Hz)');
ylabel('Ratio of compensation');


return