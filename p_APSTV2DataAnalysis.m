function varargout = p_APSTV2DataAnalysis(varargin)
%% SUBJECT LISTS
SUBJ_LIST.apstv   = get_APSTV_subjIDs('apstv');
SUBJ_LIST.apstv2  = get_APSTV_subjIDs('apstv2');
SUBJ_LIST.apstv2t = get_APSTV_subjIDs('apstv2t');
                 
SUBJ_LIST.apstv2_stut_s_PFS = get_STUT_ST_subjIDs('PFS', 'S');
SUBJ_LIST.apstv2_stut_s_PWS = get_STUT_ST_subjIDs('PWS', 'S');                 

SUBJ_LIST.apstv2_stut_t_PFS = get_STUT_ST_subjIDs('PFS', 'T');
SUBJ_LIST.apstv2_stut_t_PWS = get_STUT_ST_subjIDs('PWS', 'T');

%% Config
MAX_TRIALS_PER_REP = 8;
N_REPS = 20;

FTN_mode = 'linear'; % 'linear', 'nearest', or 'spline'
if ~isequal(FTN_mode, 'linear')
    bGo = input(sprintf('Using non-linear interpolation method: %s. Are you sure you want to proceed? (0/1): ', FTN_mode));
else
    bGo = 1;
end
if bGo == 0
    return
end

if ~isempty(fsic(varargin,'apstv'))
    subjIDs=SUBJ_LIST.apstv;
    studyMode='apstv';    
elseif ~isempty(fsic(varargin,'apstv2'))
    subjIDs=SUBJ_LIST.apstv2;
    studyMode='apstv2';
elseif ~isempty(fsic(varargin,'apstv2t'))
    subjIDs=SUBJ_LIST.apstv2t;
    studyMode='apstv2t';
elseif ~isempty(fsic(varargin,'apstv2_stut_s.PFS'))
    subjIDs=SUBJ_LIST.apstv2_stut_s_PFS;
    studyMode = 'apstv2_stut_s';   % CAUTION
elseif ~isempty(fsic(varargin,'apstv2_stut_s.PWS'))
    subjIDs=SUBJ_LIST.apstv2_stut_s_PWS;
    studyMode = 'apstv2_stut_s';   % CAUTION    
elseif ~isempty(fsic(varargin,'apstv2_stut_t.PFS'))
    subjIDs=SUBJ_LIST.apstv2_stut_t_PFS;
    studyMode = 'apstv2_stut_t';   % CAUTION
elseif ~isempty(fsic(varargin,'apstv2_stut_t.PWS'))
    subjIDs=SUBJ_LIST.apstv2_stut_t_PWS;
    studyMode = 'apstv2_stut_t';   % CAUTION
elseif ~isempty(fsic(varargin,'apstv_expanded'))
    subjIDs=[SUBJ_LIST.apstv, SUBJ_LIST.apstv2_stut_s_PFS];
    studyMode='apstv';
elseif ~isempty(fsic(varargin,'apstv2t_expanded'))
    subjIDs=[SUBJ_LIST.apstv2t, SUBJ_LIST.apstv2_stut_t_PFS];
    studyMode='apstv2t';
end

fontSize=12;

config.normT1_bw_iu_ratio=[0,1/4,2/4,3/4,1];
config.normT1_bw_uy_ratio=[1/4,2/4,3/4,1];
config.normT2_bw_uy_ratio=[1/4,2/4,3/4,1];
config.normT3_bw_yu_ratio=[1/4,2/4,3/4,1];

colors.none='k';
colors.accel=[1,0.5,0];
colors.decel=[0,0.5,0];
colors.up=[1,0,0];
colors.down=[0,0.4,1];
% colors.down='b';
% colors.up='r';
% colors.down=[0.3,0.3,0.3];
% colors.up=[0.6,0.6,0.6];

doPOA=0;
doIOA=0;
doPOA_TN=0;
doIOA_TN=0;
toWriteXls=0;

opt_OA='';
if ~isempty(fsic(varargin,'showPert'))
    opt_OA='showPert';
end

frameDur=16/12000;

posthocOpt='tukey';         % 'tukey', and 'lsd', 'scheffe' (not yet implemented), 'fwe'
multiSubOpt='';
if ~isempty(fsic(varargin,'noMultiSub'))
    multiSubOpt='noMultiSub';
end

opt_avg='all_trajs';
if ~isempty(fsic(varargin,'ind_avg_trajs'))
    opt_avg='ind_avg_trajs';
end

opt_pertBar='pertBar';
if ~isempty(fsic(varargin,'showPert')) || ~isempty(fsic(varargin,'noPertBar'))
    opt_pertBar='';
end

% ana_winWidth=280;
% ana_winHeight=200;
if isequal(studyMode,'apstv') || isequal(studyMode,'apstv2_stut_s')
    ana_winWidth=160;
    ana_winHeight=180;
else
    ana_winWidth=160;
    ana_winHeight=150;
end
if ~isempty(fsic(varargin,'ana_winWidth'))
    ana_winWidth=varargin{fsic(varargin,'ana_winWidth')+1};
end
if ~isempty(fsic(varargin,'ana_winHeight'))
    ana_winHeight=varargin{fsic(varargin,'ana_winHeight')+1};
end

if ~isempty(fsic(varargin,'doIOA'))
    doIOA=varargin{fsic(varargin,'doIOA')+1};
end


%% Custom number of subjects
if ~isempty(fsic(varargin, 'nSubjs'))
    nSubjs = varargin{fsic(varargin, 'nSubjs') + 1};
    subjIDs = subjIDs(1 : nSubjs);
end



%% Space/slot creation 
protoStructArray=struct; 
if isequal(studyMode,'apstv') || isequal(studyMode,'apstv2_stut_s')
    protoStructArray.up=[];
    protoStructArray.down=[];
    protoStructArray.none=[];
    
    protoStructArray_AT.aft_up=[];
    protoStructArray_AT.aft_down=[];
    protoStructArray_AT.aft_none=[];
elseif isequal(studyMode,'apstv2')
    protoStructArray.accel=[];
    protoStructArray.decel=[];
    protoStructArray.up=[];
    protoStructArray.down=[];
    protoStructArray.none=[];
    
    protoStructArray_AT.aft_up=[];
    protoStructArray_AT.aft_down=[];
    protoStructArray_AT.aft_none=[];
    protoStructArray_AT.aft_accel=[];
    protoStructArray_AT.aft_decel=[];
elseif isequal(studyMode,'apstv2t') || isequal(studyMode,'apstv2_stut_t')
    protoStructArray.accel=[];
    protoStructArray.decel=[];
    protoStructArray.none=[];
    
    protoStructArray_AT.aft_none=[];
    protoStructArray_AT.aft_accel=[];
    protoStructArray_AT.aft_decel=[];
end
changeYF2=protoStructArray;
changeIF2=protoStructArray;
changeIUInt=protoStructArray;
changeUYInt=protoStructArray;
changeIYInt=protoStructArray;
changeIU2Int=protoStructArray;
changeIY2Int=protoStructArray;
changeIU3Int=protoStructArray;
changeIY3Int=protoStructArray;

a_IUInt=[];         chg_IUInt=[];
a_IYInt=[];         chg_IYInt=[];

a_iF2=[];
a_uF2=[];           chg_uF2=[];
a_yF2=[];           chg_yF2=[];
a_uyMidF2=[];
a_yu2MidF2=[];
a_u2F2=[];
a_y2F2=[];
a_u3F2=[];
a_y3F2=[];
a_IU2Int=[];        chg_IU2Int=[];
a_IY2Int=[];        chg_IY2Int=[];
a_IU3Int=[];        chg_IU3Int=[];
a_IY3Int=[];        chg_IY3Int=[];
a_I_UYMid_Int=[];   chg_I_UYMid_Int=[];
a_I_YU2Mid_Int=[];  chg_I_YU2Mid_Int=[];

a_Y1Y2Int=[];       chg_Y1Y2Int=[];
a_Y2Y3Int=[];       chg_Y2Y3Int=[]; 

a_I_UYMidF_Int=[];  chg_I_UYMidF_Int=[];
a_I_YU2MidF_Int=[]; chg_I_YU2MidF_Int=[];

a_tShift_iuMidF.accel=[];
a_tShift_iuMidF.decel=[];

change_bw_iu_f2=cell(1,length(config.normT1_bw_iu_ratio));
change_bw_uy_f2=cell(1,length(config.normT1_bw_uy_ratio));
change_bw_yu_f2=cell(1,length(config.normT3_bw_yu_ratio));
normT1_bw_iu_f2=cell(1,length(config.normT1_bw_iu_ratio));
% normT1_bw_uy_f2=cell(1,length(config.normT1_bw_uy_ratio));
normT2_bw_uy_f2=cell(1,length(config.normT2_bw_uy_ratio));
normT3_bw_yu_f2=cell(1,length(config.normT3_bw_yu_ratio));

for i1=1:length(config.normT1_bw_iu_ratio)
    change_bw_iu_f2{i1}=[];
    normT1_bw_iu_f2{i1}=[];
end
for i1=1:length(config.normT2_bw_uy_ratio)
    change_bw_uy_f2{i1}=[];
    normT2_bw_uy_f2{i1}=[];
end
for i1=1:length(config.normT3_bw_yu_ratio)
    change_bw_yu_f2{i1}=[];
    normT3_bw_yu_f2{i1}=[];
end

trajF2_POA=protoStructArray;
velF2_POA=protoStructArray;
trajF1_POA=protoStructArray;
velF1_POA=protoStructArray;

trajF2_POA_TN=protoStructArray; % Time normalized POA: normalized by the average onset-yF2 duration
trajF1_POA_TN=protoStructArray;

trajF2_IOA=protoStructArray;
velF2_IOA=protoStructArray;
trajF1_IOA=protoStructArray;
velF1_IOA=protoStructArray;
trajF2_IOA_TN=protoStructArray;
trajF1_IOA_TN=protoStructArray;

trajF2_IOA_NN=protoStructArray; % non-normalized

trajF2_pert_IOA=protoStructArray;
trajF2_pert_IOA_FN=protoStructArray;

trajF2_POA_TN2=protoStructArray;
trajF1_POA_TN2=protoStructArray;

chgTrajF2_POA=protoStructArray;
chgTrajF1_POA=protoStructArray;
chgTrajF2_POA_TN=protoStructArray;
chgTrajF1_POA_TN=protoStructArray;
chgTrajF2_POA_TN_FN=protoStructArray;
chgTrajF1_POA_TN_FN=protoStructArray;
chgTrajF2_IOA=protoStructArray;     

chgTrajF2_POA.contrast = {};

semChgTrajF2_IOA = protoStructArray;
semChgTrajF2_FTN = protoStructArray;

chgTrajF1_IOA=protoStructArray;
chgTrajF2_IOA_FN=protoStructArray;
chgTrajF1_IOA_FN=protoStructArray;
chgTrajF2_IOA_TN=protoStructArray;
chgTrajF1_IOA_TN=protoStructArray;
chgTrajF2_IOA_TN_FN=protoStructArray;
chgTrajF1_IOA_TN_FN=protoStructArray;

chgTrajF2_IOA.contrast = {};

pertShiftF2_IOA=protoStructArray;
pertShiftF1_IOA=protoStructArray;
pertShiftF2_IOA_TN=protoStructArray;
pertShiftF1_IOA_TN=protoStructArray;
pertShiftF2_IOA_TN_FN=protoStructArray;
pertShiftF1_IOA_TN_FN=protoStructArray;

chgTrajF2_FTN=protoStructArray; % FTN stands for fullTNorm, or full time-normalization
pertShiftF2_FTN=protoStructArray; % FTN stands for fullTNorm, or full time-normalization

chgTrajF2_FTN.contrast = {};

avgTrajF2_FTN=protoStructArray;
avgTrajF2_FTN_FN=protoStructArray;
avgTrajF2_pert_FTN_FN=protoStructArray;

% fe_trajF2_FTN=protoStructArray;  % fe stands for "fixed effects"
fe_chgTrajF2_FTN=protoStructArray;

chgTrajF2_FTN_TA=protoStructArray_AT; % TA stands for fullTNorm

ioaTraj_tInts = protoStructArray;
ioaTraj_f2Vals = protoStructArray;

trajRMS_IOA = protoStructArray;

%% Conformation vectors: conformation of individual data patterns to pop. trend
%   -1: significantly counter-trend
%   -0.5: non-significantly counter-trend
%   0.5: non-significantly trend
%   1: significantly trend
if isequal(studyMode, 'apstv') || isequal(studyMode, 'apstv2_stut_s')
    cfmVect.iuInt = nan(numel(subjIDs), 3); % noPert-Down, Up-noPert, Up-Down
    cfmVect.uF2 = nan(numel(subjIDs), 3);   % Down-noPert, noPert-Up, Down-Up
    cfmVect.uyMidF2 = nan(numel(subjIDs), 3);
    cfmVect.yF2 = nan(numel(subjIDs), 3);   % Down-noPert, noPert-Up, Down-Up   
    cfmVect.y2u2MidF2 = nan(numel(subjIDs), 3);
    
    ssMult.iuInt = nan(numel(subjIDs), 3);  % Sample-size multiplier
    ssMult.uF2 = nan(numel(subjIDs), 3);
    ssMult.yF2 = nan(numel(subjIDs), 3);
elseif isequal(studyMode, 'apstv2t') || isequal(studyMode, 'apstv2_stut_t')
    cfmVect.iuInt = nan(numel(subjIDs), 3);
    cfmVect.iyInt = nan(numel(subjIDs), 3);
    
    ssMult.iuInt = nan(numel(subjIDs), 3);
    ssMult.iyInt = nan(numel(subjIDs), 3);
end
cfmVect_byRep = cfmVect;

%% Data extraction 
nTotTrials = 0;
nDscd_all = 0;
nDscd_prodErr = 0;

for i1=1:length(subjIDs)    
    load(fullfile(get_dacacheDir(subjIDs{i1},SUBJ_LIST), [subjIDs{i1},'.mat']));        % gives pdata
    fprintf('Loading data from subject %s\n',subjIDs{i1});
    
%     if isequal(studyMode,'apstv2');
    if isequal(studyMode,'apstv');
        pdata=format_apstv_pdata(pdata);
    end
    
    t_IF2=protoStructArray;
    t_UF2=protoStructArray;
    t_YF2=protoStructArray;
    t_IF2=protoStructArray; 
    t_U2F2=protoStructArray;
    t_U2Y2MidF2=protoStructArray;
    t_Y2F2=protoStructArray;
    t_U3F2=protoStructArray;
    t_Y3F2=protoStructArray;
    
    t_IUInt=protoStructArray;
    t_UYInt=protoStructArray;
    t_IYInt=protoStructArray;
    t_IU2Int=protoStructArray;
    t_IY2Int=protoStructArray;
    t_IU3Int=protoStructArray;
    t_IY3Int=protoStructArray;
    
    t_I_UYMid_Int=protoStructArray;
    t_I_YU2Mid_Int=protoStructArray;
    
    t_I_UYMidF_Int=protoStructArray;
    t_I_YU2MidF_Int=protoStructArray;
    
%     t_trajF2_POA_TN2=protoStructArray;
    t_trajF2=protoStructArray;
    t_trajF1=protoStructArray;
    
    flds=fields(t_YF2);
    
    if ~isempty(fsic(flds, 'down')) && isempty(fsic(flds, 'accel'))
        flds_cntrst = {'down', 'up'};
    elseif ~isempty(fsic(flds, 'accel')) && isempty(fsic(flds, 'down'))
        flds_cntrst = {'accel', 'decel'};
    end
    
    for j1 = 1 : numel(flds)
        m_UF2.(flds{j1}) = nan(MAX_TRIALS_PER_REP, N_REPS);
        m_YF2.(flds{j1}) = nan(MAX_TRIALS_PER_REP, N_REPS);
        m_IUInt.(flds{j1}) = nan(MAX_TRIALS_PER_REP, N_REPS);
        m_IYInt.(flds{j1}) = nan(MAX_TRIALS_PER_REP, N_REPS);
    end        
    
    for i2=1:length(flds)
        fld=flds{i2};
        
        if isequal(studyMode,'apstv') || isequal(studyMode,'apstv2_stut_s')
            if isequal(fld,'accel') || isequal(fld,'decel')
                continue;
            end
        elseif isequal(studyMode,'apstv2t') || isequal(studyMode,'apstv2_stut_t')
            if isequal(fld,'down') || isequal(fld,'up')
                continue;
            end
        end
        
        for i3=1:length(pdata.utters.(flds{i2}))
            if isempty(pdata.utters.(flds{i2}){i3})
                continue;
            end
            
            nTotTrials = nTotTrials + 1;
            
            bContinue = 0;
            if ~isfield(pdata.utters.(flds{i2}){i3},'bDiscard')
                bContinue = 1;                
            end
            if isfield(pdata.utters.(flds{i2}){i3},'bDiscard') 
                if pdata.utters.(flds{i2}){i3}.bDiscard==1
                    bContinue = 1;
                end
            end
            if isfield(pdata.utters.(flds{i2}){i3},'rating') && pdata.utters.(flds{i2}){i3}.rating==0
                bContinue = 1;
            end
            if isfield(pdata.utters.(flds{i2}){i3},'ratingComments') && ~isempty(strfind(pdata.utters.(flds{i2}){i3}.ratingComments,'Discarded.'))
                bContinue = 1;
            end
            
            if isfield(pdata.utters.(flds{i2}){i3}, 'bPertOkay')
                if pdata.utters.(flds{i2}){i3}.bPertOkay==0 || isequal(pdata.utters.(flds{i2}){i3}.bPertOkay,'0')                   
                    bConitnue = 1;
                end
            end
            
            if bContinue == 1
                nDscd_all = nDscd_all + 1;
                if isfield(pdata.utters.(flds{i2}){i3},'ratingComments') && ...
                   (~isempty(strfind(pdata.utters.(flds{i2}){i3}.ratingComments,'ysfluency')) || ...
                    ~isempty(strfind(pdata.utters.(flds{i2}){i3}.ratingComments,'speech error')) || ...
                    ~isempty(strfind(pdata.utters.(flds{i2}){i3}.ratingComments,'production error')) || ...
                    ~isempty(strfind(pdata.utters.(flds{i2}){i3}.ratingComments,'erroneous prod')))
                    nDscd_prodErr  = nDscd_prodErr + 1;
                end
                continue;
            end
            
            rawDataFN = pdata.utters.(flds{i2}){i3}.rawDataFN;
            [fp1, foo] = fileparts(rawDataFN);
            [foo, fp2] = fileparts(fp1);
            repN = str2num(strrep(fp2, 'rep', ''));

            t_IF2.(flds{i2})(end+1)=pdata.utters.(flds{i2}){i3}.iF2;
            t_UF2.(flds{i2})(end+1)=pdata.utters.(flds{i2}){i3}.uF2;
            t_YF2.(flds{i2})(end+1)=pdata.utters.(flds{i2}){i3}.yF2;
            t_U2F2.(flds{i2})(end+1)=pdata.utters.(flds{i2}){i3}.uF2_you;            
            t_Y2F2.(flds{i2})(end+1)=pdata.utters.(flds{i2}){i3}.yF2_yo1;
            t_U3F2.(flds{i2})(end+1)=pdata.utters.(flds{i2}){i3}.uF2_yo1;
            t_Y3F2.(flds{i2})(end+1)=pdata.utters.(flds{i2}){i3}.yF2_yo2;

            idx_u2 = round((pdata.utters.(flds{i2}){i3}.uayOnset - pdata.utters.(flds{i2}){i3}.iouOnset) / frameDur) + 1;
            idx_y2 = round((pdata.utters.(flds{i2}){i3}.yo1Onset - pdata.utters.(flds{i2}){i3}.iouOnset) / frameDur) + 1;
            idx_u2y2Mid = floor((idx_u2 + idx_y2) / 2);
            idx_u2y2Mid_frac = (idx_u2 + idx_y2) / 2 - idx_u2y2Mid;
            t_U2Y2MidF2.(flds{i2})(end+1) = pdata.utters.(flds{i2}){i3}.traj_F2(idx_u2y2Mid) + ...
                (pdata.utters.(flds{i2}){i3}.traj_F2(idx_u2y2Mid+1) - pdata.utters.(flds{i2}){i3}.traj_F2(idx_u2y2Mid)) * idx_u2y2Mid_frac;

            t_IUInt.(flds{i2})(end+1)=pdata.utters.(flds{i2}){i3}.uTime-pdata.utters.(flds{i2}){i3}.iouOnset;
            t_UYInt.(flds{i2})(end+1)=pdata.utters.(flds{i2}){i3}.youOnset-pdata.utters.(flds{i2}){i3}.uTime;
            t_IYInt.(flds{i2})(end+1)=pdata.utters.(flds{i2}){i3}.youOnset-pdata.utters.(flds{i2}){i3}.iouOnset;
            t_IU2Int.(flds{i2})(end+1)=pdata.utters.(flds{i2}){i3}.uayOnset-pdata.utters.(flds{i2}){i3}.iouOnset;
            t_IY2Int.(flds{i2})(end+1)=pdata.utters.(flds{i2}){i3}.yo1Onset-pdata.utters.(flds{i2}){i3}.iouOnset;
            t_IU3Int.(flds{i2})(end+1)=pdata.utters.(flds{i2}){i3}.yo1End-pdata.utters.(flds{i2}){i3}.iouOnset;
            t_IY3Int.(flds{i2})(end+1)=pdata.utters.(flds{i2}){i3}.yo2Onset-pdata.utters.(flds{i2}){i3}.iouOnset;

            t_I_UYMid_Int.(flds{i2})(end+1)=(pdata.utters.(flds{i2}){i3}.youOnset+pdata.utters.(flds{i2}){i3}.uTime)/2-pdata.utters.(flds{i2}){i3}.iouOnset;
            t_I_YU2Mid_Int.(flds{i2})(end+1)=(pdata.utters.(flds{i2}){i3}.youOnset+pdata.utters.(flds{i2}){i3}.uayOnset)/2-pdata.utters.(flds{i2}){i3}.iouOnset;
            
            this_utter=pdata.utters.(flds{i2}){i3};
            taxis1=this_utter.iouOnset:frameDur:this_utter.iouOnset+frameDur*(length(this_utter.traj_F2)-1);

            f2Seg=this_utter.traj_F2(taxis1>=this_utter.uTime & taxis1<this_utter.youOnset);
            f2Mid=mean([this_utter.uF2,this_utter.yF2]);
            t_I_UYMidF_Int.(flds{i2})(end+1)=get_midF(f2Seg,f2Mid,this_utter.uTime,frameDur,'up')-this_utter.iouOnset;
            
            f2Seg=this_utter.traj_F2(taxis1>=this_utter.youOnset & taxis1<this_utter.uayOnset);
            f2Mid=mean([this_utter.yF2,this_utter.uF2_you]);           
            t_I_YU2MidF_Int.(flds{i2})(end+1)=get_midF(f2Seg,f2Mid,this_utter.youOnset,frameDur,'down')-this_utter.iouOnset;
            
            t_trajF2.(flds{i2}){end+1}=this_utter.traj_F2;
            t_trajF1.(flds{i2}){end+1}=this_utter.traj_F1;
                        
            m_UF2.(fld)(min(find(isnan(m_UF2.(fld)(:, repN)))), repN) = pdata.utters.(flds{i2}){i3}.uF2;
            m_YF2.(fld)(min(find(isnan(m_YF2.(fld)(:, repN)))), repN) = pdata.utters.(flds{i2}){i3}.yF2;
            m_IUInt.(fld)(min(find(isnan(m_IUInt.(fld)(:, repN)))), repN) = ...
                pdata.utters.(flds{i2}){i3}.uTime-pdata.utters.(flds{i2}){i3}.iouOnset;
            m_IYInt.(fld)(min(find(isnan(m_IYInt.(fld)(:, repN)))), repN) = ...
                pdata.utters.(flds{i2}){i3}.youOnset-pdata.utters.(flds{i2}){i3}.iouOnset;
        end
        
        if isequal(fld,'none') && ~isfield(pdata.stage2.normT1_bw_iu_f2{1},'none')
            fld1='noPert';
        else
            fld1=fld;
        end
        
        % === Extract measures from the IOA trajectories (for comparison) ===
        [iuInt, iyInt, iu2Int] = extract_st_meas(pdata.stage2.avg_f2Traj.(fld1)(:, 1), frameDur);
        ioaTraj_tInts.(fld) = [ioaTraj_tInts.(fld); [iuInt, iyInt, iu2Int]];
        % === ~Extract measures from the IOA trajectories (for comparison) ===
        
        
    end
    
    mat_UF2 = [nanmean(m_UF2.(flds{1})); nanmean(m_UF2.(flds{2})); nanmean(m_UF2.(flds{3}))];
    mat_YF2 = [nanmean(m_YF2.(flds{1})); nanmean(m_YF2.(flds{2})); nanmean(m_YF2.(flds{3}))];
    mat_IUInt = [nanmean(m_IUInt.(flds{1})); nanmean(m_IUInt.(flds{2})); nanmean(m_IUInt.(flds{3}))];
    mat_IYInt = [nanmean(m_IYInt.(flds{1})); nanmean(m_IYInt.(flds{2})); nanmean(m_IYInt.(flds{3}))];
    
    % === Extract the F2 values at the landmarks (F2 extrema and their quarter/mid points) ===
    for i3=1:length(config.normT1_bw_iu_ratio)
        if isequal(studyMode,'apstv2')
            normT1_bw_iu_f2{i3}=[normT1_bw_iu_f2{i3};[mean(pdata.stage2.normT1_bw_iu_f2{i3}.noPert),mean(pdata.stage2.normT1_bw_iu_f2{i3}.accel),...
                mean(pdata.stage2.normT1_bw_iu_f2{i3}.decel),mean(pdata.stage2.normT1_bw_iu_f2{i3}.down),mean(pdata.stage2.normT1_bw_iu_f2{i3}.up)]];
        elseif isequal(studyMode,'apstv') || isequal(studyMode,'apstv2_stut_s')
            normT1_bw_iu_f2{i3}=[normT1_bw_iu_f2{i3};[mean(pdata.stage2.normT1_bw_iu_f2{i3}.noPert),...
                mean(pdata.stage2.normT1_bw_iu_f2{i3}.down),mean(pdata.stage2.normT1_bw_iu_f2{i3}.up)]];
        elseif isequal(studyMode,'apstv2t') || isequal(studyMode,'apstv2_stut_t')
            normT1_bw_iu_f2{i3}=[normT1_bw_iu_f2{i3};[mean(pdata.stage2.normT1_bw_iu_f2{i3}.noPert),...
                mean(pdata.stage2.normT1_bw_iu_f2{i3}.accel),mean(pdata.stage2.normT1_bw_iu_f2{i3}.decel)]];
        end
    end
    for i3=1:length(config.normT2_bw_uy_ratio)
        if isequal(studyMode,'apstv2')
            normT2_bw_uy_f2{i3}=[normT2_bw_uy_f2{i3};[mean(pdata.stage2.normT2_bw_uy_f2{i3}.noPert),mean(pdata.stage2.normT2_bw_uy_f2{i3}.accel),...
                mean(pdata.stage2.normT2_bw_uy_f2{i3}.decel),mean(pdata.stage2.normT2_bw_uy_f2{i3}.down),mean(pdata.stage2.normT2_bw_uy_f2{i3}.up)]];
        elseif isequal(studyMode,'apstv') || isequal(studyMode,'apstv2_stut_s')
            normT2_bw_uy_f2{i3}=[normT2_bw_uy_f2{i3};[mean(pdata.stage2.normT2_bw_uy_f2{i3}.noPert),...
                mean(pdata.stage2.normT2_bw_uy_f2{i3}.down),mean(pdata.stage2.normT2_bw_uy_f2{i3}.up)]];
        elseif isequal(studyMode,'apstv2t') || isequal(studyMode,'apstv2_stut_t')
            normT2_bw_uy_f2{i3}=[normT2_bw_uy_f2{i3};[mean(pdata.stage2.normT2_bw_uy_f2{i3}.noPert),...
                mean(pdata.stage2.normT2_bw_uy_f2{i3}.accel),mean(pdata.stage2.normT2_bw_uy_f2{i3}.decel)]];
        end
    end
    for i3=1:length(config.normT3_bw_yu_ratio)
        if isequal(studyMode,'apstv2')
            normT3_bw_yu_f2{i3}=[normT3_bw_yu_f2{i3};[mean(pdata.stage2.normT3_bw_yu_f2{i3}.noPert),mean(pdata.stage2.normT3_bw_yu_f2{i3}.accel),...
                mean(pdata.stage2.normT3_bw_yu_f2{i3}.decel),mean(pdata.stage2.normT3_bw_yu_f2{i3}.down),mean(pdata.stage2.normT3_bw_yu_f2{i3}.up)]];
        elseif isequal(studyMode,'apstv') || isequal(studyMode,'apstv2_stut_s')
            normT3_bw_yu_f2{i3}=[normT3_bw_yu_f2{i3};[mean(pdata.stage2.normT3_bw_yu_f2{i3}.noPert),...
                mean(pdata.stage2.normT3_bw_yu_f2{i3}.down),mean(pdata.stage2.normT3_bw_yu_f2{i3}.up)]];
        elseif isequal(studyMode,'apstv2t') || isequal(studyMode,'apstv2_stut_t')
            normT3_bw_yu_f2{i3}=[normT3_bw_yu_f2{i3};[mean(pdata.stage2.normT3_bw_yu_f2{i3}.noPert),...
                mean(pdata.stage2.normT3_bw_yu_f2{i3}.accel),mean(pdata.stage2.normT3_bw_yu_f2{i3}.decel)]];
        end
    end
    % === ~Extract the F2 values at the landmarks (F2 extrema and their quarter/mid points) ===
    
    % === Extract differences in avg. F1/F2 trajectories b/w the pert conditions and the noPert condition ===
    flds0=fields(pdata.stage2.avg_f1Traj_POA);
    idxBL=fsic(flds0,'noPert');
    if isempty(idxBL)
        idxBL=fsic(flds0,'none');
    end
    for k1=1:length(flds0)
        if isequal(flds{k1},'none')
            fld1='noPert';
        else
            fld1=flds{k1};
        end
        trajF2_IOA_NN.(flds{k1}){end+1}=pdata.stage2.avg_f2Traj.(fld1)(:,1);
        
        t_mat=pdata.stage2.fullTNormF2.(fld1);
        idx_ok=find(~isnan(sum(t_mat')));
%         fe_trajF2_FTN.(flds{k1})=[fe_trajF2_FTN.(flds{k1}); t_mat(idx_ok,:)];
        
        if isequal(FTN_mode, 'linear')
            avgTrajF2_FTN.(flds{k1}){end+1}=pdata.stage2.avg_fullTNormF2.(fld1)(:,1);            
        elseif isequal(FTN_mode, 'nearest')
            avgTrajF2_FTN.(flds{k1}){end+1}=pdata.stage2.avg_fullTNormF2_nn.(fld1)(:,1);
        elseif isequal(FTN_mode, 'spline')
            avgTrajF2_FTN.(flds{k1}){end+1}=pdata.stage2.avg_fullTNormF2_cs.(fld1)(:,1);
        end
        avgTrajF2_FTN_FN.(flds{k1}){end+1}=pdata.stage2.avg_fullTNormF2_FN.(fld1)(:,1);
        
%         fe_trajF2_FTN.(flds{k1})=[fe_trajF2_FTN.(flds{k1}); pdata.stage2.avg_fullTNormF2.(fld1)(:,1)'];
        
        if k1==idxBL
            continue;            
        end
        
        avgTrajF2_pert_FTN_FN.(flds0{k1})=[avgTrajF2_pert_FTN_FN.(flds0{k1}); pdata.stage2.avg_fullTNormF2_pert_FN.(flds0{k1})(:,1)'];
        
        len=min([length(pdata.stage2.avg_f1Traj_POA.(flds0{k1})(:,1)),length(pdata.stage2.avg_f1Traj_POA.(flds0{idxBL})(:,1))]);
        chgTrajF1_POA.(flds0{k1}){end+1}=pdata.stage2.avg_f1Traj_POA.(flds0{k1})(1:len,1)-pdata.stage2.avg_f1Traj_POA.(flds0{idxBL})(1:len,1);
        len=min([length(pdata.stage2.avg_f2Traj_POA.(flds0{k1})(:,1)),length(pdata.stage2.avg_f2Traj_POA.(flds0{idxBL})(:,1))]);
        chgTrajF2_POA.(flds0{k1}){end+1}=pdata.stage2.avg_f2Traj_POA.(flds0{k1})(1:len,1)-pdata.stage2.avg_f2Traj_POA.(flds0{idxBL})(1:len,1);
        
        len=min([length(pdata.stage2.avg_f1Traj_POA_TN.(flds0{k1})(:,1)),length(pdata.stage2.avg_f1Traj_POA_TN.(flds0{idxBL})(:,1))]);
        chgTrajF1_POA_TN.(flds0{k1}){end+1}=pdata.stage2.avg_f1Traj_POA_TN.(flds0{k1})(1:len,1)-pdata.stage2.avg_f1Traj_POA_TN.(flds0{idxBL})(1:len,1);
        len=min([length(pdata.stage2.avg_f2Traj_POA_TN.(flds0{k1})(:,1)),length(pdata.stage2.avg_f2Traj_POA_TN.(flds0{idxBL})(:,1))]);
        chgTrajF2_POA_TN.(flds0{k1}){end+1}=pdata.stage2.avg_f2Traj_POA_TN.(flds0{k1})(1:len,1)-pdata.stage2.avg_f2Traj_POA_TN.(flds0{idxBL})(1:len,1);                
        
        len=min([length(pdata.stage2.avg_f1Traj_POA_TN_FN.(flds0{k1})(:,1)),length(pdata.stage2.avg_f1Traj_POA_TN_FN.(flds0{idxBL})(:,1))]);
        chgTrajF1_POA_TN_FN.(flds0{k1}){end+1}=pdata.stage2.avg_f1Traj_POA_TN_FN.(flds0{k1})(1:len,1)-pdata.stage2.avg_f1Traj_POA_TN_FN.(flds0{idxBL})(1:len,1);
        len=min([length(pdata.stage2.avg_f2Traj_POA_TN_FN.(flds0{k1})(:,1)),length(pdata.stage2.avg_f2Traj_POA_TN_FN.(flds0{idxBL})(:,1))]);
        chgTrajF2_POA_TN_FN.(flds0{k1}){end+1}=pdata.stage2.avg_f2Traj_POA_TN_FN.(flds0{k1})(1:len,1)-pdata.stage2.avg_f2Traj_POA_TN_FN.(flds0{idxBL})(1:len,1);
        
        len=min([length(pdata.stage2.avg_f1Traj.(flds0{k1})(:,1)),length(pdata.stage2.avg_f1Traj.(flds0{idxBL})(:,1))]);
        chgTrajF1_IOA.(flds0{k1}){end+1}=pdata.stage2.avg_f1Traj.(flds0{k1})(1:len,1)-pdata.stage2.avg_f1Traj.(flds0{idxBL})(1:len,1);
        len=min([length(pdata.stage2.avg_f2Traj.(flds0{k1})(:,1)),length(pdata.stage2.avg_f2Traj.(flds0{idxBL})(:,1))]);
        chgTrajF2_IOA.(flds0{k1}){end+1}=pdata.stage2.avg_f2Traj.(flds0{k1})(1:len,1)-pdata.stage2.avg_f2Traj.(flds0{idxBL})(1:len,1);
        semChgTrajF2_IOA.(flds0{k1}){end+1}=pdata.stage2.avg_f2Traj.(flds0{k1})(1:len,2) ./ sqrt(pdata.stage2.avg_f2Traj.(flds0{k1})(1:len,3));
        
        if isequal(FTN_mode, 'linear')
            semChgTrajF2_FTN.(flds0{k1}){end+1}=pdata.stage2.avg_fullTNormF2.(flds0{k1})(: ,2) ./ sqrt(pdata.stage2.avg_fullTNormF2.(flds0{k1})(: ,3));
        elseif isequal(FTN_mode, 'nearest')
            semChgTrajF2_FTN.(flds0{k1}){end+1}=pdata.stage2.avg_fullTNormF2_nn.(flds0{k1})(: ,2) ./ sqrt(pdata.stage2.avg_fullTNormF2_nn.(flds0{k1})(: ,3));
        elseif isequal(FTN_mode, 'spline')
            semChgTrajF2_FTN.(flds0{k1}){end+1}=pdata.stage2.avg_fullTNormF2_cs.(flds0{k1})(: ,2) ./ sqrt(pdata.stage2.avg_fullTNormF2_cs.(flds0{k1})(: ,3));
        end
        
        if isequal(FTN_mode, 'linear')
            chgTrajF2_FTN.(flds0{k1}){end+1}=pdata.stage2.avg_fullTNormF2.(flds0{k1})(:,1)-pdata.stage2.avg_fullTNormF2.(flds0{idxBL})(:,1);
        elseif isequal(FTN_mode, 'nearest')
            chgTrajF2_FTN.(flds0{k1}){end+1}=pdata.stage2.avg_fullTNormF2_nn.(flds0{k1})(:,1)-pdata.stage2.avg_fullTNormF2_nn.(flds0{idxBL})(:,1);
        elseif isequal(FTN_mode, 'spline')
            chgTrajF2_FTN.(flds0{k1}){end+1}=pdata.stage2.avg_fullTNormF2_cs.(flds0{k1})(:,1)-pdata.stage2.avg_fullTNormF2_cs.(flds0{idxBL})(:,1);
        end
        
        frameDur_FTN=pdata.stage2.fullTNormF2_tAxis(2)-pdata.stage2.fullTNormF2_tAxis(1);
        
        chgTrajF2_FTN_TA.(['aft_',flds0{k1}]){end+1}=pdata.stage2.avg_fullTNormF2_AT.(['aft_',flds0{k1}])(:,1)-pdata.stage2.avg_fullTNormF2_AT.aft_noPert(:,1);    % Adaptation analysis
        
        
        if k1 == 2
            len1 = min([size(pdata.stage2.avg_f2Traj_POA.(flds_cntrst{1}), 1), size(pdata.stage2.avg_f2Traj_POA.(flds_cntrst{2}), 1)]);
            chgTrajF2_POA.contrast{end + 1} = pdata.stage2.avg_f2Traj_POA.(flds_cntrst{1})(1:len1,1) - pdata.stage2.avg_f2Traj_POA.(flds_cntrst{2})(1:len1,1);
            
            len1 = min([size(pdata.stage2.avg_f2Traj.(flds_cntrst{1}), 1), size(pdata.stage2.avg_f2Traj.(flds_cntrst{2}), 1)]);
            chgTrajF2_IOA.contrast{end + 1} = pdata.stage2.avg_f2Traj.(flds_cntrst{1})(1:len1,1) - pdata.stage2.avg_f2Traj.(flds_cntrst{2})(1:len1,1);
            
            if isequal(FTN_mode, 'linear')
                chgTrajF2_FTN.contrast{end + 1} = pdata.stage2.avg_fullTNormF2.(flds_cntrst{1})(:, 1) - pdata.stage2.avg_fullTNormF2.(flds_cntrst{2})(:, 1);
            elseif isequal(FTN_mode, 'nearest')
                chgTrajF2_FTN.contrast{end + 1} = pdata.stage2.avg_fullTNormF2_nn.(flds_cntrst{1})(:, 1) - pdata.stage2.avg_fullTNormF2_nn.(flds_cntrst{2})(:, 1);
            elseif isequal(FTN_mode, 'spline')
                chgTrajF2_FTN.contrast{end + 1} = pdata.stage2.avg_fullTNormF2_cs.(flds_cntrst{1})(:, 1) - pdata.stage2.avg_fullTNormF2_cs.(flds_cntrst{2})(:, 1);
            end
        end
        
        t_mat=pdata.stage2.fullTNormF2.(flds0{k1});
        idx_ok=find(~isnan(sum(t_mat')));
        
        fe_chgTrajF2_FTN.(flds0{k1})=[fe_chgTrajF2_FTN.(flds0{k1}); t_mat(idx_ok,:)-repmat(pdata.stage2.avg_fullTNormF2.(flds0{idxBL})(:,1)',numel(idx_ok),1)];
        
        if ~(isequal(flds0{k1},'none') || isequal(flds0{k1},'noPert'))
            pertShiftF2_FTN.(flds0{k1}){end+1}=pdata.stage2.avg_fullTNormPertShiftF2.(flds0{k1})(:,1);
        end
        
        
        
    end
    % === ~Extract differences in avg. F1/F2 trajectories b/w the pert conditions and the noPert condition ===
    
    % === Extract RMS trajectory, if availabel === 
    if isfield(pdata.stage2, 'avg_rmsTraj')
        for k0 = 1 : numel(flds)
            if isequal(flds{k0}, 'none') && isfield(pdata.stage2.avg_rmsTraj, 'noPert')
                flds00 = 'noPert';
            else
                flds00 = flds{k0};
            end
            trajRMS_IOA.(flds{k0}){end + 1} = pdata.stage2.avg_rmsTraj.(flds00)(:, 1);
        end
%         trajRMS_IOA.
    end
    % === ~Extract RMS trajectory, if availabel === 
    
    if isequal(studyMode,'apstv2') 
        a_IUInt=[a_IUInt;[mean(t_IUInt.none),mean(t_IUInt.accel),mean(t_IUInt.decel),mean(t_IUInt.down),mean(t_IUInt.up)]];
        a_IYInt=[a_IYInt;[mean(t_IYInt.none),mean(t_IYInt.accel),mean(t_IYInt.decel),mean(t_IYInt.down),mean(t_IYInt.up)]];
        
        a_iF2=[a_iF2;[mean(t_IF2.none),mean(t_IF2.accel),mean(t_IF2.decel),mean(t_IF2.down),mean(t_IF2.up)]];
        a_uF2=[a_uF2;[mean(t_UF2.none),mean(t_UF2.accel),mean(t_UF2.decel),mean(t_UF2.down),mean(t_UF2.up)]];
        a_yF2=[a_yF2;[mean(t_YF2.none),mean(t_YF2.accel),mean(t_YF2.decel),mean(t_YF2.down),mean(t_YF2.up)]];
        a_u2F2=[a_u2F2;[mean(t_U2F2.none),mean(t_U2F2.accel),mean(t_U2F2.decel),mean(t_U2F2.down),mean(t_U2F2.up)]];
        a_y2F2=[a_y2F2;[mean(t_Y2F2.none),mean(t_Y2F2.accel),mean(t_Y2F2.decel),mean(t_Y2F2.down),mean(t_Y2F2.up)]];
        a_u2y2MidF2=[a_u2F2;[mean(t_U2Y2MidF2.none),mean(t_U2Y2MidF2.accel),mean(t_U2Y2MidF2.decel),mean(t_U2Y2MidF2.down),mean(t_U2Y2MidF2.up)]];
        a_u3F2=[a_u3F2;[mean(t_U3F2.none),mean(t_U3F2.accel),mean(t_U3F2.decel),mean(t_U3F2.down),mean(t_U3F2.up)]];
        a_y3F2=[a_y3F2;[mean(t_Y3F2.none),mean(t_Y3F2.accel),mean(t_Y3F2.decel),mean(t_Y3F2.down),mean(t_Y3F2.up)]];
        
        a_IU2Int=[a_IU2Int;[mean(t_IU2Int.none),mean(t_IU2Int.accel),mean(t_IU2Int.decel),mean(t_IU2Int.down),mean(t_IU2Int.up)]];
        a_IY2Int=[a_IY2Int;[mean(t_IY2Int.none),mean(t_IY2Int.accel),mean(t_IY2Int.decel),mean(t_IY2Int.down),mean(t_IY2Int.up)]];
        a_IU3Int=[a_IU3Int;[mean(t_IU3Int.none),mean(t_IU3Int.accel),mean(t_IU3Int.decel),mean(t_IU3Int.down),mean(t_IU3Int.up)]];
        a_IY3Int=[a_IY3Int;[mean(t_IY3Int.none),mean(t_IY3Int.accel),mean(t_IY3Int.decel),mean(t_IY3Int.down),mean(t_IY3Int.up)]];
        
        a_Y1Y2Int=[a_Y1Y2Int;a_IY2Int(end,:)-a_IYInt(end,:)];
        a_Y2Y3Int=[a_Y2Y3Int;a_IY3Int(end,:)-a_IY2Int(end,:)];
    
        a_I_UYMid_Int=[a_I_UYMid_Int;[mean(t_I_UYMid_Int.none),mean(t_I_UYMid_Int.accel),mean(t_I_UYMid_Int.decel),mean(t_I_UYMid_Int.down),mean(t_I_UYMid_Int.up)]];
        a_I_YU2Mid_Int=[a_I_YU2Mid_Int;[mean(t_I_YU2Mid_Int.none),mean(t_I_YU2Mid_Int.accel),mean(t_I_YU2Mid_Int.decel),mean(t_I_YU2Mid_Int.down),mean(t_I_YU2Mid_Int.up)]];
    
        a_I_UYMidF_Int=[a_I_UYMidF_Int;[mean(t_I_UYMidF_Int.none),mean(t_I_UYMidF_Int.accel),mean(t_I_UYMidF_Int.decel),mean(t_I_UYMidF_Int.down),mean(t_I_UYMidF_Int.up)]];
        a_I_YU2MidF_Int=[a_I_YU2MidF_Int;[mean(t_I_YU2MidF_Int.none),mean(t_I_YU2MidF_Int.accel),mean(t_I_YU2MidF_Int.decel),mean(t_I_YU2MidF_Int.down),mean(t_I_YU2MidF_Int.up)]];
    elseif isequal(studyMode,'apstv') || isequal(studyMode,'apstv2_stut_s')
        a_IUInt=[a_IUInt;[mean(t_IUInt.none),mean(t_IUInt.down),mean(t_IUInt.up)]];
        a_IYInt=[a_IYInt;[mean(t_IYInt.none),mean(t_IYInt.down),mean(t_IYInt.up)]];
        
        a_iF2=[a_iF2;[mean(t_IF2.none),mean(t_IF2.down),mean(t_IF2.up)]];
        a_uF2=[a_uF2;[mean(t_UF2.none),mean(t_UF2.down),mean(t_UF2.up)]];
        a_yF2=[a_yF2;[mean(t_YF2.none),mean(t_YF2.down),mean(t_YF2.up)]];
        a_u2F2=[a_u2F2;[mean(t_U2F2.none),mean(t_U2F2.down),mean(t_U2F2.up)]];
        a_y2F2=[a_y2F2;[mean(t_Y2F2.none),mean(t_Y2F2.down),mean(t_Y2F2.up)]];
        a_u2y2MidF2=[a_u2F2;[mean(t_U2Y2MidF2.none),mean(t_U2Y2MidF2.down),mean(t_U2Y2MidF2.up)]];
        a_u3F2=[a_u3F2;[mean(t_U3F2.none),mean(t_U3F2.down),mean(t_U3F2.up)]];
        a_y3F2=[a_y3F2;[mean(t_Y3F2.none),mean(t_Y3F2.down),mean(t_Y3F2.up)]];
        
        a_IU2Int=[a_IU2Int;[mean(t_IU2Int.none),mean(t_IU2Int.down),mean(t_IU2Int.up)]];
        a_IY2Int=[a_IY2Int;[mean(t_IY2Int.none),mean(t_IY2Int.down),mean(t_IY2Int.up)]];
        a_IU3Int=[a_IU3Int;[mean(t_IU3Int.none),mean(t_IU3Int.down),mean(t_IU3Int.up)]];
        a_IY3Int=[a_IY3Int;[mean(t_IY3Int.none),mean(t_IY3Int.down),mean(t_IY3Int.up)]];
        
        a_Y1Y2Int=[a_Y1Y2Int;a_IY2Int(end,:)-a_IYInt(end,:)];
        a_Y2Y3Int=[a_Y2Y3Int;a_IY3Int(end,:)-a_IY2Int(end,:)];

        a_I_UYMid_Int=[a_I_UYMid_Int;[mean(t_I_UYMid_Int.none),mean(t_I_UYMid_Int.down),mean(t_I_UYMid_Int.up)]];
        a_I_YU2Mid_Int=[a_I_YU2Mid_Int;[mean(t_I_YU2Mid_Int.none),mean(t_I_YU2Mid_Int.down),mean(t_I_YU2Mid_Int.up)]];
    
        a_I_UYMidF_Int=[a_I_UYMidF_Int;[nanmean(t_I_UYMidF_Int.none),nanmean(t_I_UYMidF_Int.down),nanmean(t_I_UYMidF_Int.up)]];
        a_I_YU2MidF_Int=[a_I_YU2MidF_Int;[nanmean(t_I_YU2MidF_Int.none),nanmean(t_I_YU2MidF_Int.down),nanmean(t_I_YU2MidF_Int.up)]];
    elseif isequal(studyMode,'apstv2t') || isequal(studyMode,'apstv2_stut_t')
        a_IUInt=[a_IUInt;[mean(t_IUInt.none),mean(t_IUInt.accel),mean(t_IUInt.decel)]];
        a_IYInt=[a_IYInt;[mean(t_IYInt.none),mean(t_IYInt.accel),mean(t_IYInt.decel)]];
        
        a_iF2=[a_iF2;[mean(t_IF2.none),mean(t_IF2.accel),mean(t_IF2.decel)]];
        a_uF2=[a_uF2;[mean(t_UF2.none),mean(t_UF2.accel),mean(t_UF2.decel)]];
        a_yF2=[a_yF2;[mean(t_YF2.none),mean(t_YF2.accel),mean(t_YF2.decel)]];
        a_u2F2=[a_u2F2;[mean(t_U2F2.none),mean(t_U2F2.accel),mean(t_U2F2.decel)]];
        a_y2F2=[a_y2F2;[mean(t_Y2F2.none),mean(t_Y2F2.accel),mean(t_Y2F2.decel)]];
        a_u2y2MidF2=[a_u2F2;[mean(t_U2Y2MidF2.none),mean(t_U2Y2MidF2.accel),mean(t_U2Y2MidF2.decel)]];
        a_u3F2=[a_u3F2;[mean(t_U3F2.none),mean(t_U3F2.accel),mean(t_U3F2.decel)]];
        a_y3F2=[a_y3F2;[mean(t_Y3F2.none),mean(t_Y3F2.accel),mean(t_Y3F2.decel)]];
        
        a_IU2Int=[a_IU2Int;[mean(t_IU2Int.none),mean(t_IU2Int.accel),mean(t_IU2Int.decel)]];
        a_IY2Int=[a_IY2Int;[mean(t_IY2Int.none),mean(t_IY2Int.accel),mean(t_IY2Int.decel)]];
        a_IU3Int=[a_IU3Int;[mean(t_IU3Int.none),mean(t_IU3Int.accel),mean(t_IU3Int.decel)]];
        a_IY3Int=[a_IY3Int;[mean(t_IY3Int.none),mean(t_IY3Int.accel),mean(t_IY3Int.decel)]];
        
        a_Y1Y2Int=[a_Y1Y2Int;a_IY2Int(end,:)-a_IYInt(end,:)];
        a_Y2Y3Int=[a_Y2Y3Int;a_IY3Int(end,:)-a_IY2Int(end,:)];

        a_I_UYMid_Int=[a_I_UYMid_Int;[mean(t_I_UYMid_Int.none),mean(t_I_UYMid_Int.accel),mean(t_I_UYMid_Int.decel)]];
        a_I_YU2Mid_Int=[a_I_YU2Mid_Int;[mean(t_I_YU2Mid_Int.none),mean(t_I_YU2Mid_Int.accel),mean(t_I_YU2Mid_Int.decel)]];
    
        a_I_UYMidF_Int=[a_I_UYMidF_Int;[mean(t_I_UYMidF_Int.none),mean(t_I_UYMidF_Int.accel),mean(t_I_UYMidF_Int.decel)]];
        a_I_YU2MidF_Int=[a_I_YU2MidF_Int;[mean(t_I_YU2MidF_Int.none),mean(t_I_YU2MidF_Int.accel),mean(t_I_YU2MidF_Int.decel)]];
    end
    
    % === Trend conformation analysis ===
    if isequal(studyMode, 'apstv')
        [cfmVect.iuInt(i1, :), ssMult.iuInt(i1, :), cfmVect_byRep.iuInt(i1, :)] = getCfmVect(t_IUInt, 'iuInt', studyMode, mat_IUInt, flds);
        [cfmVect.uF2(i1, :), ssMult.uF2(i1, :), cfmVect_byRep.uF2(i1, :)] = getCfmVect(t_UF2, 'uF2', studyMode, mat_UF2, flds);
%         cfmVect.uyMidF2(i1, :) = getCfmVect(t_UYMidF2, 'y2u2MidF2', studyMode);
        [cfmVect.yF2(i1, :), ssMult.yF2(i1, :), cfmVect_byRep.yF2(i1, :)] = getCfmVect(t_YF2, 'yF2', studyMode, mat_YF2, flds);
%         cfmVect.y2u2MidF2(i1, :) = getCfmVect(t_U2Y2MidF2, 'y2u2MidF2', studyMode);
    else
        [cfmVect.iuInt(i1, :), ssMult.iuInt(i1, :), cfmVect_byRep.iuInt(i1, :)] = getCfmVect(t_IUInt, 'iuInt', studyMode, mat_IUInt, flds);
        [cfmVect.iyInt(i1, :), ssMult.iyInt(i1, :), cfmVect_byRep.iyInt(i1, :)] = getCfmVect(t_IYInt, 'iuInt', studyMode, mat_IYInt, flds);
    end
    % === ~Trend conformation analysis ===
    
    if ~isequal(studyMode,'apstv') || isequal(studyMode,'apstv2_stut_s')
        a_tShift_iuMidF.accel(end+1)=nanmean(pdata.stage2.timeShifts.iuMidF.accel);
        a_tShift_iuMidF.decel(end+1)=nanmean(pdata.stage2.timeShifts.iuMidF.decel);    
    end
%     for i3=1:length(config.normT1_bw_iu_ratio)
%         change_bw_iu_f2{i3}.accel=[change_bw_iu_f2{i3}.accel,normT1_bw_iu_f2{i3}.accel(end)-normT1_bw_iu_f2{i3}.none(end)];
%         change_bw_iu_f2{i3}.decel=[change_bw_iu_f2{i3}.decel,normT1_bw_iu_f2{i3}.decel(end)-normT1_bw_iu_f2{i3}.none(end)];
%     end
%     for i3=1:length(config.normT1_bw_uy_ratio)
%         change_bw_uy_f2{i3}.accel=[change_bw_uy_f2{i3}.accel,normT2_bw_uy_f2{i3}.accel(end)-normT2_bw_uy_f2{i3}.none(end)];
%         change_bw_uy_f2{i3}.decel=[change_bw_uy_f2{i3}.decel,normT2_bw_uy_f2{i3}.decel(end)-normT2_bw_uy_f2{i3}.none(end)];
%     end
%     for i3=1:length(config.normT3_bw_yu_ratio)
%         change_bw_yu_f2{i3}.accel=[change_bw_yu_f2{i3}.accel,normT3_bw_yu_f2{i3}.accel(end)-normT3_bw_yu_f2{i3}.none(end)];
%         change_bw_yu_f2{i3}.decel=[change_bw_yu_f2{i3}.decel,normT3_bw_yu_f2{i3}.decel(end)-normT3_bw_yu_f2{i3}.none(end)];
%     end
    
    flds=fields(colors);
    
    for i2=1:length(flds)
        fld=flds{i2};
        
        if isequal(studyMode,'apstv')  || isequal(studyMode,'apstv2_stut_s')
            if isequal(fld,'accel') || isequal(fld,'decel')
                continue;
            end
        elseif isequal(studyMode,'apstv2t') || isequal(studyMode,'apstv2_stut_t')
            if isequal(fld,'down') || isequal(fld,'up')
                continue;
            end
        end
        
        if isequal(fld,'none')
            if ~isfield(pdata.stage2.avg_f2Traj_POA,'none')
                fld1='noPert';
            else
                fld1=fld;
            end            
            
            % --- POA ---
            t_traj_f2=pdata.stage2.avg_f2Traj_POA.(fld1)(:,1);
            diff_t_traj_f2=diff(t_traj_f2);
            idx_turn=find(diff_t_traj_f2(1:end-1).*diff_t_traj_f2(2:end)<0);
            f2LB=t_traj_f2(idx_turn(1)+1);
            f2UB=t_traj_f2(idx_turn(2)+1);
            
            d_t_traj_f2=diff(t_traj_f2);
            idx_peak_1=find(d_t_traj_f2(1:end-1)>0 & d_t_traj_f2(2:end)<0);
            
            idxNormInt=idx_peak_1(1)+1;
            
            t_traj_f1=pdata.stage2.avg_f1Traj_POA.(fld1)(:,1);
            diff_t_traj_f1=diff(t_traj_f1);
            idx_maxima=find(diff_t_traj_f1(1:end-1).*diff_t_traj_f1(2:end)<0 & diff_t_traj_f1(1:end-1)>0);
            idx_maxima=idx_maxima(1);            
            idx_minima=find(diff_t_traj_f1(1:end-1).*diff_t_traj_f1(2:end)<0 & diff_t_traj_f1(1:end-1)<0);
            idx_minima=idx_minima(idx_minima>idx_maxima);            
            f1LB=t_traj_f1(idx_minima(1)+1);
            f1UB=t_traj_f1(idx_maxima+1);
            
            % --- IOA ---
            t_traj_f2=pdata.stage2.avg_f2Traj.(fld1)(:,1);
            diff_t_traj_f2=diff(t_traj_f2);
            idx_turn=find(diff_t_traj_f2(1:end-1).*diff_t_traj_f2(2:end)<0);
            f2LB_IOA=t_traj_f2(idx_turn(1)+1);
            f2UB_IOA=t_traj_f2(idx_turn(2)+1);
            
            d_t_traj_f2=diff(t_traj_f2);
            idx_peak_1=find(d_t_traj_f2(1:end-1)>0 & d_t_traj_f2(2:end)<0);
            
            idxNormInt_IOA=idx_peak_1(1)+1;
            
            t_traj_f1=pdata.stage2.avg_f1Traj_POA.(fld1)(:,1);
            diff_t_traj_f1=diff(t_traj_f1);
            idx_maxima=find(diff_t_traj_f1(1:end-1).*diff_t_traj_f1(2:end)<0 & diff_t_traj_f1(1:end-1)>0);
            idx_maxima=idx_maxima(1);            
            idx_minima=find(diff_t_traj_f1(1:end-1).*diff_t_traj_f1(2:end)<0 & diff_t_traj_f1(1:end-1)<0);
            idx_minima=idx_minima(idx_minima>idx_maxima);            
            f1LB_IOA=t_traj_f1(idx_minima(1)+1);
            f1UB_IOA=t_traj_f1(idx_maxima+1);
        end
        if isequal(fld,'none') && ~isfield(pdata.stage2.avg_f2Traj_POA,'none')
            fld1='noPert';
        else
            fld1=fld;
        end
        % --- POA_TN ---
        for i3=1:length(pdata.stage2.trajF2_POA.(fld1))            
            t_ind_traj_f2=pdata.stage2.trajF2_POA.(fld1){i3};
            t_ind_traj_f1=pdata.stage2.trajF1_POA.(fld1){i3};

            fNormF2=(t_ind_traj_f2-f2LB)/(f2UB-f2LB);
            fNormF1=(t_ind_traj_f1-f1LB)/(f1UB-f1LB);
            
            trajF2_POA.(fld){end+1}=fNormF2;
            velF2_POA.(fld){end+1}=diff(trajF2_POA.(fld){end});
            trajF1_POA.(fld){end+1}=fNormF1;
            velF1_POA.(fld){end+1}=diff(trajF1_POA.(fld){end});
            
            oldTAxis=1:length(fNormF2);
            newTAxis=linspace(1,idxNormInt,400);
            dt=newTAxis(2)-newTAxis(1);
            newTAxis=[newTAxis,newTAxis(end)+dt:dt:oldTAxis(end)];
            trajF2_POA_TN.(fld){end+1}=interp1(oldTAxis,fNormF2,newTAxis);
            trajF1_POA_TN.(fld){end+1}=interp1(oldTAxis,fNormF1,newTAxis);
        end
        % --- IOA_TN ---
        for i3=1:length(t_trajF2.(fld))
            t_ind_traj_f2=t_trajF2.(fld){i3};
            t_ind_traj_f1=t_trajF1.(fld){i3};

            fNormF2=(t_ind_traj_f2-f2LB_IOA)/(f2UB_IOA-f2LB_IOA);
            fNormF1=(t_ind_traj_f1-f1LB_IOA)/(f1UB_IOA-f1LB_IOA);
            
            if isequal(opt_avg,'all_trajs')
                trajF2_IOA.(fld){end+1}=fNormF2;
                trajF1_IOA.(fld){end+1}=fNormF1;
            end
            
            
            velF2_IOA.(fld){end+1}=diff(fNormF2);
            velF1_IOA.(fld){end+1}=diff(fNormF1);
            
            oldTAxis=1:length(fNormF2);
            newTAxis=linspace(1,idxNormInt_IOA,400);
            dt=newTAxis(2)-newTAxis(1);
            newTAxis=[newTAxis,newTAxis(end)+dt:dt:oldTAxis(end)];
            trajF2_IOA_TN.(fld){end+1}=interp1(oldTAxis,fNormF2,newTAxis);
            trajF1_IOA_TN.(fld){end+1}=interp1(oldTAxis,fNormF1,newTAxis);
        end
        
        if isequal(opt_avg,'ind_avg_trajs')
            trajF2_IOA.(fld){end+1}=(pdata.stage2.avg_f2Traj.(fld1)-f2LB_IOA)/(f2UB_IOA-f2LB_IOA);
            trajF1_IOA.(fld){end+1}=(pdata.stage2.avg_f1Traj.(fld1)-f1LB_IOA)/(f1UB_IOA-f1LB_IOA);
        end
        if ~isequal(fld,'none')
            trajF2_pert_IOA.(fld){end+1}=pdata.stage2.avg_f2Traj_pert.(fld);
            trajF2_pert_IOA_FN.(fld){end+1}=(pdata.stage2.avg_f2Traj_pert.(fld)-f2LB_IOA)/(f2UB_IOA-f2LB_IOA);
        end
    end
    
     % === Extract differences in avg. F1/F2 (TN) trajectories b/w the pert conditions and the noPert condition ===
    flds0=fields(pdata.stage2.avg_f2Traj_IOA_TN);
    idxBL=fsic(flds0,'noPert');
    if isempty(idxBL)
        idxBL=fsic(flds0,'none');
    end
    for k1=1:length(flds0)
        if k1==idxBL
            continue;            
        end
        fld=flds0{k1};

        len=min([length(pdata.stage2.avg_f2Traj_IOA_TN.(fld)(:,1)),length(pdata.stage2.avg_f2Traj_IOA_TN.(flds0{idxBL})(:,1))]);
        chgTrajF2_IOA_TN.(flds0{k1}){end+1}=pdata.stage2.avg_f2Traj_IOA_TN.(fld)(1:len,1)-...
            pdata.stage2.avg_f2Traj_IOA_TN.(flds0{idxBL})(1:len,1);
        len=min([length(pdata.stage2.avg_f1Traj_IOA_TN.(fld)(:,1)),length(pdata.stage2.avg_f1Traj_IOA_TN.(flds0{idxBL})(:,1))]);
        chgTrajF1_IOA_TN.(flds0{k1}){end+1}=pdata.stage2.avg_f1Traj_IOA_TN.(fld)(1:len,1)-...
            pdata.stage2.avg_f1Traj_IOA_TN.(flds0{idxBL})(1:len,1);
        
        len=min([length(pdata.stage2.avg_f2Traj_IOA_FN.(fld)(:,1)),length(pdata.stage2.avg_f2Traj_IOA_FN.(flds0{idxBL})(:,1))]);
        chgTrajF2_IOA_FN.(flds0{k1}){end+1}=pdata.stage2.avg_f2Traj_IOA_FN.(fld)(1:len,1)-...
            pdata.stage2.avg_f2Traj_IOA_FN.(flds0{idxBL})(1:len,1);
        len=min([length(pdata.stage2.avg_f1Traj_IOA_FN.(fld)(:,1)),length(pdata.stage2.avg_f1Traj_IOA_FN.(flds0{idxBL})(:,1))]);
        chgTrajF1_IOA_FN.(flds0{k1}){end+1}=pdata.stage2.avg_f1Traj_IOA_FN.(fld)(1:len,1)-...
            pdata.stage2.avg_f1Traj_IOA_FN.(flds0{idxBL})(1:len,1);
        
        len=min([length(pdata.stage2.avg_f2Traj_IOA_TN_FN.(fld)(:,1)),length(pdata.stage2.avg_f2Traj_IOA_TN_FN.(flds0{idxBL})(:,1))]);
        chgTrajF2_IOA_TN_FN.(flds0{k1}){end+1}=pdata.stage2.avg_f2Traj_IOA_TN_FN.(fld)(1:len,1)-...
            pdata.stage2.avg_f2Traj_IOA_TN_FN.(flds0{idxBL})(1:len,1);
        len=min([length(pdata.stage2.avg_f1Traj_IOA_TN_FN.(fld)(:,1)),length(pdata.stage2.avg_f1Traj_IOA_TN_FN.(flds0{idxBL})(:,1))]);
        chgTrajF1_IOA_TN_FN.(flds0{k1}){end+1}=pdata.stage2.avg_f1Traj_IOA_TN_FN.(fld)(1:len,1)-...
            pdata.stage2.avg_f1Traj_IOA_TN_FN.(flds0{idxBL})(1:len,1);
        
        pertShiftF2_IOA.(fld){end+1}=pdata.stage2.avg_f2_pertShift.(fld)(:,1);
        pertShiftF1_IOA.(fld){end+1}=pdata.stage2.avg_f1_pertShift.(fld)(:,1);
        pertShiftF2_IOA_TN.(fld){end+1}=pdata.stage2.avg_f2_pertShift_IOA_TN(:,1).(fld)(:,1);
        pertShiftF1_IOA_TN.(fld){end+1}=pdata.stage2.avg_f1_pertShift_IOA_TN(:,1).(fld)(:,1);
        pertShiftF2_IOA_TN_FN.(fld){end+1}=pdata.stage2.avg_f2_pertShift_IOA_TN_FN(:,1).(fld)(:,1);
        pertShiftF1_IOA_TN_FN.(fld){end+1}=pdata.stage2.avg_f1_pertShift_IOA_TN_FN(:,1).(fld)(:,1);
               
    end
    
    % === ~Extract differences in avg. F1/F2 trajectories b/w the pert conditions and the noPert
    % condition ===
    
end

a_uyMidF2 = normT2_bw_uy_f2{find(config.normT2_bw_uy_ratio == 1/2)};
a_yu2MidF2 = normT3_bw_yu_f2{find(config.normT3_bw_yu_ratio == 1/2)};

%% Print the table of time intervals (noPert)
tab_intervals = {'a_IUInt', 'a_IYInt', 'a_IU2Int', ...
                 'a_IY2Int', 'a_IU3Int', 'a_IY3Int'};
fprintf('\n');
for i1 = 1 : numel(tab_intervals)
    eval(['t_ints = ', tab_intervals{i1}, ';']);
    
    t_int_mean = mean(t_ints(:, 1));
    t_int_std = std(t_ints(:, 1));
    fprintf('%s\t\t%.1f+/-%.1f ms\n', tab_intervals{i1}, ...
            t_int_mean * 1e3, t_int_std * 1e3);
end
fprintf('\n');

%% Print the table of F2 values (noPert)
tab_F2s = {'a_iF2', 'a_uF2', 'a_yF2', 'a_u2F2', ...
           'a_y2F2', 'a_u3F2', 'a_y3F2'};
fprintf('\n');
for i1 = 1 : numel(tab_F2s)
    eval(['t_F2s = ', tab_F2s{i1}, ';']);
    
    t_F2_mean = mean(t_F2s(:, 1));
    t_F2_std = std(t_F2s(:, 1));
    fprintf('%s\t\t%.0f+/-%.0f Hz\n', tab_F2s{i1}, t_F2_mean, t_F2_std);
end
fprintf('\n');

%%
flds=fields(avgTrajF2_pert_FTN_FN);
for i1=1:numel(flds)
    if ~isempty(avgTrajF2_pert_FTN_FN.(flds{i1}))
        t_avg=nanmean(avgTrajF2_pert_FTN_FN.(flds{i1}));
        t_std=nanstd(avgTrajF2_pert_FTN_FN.(flds{i1}));
        t_n=size(avgTrajF2_pert_FTN_FN.(flds{i1}),1)*ones(size(avgTrajF2_pert_FTN_FN.(flds{i1})));
    
        avgTrajF2_pert_FTN_FN.(flds{i1})=[t_avg',t_std',t_n'];
    end
end

if isequal(studyMode,'apstv2');
    nCols=5;
else
    nCols=3;
end

chg_uF2=a_uF2(:,2:nCols)-repmat(a_uF2(:,1),1,nCols-1);
chg_yF2=a_yF2(:,2:nCols)-repmat(a_yF2(:,1),1,nCols-1);
chg_u2F2=a_u2F2(:,2:nCols)-repmat(a_u2F2(:,1),1,nCols-1);
chg_u2y2MidF2=a_u2y2MidF2(:,2:nCols)-repmat(a_u2y2MidF2(:,1),1,nCols-1);
chg_y2F2=a_y2F2(:,2:nCols)-repmat(a_y2F2(:,1),1,nCols-1);

chg_uyMidF2 = a_uyMidF2(:,2:nCols)-repmat(a_uyMidF2(:,1),1,nCols-1);
chg_yu2MidF2 = a_yu2MidF2(:,2:nCols)-repmat(a_yu2MidF2(:,1),1,nCols-1);

chg_IUInt=a_IUInt(:,2:nCols)-repmat(a_IUInt(:,1),1,nCols-1);
chg_IYInt=a_IYInt(:,2:nCols)-repmat(a_IYInt(:,1),1,nCols-1);
chg_IU2Int=a_IU2Int(:,2:nCols)-repmat(a_IU2Int(:,1),1,nCols-1);
chg_IY2Int=a_IY2Int(:,2:nCols)-repmat(a_IY2Int(:,1),1,nCols-1);
chg_IU3Int=a_IU3Int(:,2:nCols)-repmat(a_IU3Int(:,1),1,nCols-1);
chg_IY3Int=a_IY3Int(:,2:nCols)-repmat(a_IY3Int(:,1),1,nCols-1);

chg_Y1Y2Int=a_Y1Y2Int(:,2:nCols)-repmat(a_Y1Y2Int(:,1),1,nCols-1);
chg_Y2Y3Int=a_Y2Y3Int(:,2:nCols)-repmat(a_Y2Y3Int(:,1),1,nCols-1);

chg_I_UYMid_Int=a_I_UYMid_Int(:,2:nCols)-repmat(a_I_UYMid_Int(:,1),1,nCols-1);
chg_I_YU2Mid_Int=a_I_YU2Mid_Int(:,2:nCols)-repmat(a_I_YU2Mid_Int(:,1),1,nCols-1);

chg_I_UYMidF_Int=a_I_UYMidF_Int(:,2:nCols)-repmat(a_I_UYMidF_Int(:,1),1,nCols-1);
chg_I_YU2MidF_Int=a_I_YU2MidF_Int(:,2:nCols)-repmat(a_I_YU2MidF_Int(:,1),1,nCols-1);

for i1=1:numel(normT1_bw_iu_f2)
    change_bw_iu_f2{i1}=normT1_bw_iu_f2{i1}(:,2:nCols)-repmat(normT1_bw_iu_f2{i1}(:,1),1,nCols-1);
end
for i1=1:numel(normT2_bw_uy_f2)
    change_bw_uy_f2{i1}=normT2_bw_uy_f2{i1}(:,2:nCols)-repmat(normT2_bw_uy_f2{i1}(:,1),1,nCols-1);
end
for i1=1:numel(normT3_bw_yu_f2)
    change_bw_yu_f2{i1}=normT3_bw_yu_f2{i1}(:,2:nCols)-repmat(normT3_bw_yu_f2{i1}(:,1),1,nCols-1);
end

%% Compute and plot the average chgTrajs
flds0=fields(chgTrajF1_POA);
for i1=1:length(flds0)
    fld=flds0{i1};
    if ~isempty(chgTrajF1_POA.(fld))
        avgChgTrajF1_POA.(fld)=avgTrace1(chgTrajF1_POA.(fld));
        avgChgTrajF2_POA.(fld)=avgTrace1(chgTrajF2_POA.(fld));
        avgChgTrajF1_POA_TN.(fld)=avgTrace1(chgTrajF1_POA_TN.(fld));
        avgChgTrajF2_POA_TN.(fld)=avgTrace1(chgTrajF2_POA_TN.(fld));
        avgChgTrajF1_POA_TN_FN.(fld)=avgTrace1(chgTrajF1_POA_TN_FN.(fld));
        avgChgTrajF2_POA_TN_FN.(fld)=avgTrace1(chgTrajF2_POA_TN_FN.(fld));
        avgChgTrajF1_IOA.(fld)=avgTrace1(chgTrajF1_IOA.(fld));
        avgChgTrajF2_IOA.(fld)=avgTrace1(chgTrajF2_IOA.(fld));
        avgChgTrajF1_IOA_TN.(fld)=avgTrace1(chgTrajF1_IOA_TN.(fld));
        avgChgTrajF2_IOA_TN.(fld)=avgTrace1(chgTrajF2_IOA_TN.(fld));
        
        avgChgTrajF2_IOA_FN.(fld)=avgTrace1(chgTrajF2_IOA_FN.(fld));
        avgChgTrajF1_IOA_FN.(fld)=avgTrace1(chgTrajF1_IOA_FN.(fld));        
        
        avgChgTrajF2_IOA_TN_FN.(fld)=avgTrace1(chgTrajF2_IOA_TN_FN.(fld));
        avgChgTrajF1_IOA_TN_FN.(fld)=avgTrace1(chgTrajF1_IOA_TN_FN.(fld));
        
        avgPertShiftF1_IOA.(fld)=avgTrace1(pertShiftF1_IOA.(fld));
        avgPertShiftF2_IOA.(fld)=avgTrace1(pertShiftF2_IOA.(fld));
        avgPertShiftF1_IOA_TN.(fld)=avgTrace1(pertShiftF1_IOA_TN.(fld));
        avgPertShiftF2_IOA_TN.(fld)=avgTrace1(pertShiftF2_IOA_TN.(fld));
        avgPertShiftF1_IOA_TN_FN.(fld)=avgTrace1(pertShiftF1_IOA_TN_FN.(fld));
        avgPertShiftF2_IOA_TN_FN.(fld)=avgTrace1(pertShiftF2_IOA_TN_FN.(fld));
        
        avg_trajF2_IOA.(fld)=avgTrace1(trajF2_IOA_NN.(fld));
        
        if ~(isequal(fld,'none') || isequal(fld,'noPert'))
            avgChgTrajF2_FTN.(fld)=avgTrace1(chgTrajF2_FTN.(fld));
            avgPertShiftF2_FTN.(fld)=avgTrace1(pertShiftF2_FTN.(fld));
            
            avgChgTrajF2_FTN_TA.(['aft_',fld])=avgTrace1(chgTrajF2_FTN_TA.(['aft_',fld]));            
        end
        
    end
end
avg_trajF2_IOA.none=avgTrace1(trajF2_IOA_NN.none);

avgChgTrajF2_POA_contrast = avgTrace1(chgTrajF2_POA.contrast);
avgChgTrajF2_IOA_contrast = avgTrace1(chgTrajF2_IOA.contrast);
avgChgTrajF2_FTN_contrast = avgTrace1(chgTrajF2_FTN.contrast);

frameDur=16/12e3*1e3;
normFrameDur=1/400;
SEM_ratio=1;

% plot_traj_changes(avgChgTrajF2_POA,frameDur,colors,'avgChg in trajF2_POA','SEM_ratio',SEM_ratio,...
%     'XLim',[0,500],'YLim',[-100,100]);
% plot_traj_changes(avgChgTrajF2_POA_TN_FN,normFrameDur,colors,'avgChg in trajF2_POA_TN_FN','SEM_ratio',SEM_ratio,...
%     'XLim',[0,1.5],'YLim',[-0.12,0.12]);
% plot_traj_changes(avgChgTrajF2_IOA,frameDur,colors,'avgChg in trajF2_IOA','SEM_ratio',SEM_ratio,...
%     'XLim',[0,500],'YLim',[-50,50],'windowHeight',250,'windowWidth',500,'boxOff','fontSize',fontSize,'timeUnit','ms');
% plot_traj_changes(avgChgTrajF2_IOA_TN,normFrameDur,colors,'avgChg in trajF2_IOA_TN','SEM_ratio',SEM_ratio,...
%     'XLim',[0,1.5],'YLim',[-100,100]);
% plot_traj_changes(avgChgTrajF2_IOA_TN_FN,normFrameDur,colors,'avgChg in trajF2_IOA_TN_FN','SEM_ratio',SEM_ratio,...
%     'XLim',[0,1.5],'YLim',[-0.12,0.12]);
% plot_traj_changes(avgPertShiftF2_IOA,frameDur,colors,'Pert of F2 (IOA)','SEM_ratio',SEM_ratio,...
%     'XLim',[0,500],'YLim',[-250,250]);
% plot_traj_changes(avgPertShiftF2_IOA_TN,normFrameDur,colors,'Pert of F2 (IOA_TN)','SEM_ratio',SEM_ratio,...
%     'XLim',[0,1.5],'YLim',[-250,250]);
% plot_traj_changes(avgPertShiftF2_IOA_TN_FN,normFrameDur,colors,'Pert of F2 (IOA_TN_FN)','SEM_ratio',SEM_ratio,...
%     'XLim',[0,1.5]);

% sigplot: plot the significance values (-log(p) values) explicitly
% This subroutine will compute the ratio of compensation in a hopefully more accurate way!
if isequal(studyMode, 'apstv') || isequal(studyMode, 'apstv_expanded') || isequal(studyMode, 'apstv2_stut_s.PFS')
    prefix = 'PFS_S';
elseif isequal(studyMode, 'apstv2t') || isequal(studyMode, 'apstv2t_expanded')
    prefix = 'PFS_T';
elseif ~isempty(fsic(varargin,'apstv2_stut_s.PWS'))
    prefix = 'PWS_S';
elseif ~isempty(fsic(varargin,'apstv2_stut_t.PWS'))
    prefix = 'PWS_T';
else
    prefix = 'prefix';
end

plot_traj_change_w_pert(avgPertShiftF2_IOA,avgChgTrajF2_IOA,pertShiftF2_IOA,chgTrajF2_IOA,semChgTrajF2_IOA,...
    frameDur,colors,'F2 (IOA)',prefix,'XLim',[0,500],'YLim',[-160,160],...
    'windowWidth',800,'windowHeight',800,'sigPlot',[-260,280; -280, 260]);
tLim = round(mean(a_IY2Int(:,1)) / (frameDur/1e3));
compute_correlation(avgChgTrajF2_IOA,chgTrajF2_IOA,tLim,frameDur,subjIDs);
find_repres_subj({a_uF2,a_uyMidF2,a_yF2,a_yu2MidF2},{a_IUInt,a_IYInt},subjIDs);

plot_traj_change_w_pert(avgPertShiftF2_IOA,avgChgTrajF2_IOA,pertShiftF2_IOA,chgTrajF2_IOA,semChgTrajF2_IOA,...
    frameDur,colors,'F2 (IOA)',prefix,'XLim',[0,500],'YLim',[-180,180],...
    'windowWidth',500,'windowHeight',260);
% plot_traj_change_w_pert_2sub(avg_trajF2_IOA,avgPertShiftF2_IOA,avgPertShiftF2_IOA,avgChgTrajF2_IOA,pertShiftF2_IOA,chgTrajF2_IOA,frameDur,colors,'F2 (IOA)','XLim',[0,500],'YLim',[-180,180],'windowWidth',800,'windowHeight',600);

% FTN individual-subject plots
frameDur_FTN_1 = 7 / 1495; 
plot_traj_change_w_pert(avgPertShiftF2_FTN,avgChgTrajF2_FTN,pertShiftF2_FTN,chgTrajF2_FTN,semChgTrajF2_FTN,...
    frameDur_FTN_1,colors,'F2 (FTN)',prefix,'XLim',[0, 7],'YLim',[-160,160],...
    'windowWidth',800,'windowHeight',800,'sigPlot',[-100,250; -250, 100]);

%% Fixed -effects analysis of the fullTNormF2 trajectories
if isfield(fe_chgTrajF2_FTN,'down');
    pertFld1='down';
    pertFld2='up';
else
    pertFld1='accel';
    pertFld2='decel';
end
fe_chgTrajF2_FTN_ste.(pertFld1)=std(fe_chgTrajF2_FTN.(pertFld1))/sqrt(size(fe_chgTrajF2_FTN.(pertFld1),1));
fe_chgTrajF2_FTN_ste.(pertFld2)=std(fe_chgTrajF2_FTN.(pertFld2))/sqrt(size(fe_chgTrajF2_FTN.(pertFld2),1));

tAxisN=pdata.stage2.fullTNormF2_tAxis;

figure;
subplot('Position',[0.1,0.4,0.85,0.55]);
plot(tAxisN,mean(fe_chgTrajF2_FTN.(pertFld1)),'-','color',colors.(pertFld1)); hold on;
plot(tAxisN,mean(fe_chgTrajF2_FTN.(pertFld1))-fe_chgTrajF2_FTN_ste.(pertFld1),'--','color',colors.(pertFld1));
plot(tAxisN,mean(fe_chgTrajF2_FTN.(pertFld1))+fe_chgTrajF2_FTN_ste.(pertFld1),'--','color',colors.(pertFld1));

plot(tAxisN,mean(fe_chgTrajF2_FTN.(pertFld2)),'color',colors.(pertFld2));
plot(tAxisN,mean(fe_chgTrajF2_FTN.(pertFld2))-fe_chgTrajF2_FTN_ste.(pertFld2),'--','color',colors.(pertFld2));
plot(tAxisN,mean(fe_chgTrajF2_FTN.(pertFld2))+fe_chgTrajF2_FTN_ste.(pertFld2),'--','color',colors.(pertFld2));
set(gca,'XLim',[tAxisN(1),tAxisN(end)]);
set(gca,'XTickLabel',{});
ylabel('Average F2 trajectory chg re noPert');

fe_p_t2s=nan(1,size(fe_chgTrajF2_FTN.(pertFld1),2));
for i1=1:numel(fe_p_t2s)
    [h,fe_p_t2s(i1)]=ttest2(fe_chgTrajF2_FTN.(pertFld1)(:,i1),fe_chgTrajF2_FTN.(pertFld2)(:,i1));
end

subplot('Position',[0.1,0.075,0.85,0.325]);
plot(tAxisN,-log10(fe_p_t2s)); hold on;
set(gca,'XLim',[tAxisN(1),tAxisN(end)]);

xs=get(gca,'XLim');
plot(xs,[2,2],'--','Color',[0.5,0.5,0.5]);
xlabel('Normalized time');
ylabel('-log10(p) (F.E)');

%% AT (adaptation test): F2 trajectories
figure('Name','Adaptation test: F2');
flds1=fields(avgChgTrajF2_FTN_TA);
legend_flds={};
legend_clrs={};
for i1=1:numel(flds1)
    fld=flds1{i1};
    if ~isempty(avgChgTrajF2_FTN_TA.(fld))
        
        plot(pdata.stage2.fullTNormF2_tAxis,avgChgTrajF2_FTN_TA.(fld)(:,1),'-','color',colors.(strrep(fld,'aft_','')));
        hold on;
        plot(pdata.stage2.fullTNormF2_tAxis,avgChgTrajF2_FTN_TA.(fld)(:,1)-avgChgTrajF2_FTN_TA.(fld)(:,2)./sqrt(avgChgTrajF2_FTN_TA.(fld)(:,3)),...
            '--','color',colors.(strrep(fld,'aft_','')));
        plot(pdata.stage2.fullTNormF2_tAxis,avgChgTrajF2_FTN_TA.(fld)(:,1)+avgChgTrajF2_FTN_TA.(fld)(:,2)./sqrt(avgChgTrajF2_FTN_TA.(fld)(:,3)),...
            '--','color',colors.(strrep(fld,'aft_','')));
        
        legend_flds{end+1}=strrep(fld,'_','\_');
        legend_clrs{end+1}=colors.(strrep(fld,'aft_',''));
    end
end
xs=[pdata.stage2.fullTNormF2_tAxis(1),pdata.stage2.fullTNormF2_tAxis(end)];
ys=get(gca,'YLim');
plot(xs,[0,0],'-','Color',[0.5,0.5,0.5]);
set(gca,'XLim',xs);
xlabel('Normalized time');
ylabel('Change from aft\_noPert');
ezlegend([xs(1)+0.1*range(xs),ys(1)+0.75*range(ys),0.3*range(xs),0.2*range(ys)],0.45,legend_flds,legend_clrs,...
    [12,12],{'-','-'},legend_clrs,[1,1],[0,0]);


%% Debug: APSTV: comparison between chg_yF2 and FTN changes
chg_yF2_FTN=struct;
flds1=fields(chgTrajF2_FTN);
for i1=1:numel(flds1)
    fld=flds1{i1};
    chg_yF2_FTN.(fld)=nan(numel(chgTrajF2_FTN.(fld)),1);
    for i2=1:numel(chgTrajF2_FTN.(fld))
        [foo,idx_2]=min(abs(pdata.stage2.fullTNormF2_tAxis-2));
        chg_yF2_FTN.(fld)(i2)=chgTrajF2_FTN.(fld){i2}(idx_2);
    end
end

if isfield(chg_yF2_FTN,'down')
    pertFld1='down';
    pertFld2='up';
else
    pertFld1='accel';
    pertFld2='decel';
end

figure;
plot(chg_yF2(:,1),chg_yF2_FTN.(pertFld1),'o','color',colors.(pertFld1));
axis equal
hold on;
for i1=1:numel(chg_yF2_FTN.(pertFld1))
    text(chg_yF2(i1,1),chg_yF2_FTN.(pertFld1)(i1),strrep(subjIDs{i1},'_','\_'),'FontSize',7,'Color',colors.(pertFld1));
end
plot(chg_yF2(:,2),chg_yF2_FTN.(pertFld2),'o','color',colors.(pertFld2));
for i1=1:numel(chg_yF2_FTN.(pertFld2))
    text(chg_yF2(i1,2),chg_yF2_FTN.(pertFld2)(i1),strrep(subjIDs{i1},'_','\_'),'FontSize',7,'Color',colors.up);
end
set(gca,'XLim',[-40,40],'YLim',[-40,40]);
plot([-40,40],[-40,40],'-','Color',[0.5,0.5,0.5]);
xlabel('chg\_yF2');
ylabel('chg\_yF2\_FTN');



%% Plot the full time-normalized F2 trajectories 
if isequal(studyMode,'apstv') || isequal(studyMode,'apstv2_stut_s')
    plot_traj_change_w_pert(avgPertShiftF2_FTN,avgChgTrajF2_FTN,pertShiftF2_FTN,chgTrajF2_FTN,[],frameDur_FTN,colors,'F2_FTN','','XLim',[0,6],'YLim',[-190,190],...
        'windowWidth',500,'windowHeight',500,'tPlot',[-260,280; -280, 260],'vertLines',[1:5],'arrows',{'D','E','F','G'},'labels');
elseif isequal(studyMode,'apstv2t') || isequal(studyMode,'apstv2_stut_t')
    plot_traj_change_w_pert(avgPertShiftF2_FTN,avgChgTrajF2_FTN,pertShiftF2_FTN,chgTrajF2_FTN,[],frameDur_FTN,colors,'F2_FTN','','XLim',[0,6],'YLim',[-180,260],...
        'windowWidth',500,'windowHeight',400,'tPlot',[-260,280; -280, 260],'vertLines',[1:5],'labels');
end

plot_traj_change_w_pert(avgPertShiftF2_FTN,avgChgTrajF2_FTN,pertShiftF2_FTN,chgTrajF2_FTN,[],frameDur_FTN,colors,'F2 (IOA)','','XLim',[0,6],'YLim',[-175,175],...
    'windowWidth',500,'windowHeight',260,'vertLines',[1:5]);
set(gca, 'YTick', [-150:50:150]);

% plot_traj_change_w_pert(avgPertShiftF2_FTN,avgChgTrajF2_FTN,pertShiftF2_FTN,chgTrajF2_FTN,frameDur_FTN,colors,'F2 (IOA)','XLim',[0,6],'YLim',[-200,200],...
%     'windowWidth',800,'windowHeight',800,'sigPlot',[-260,280; -280, 260],'vertLines',[1:5]);

%% Write to excel files (for SYSTAT use)
if (toWriteXls)
    t_xls_fn=fullfile(systatDir,'IUInt.xls');
    a_IUInt_cell=num2cell(a_IUInt);
    if isequal(studyMode,'apstv2')
        a_IUInt_cell=[{'noPert','Accel','Decel','Down','Up'};a_IUInt_cell];
    elseif isequal(studyMode,'apstv2t')  || isequal(studyMode,'apstv2_stut_t')
        a_IUInt_cell=[{'noPert','Accel','Decel'};a_IUInt_cell];
    end
    status=xlswrite(t_xls_fn,a_IUInt_cell);
    if status==1
        fprintf('%s successfully written.\n',t_xls_fn);    
    else
        fprintf('Writing to %s was unsuccessful.\n',t_xls_fn);   
    end
    
    t_xls_fn=fullfile(systatDir,'IYInt.xls');
    a_IYInt_cell=num2cell(a_IYInt);
    if isequal(studyMode,'apstv2')
        a_IYInt_cell=[{'noPert','Accel','Decel','Down','Up'};a_IYInt_cell];
    elseif isequal(studyMode,'apstv2t')  || isequal(studyMode,'apstv2_stut_t')
        a_IYInt_cell=[{'noPert','Accel','Decel'};a_IYInt_cell];
    end
        
    status=xlswrite(t_xls_fn,a_IYInt_cell);    
    if status==1
        fprintf('%s successfully written.\n',t_xls_fn);    
    else
        fprintf('Writing to %s was unsuccessful.\n',t_xls_fn);   
    end
    
    t_xls_fn=fullfile(systatDir,'I_UYMid_Int.xls');
    a_I_UYMid_Int_cell=num2cell(a_I_UYMid_Int);
    if isequal(studyMode,'apstv2')
        a_I_UYMid_Int_cell=[{'noPert','Accel','Decel','Down','Up'};a_I_UYMid_Int_cell];
    elseif isequal(studyMode,'apstv2t')  || isequal(studyMode,'apstv2_stut_t')
        a_I_UYMid_Int_cell=[{'noPert','Accel','Decel'};a_I_UYMid_Int_cell];
    end
    status=xlswrite(t_xls_fn,a_I_UYMid_Int_cell);    
    if status==1
        fprintf('%s successfully written.\n',t_xls_fn);    
    else
        fprintf('Writing to %s was unsuccessful.\n',t_xls_fn);   
    end
    
    t_xls_fn=fullfile(systatDir,'I_YU2Mid_Int.xls');
    a_I_YU2Mid_Int_cell=num2cell(a_I_YU2Mid_Int);
    if isequal(studyMode,'apstv2')
        a_I_YU2Mid_Int_cell=[{'noPert','Accel','Decel','Down','Up'};a_I_YU2Mid_Int_cell];
    elseif isequal(studyMode,'apstv2t')  || isequal(studyMode,'apstv2_stut_t')
        a_I_YU2Mid_Int_cell=[{'noPert','Accel','Decel'};a_I_YU2Mid_Int_cell];
    end
    status=xlswrite(t_xls_fn,a_I_YU2Mid_Int_cell);    
    if status==1
        fprintf('%s successfully written.\n',t_xls_fn);    
    else
        fprintf('Writing to %s was unsuccessful.\n',t_xls_fn);   
    end

    t_xls_fn=fullfile(systatDir,'I_UYMidF_Int.xls');
    a_I_UYMidF_Int_cell=num2cell(a_I_UYMidF_Int);
    if isequal(studyMode,'apstv2')
        a_I_UYMidF_Int_cell=[{'noPert','Accel','Decel','Down','Up'};a_I_UYMidF_Int_cell];
    elseif isequal(studyMode,'apstv2t')  || isequal(studyMode,'apstv2_stut_t')
        a_I_UYMidF_Int_cell=[{'noPert','Accel','Decel'};a_I_UYMidF_Int_cell];
    end
    status=xlswrite(t_xls_fn,a_I_UYMidF_Int_cell);
    if status==1
        fprintf('%s successfully written.\n',t_xls_fn);
    else
        fprintf('Writing to %s was unsuccessful.\n',t_xls_fn);
    end
    
    t_xls_fn=fullfile(systatDir,'I_YU2MidF_Int.xls');
    a_I_YU2MidF_Int_cell=num2cell(a_I_YU2MidF_Int);
    if isequal(studyMode,'apstv2')
        a_I_YU2MidF_Int_cell=[{'noPert','Accel','Decel','Down','Up'};a_I_YU2MidF_Int_cell];
    elseif isequal(studyMode,'apstv2t')  || isequal(studyMode,'apstv2_stut_t')
        a_I_YU2MidF_Int_cell=[{'noPert','Accel','Decel'};a_I_YU2MidF_Int_cell];
    end
        
    status=xlswrite(t_xls_fn,a_I_YU2MidF_Int_cell);
    if status==1
        fprintf('%s successfully written.\n',t_xls_fn);
    else
        fprintf('Writing to %s was unsuccessful.\n',t_xls_fn);
    end
    
%     t_xls_fn=fullfile(systatDir,'IUInt_Changes.xls');
%     status=xlswrite(t_xls_fn,[transpose(changeIUInt.accel),transpose(changeIUInt.decel)]);
%     if status==1
%         fprintf('%s successfully written.\n',t_xls_fn);    
%     else
%         fprintf('Writing to %s was unsuccessful.\n',t_xls_fn);
%     end
end

%% The POA trajectories
if doPOA
    metaTracePlot_(trajF2_POA,colors,'trajF2_POA','F2 (normalized)','Time (sec)',studyMode,...
        'XLim',[0,0.45],'YLim',[-0.1,1.1]);
%     metaTracePlot_(trajF1_POA,colors,'trajF1_POA','F1 (normalized)','Time (sec)','XLim',[0,0.45],'YLim',[-0.2,1.2]);
% meta%racePlot_(velF2_POA,colors,'velF2_POA','F2 (Hz)','Time (sec)','XLim',[0,0.35]);    
else
    fprintf('*** Skipped POA analysis. ***\n');
end

%% The IOA trajctories
flds=fields(trajF2_pert_IOA);
for i1=1:numel(flds)
    fld=flds{i1};
    if ~isempty(trajF2_pert_IOA.(fld))
        avg_trajF2_pert_IOA.(fld)=avgTrace1(trajF2_pert_IOA.(fld));
        avg_trajF2_pert_IOA_FN.(fld)=avgTrace1(trajF2_pert_IOA_FN.(fld));
    end
end

if doIOA
    metaTracePlot_(trajF2_IOA,colors,'trajF2_IOA','F2 (normalized)','Time (sec)',studyMode,...
        'XLim',[0,0.45],'YLim',[-0.1,1.1],'xlabel','Time from [i] (ms)',opt_pertBar,[10,280],opt_OA,avg_trajF2_pert_IOA_FN,...
        'showArrows',{'(C)','(D)','(E)','(F)'}, 'showPert', avg_trajF2_pert_IOA_FN);
    metaTracePlot_(trajF2_IOA,colors,'trajF2_IOA','F2 (normalized)','Time (sec)',studyMode,...
        'XLim',[0,0.45],'YLim',[-0.1,1.1],'xlabel','Time from [i] (ms)',opt_pertBar,[10,280],opt_OA,avg_trajF2_pert_IOA_FN,...
        'showArrows',{'(C)','(D)','(E)','(F)'});
%     metaTracePlot_(trajF1_IOA,colors,'trajF1_IOA','F1 (normalized)','Time (sec)',studyMode,...
%         'XLim',[0,0.45],'YLim',[-0.1,1.1],'xlabel','Time from [i] (ms)',opt_pertBar,[10,280],opt_OA);
    metaTracePlot_(avgTrajF2_FTN_FN,colors,'trajF2_FTN','F2 (normalized)','Time (sec)',studyMode,...
        'lw',0.5,'XLim',[0,3],'YLim',[-0.18,1.04],'xlabel','Normalized time',opt_pertBar,[10,280],opt_OA,avgTrajF2_pert_FTN_FN);
metaTracePlot_(trajF2_IOA,colors,'trajF2_IOA','F2 (normalized)','Time (sec)',studyMode,...
        'XLim',[0,0.45],'YLim',[-0.1,1.1],'xlabel','Time from [i] (ms)',opt_pertBar,[10,280],opt_OA,avg_trajF2_pert_IOA_FN);
else
    fprintf('*** Skipped IOA analysis. ***\n');
end

%% POA_TN trajectories
if doPOA_TN
    metaTracePlot_(trajF2_POA_TN,colors,'trajF2_POA_TN','F2 (normalized)','Time (normalized)',studyMode,...
        'XLim',[0,1.75],'YLim',[-0.1,1.1]);
else
    fprintf('*** Skipped POA_TN analysis. ***\n');
end

%% FTN_FN trajectories;
metaTracePlot_(avgTrajF2_FTN_FN,colors,'trajF2_IOA_FTN (RE)','F2','Time (normalized)',studyMode,...
    'XLim',[0,3],'YLim',[-0.05,1.1],'lw',0.5); % FTN stands for full time normalizatoin

%% IOA_TN trajectories
if doIOA_TN
    metaTracePlot_(trajF2_IOA_TN,colors,'trajF2_IOA_TN (FE)','F2 (normalized)','Time (normalized)',studyMode,...
        'XLim',[0,1.75],'YLim',[-0.1,1.1]);
%     metaTracePlot_(avgTrajF2_FTN,colors,'trajF2_IOA_TN (RE)','F2','Time (normalized)',studyMode,...
%         'XLim',[0,1.75]); % FTN stands for full time normalizatoin

else
    fprintf('*** Skipped POA_TN analysis. ***\n');
end

%% Show the RMS trajectories
% metaTracePlot_(trajRMS_IOA,colors,'trajRMS_IOA','RMS (demeaned, normed)','Time (sec)',studyMode);
% set(gca, 'YLim', [-1, 1], 'XLim', [0, 500]);
% xlabel('Time re [i] (ms)');

%% Output and exit (if applicable)
if ~isempty(fsic(varargin,'output'))
    outputItems=varargin{fsic(varargin,'output')+1};
    
    if ~iscell(outputItems)
        outputItems={outputItems};
    end
    for i1=1:numel(outputItems)
        eval(sprintf('varargout{%d}=%s;',i1,outputItems{i1}));
    end
    close all;
    return;
end


%% Percentage of compensation: 
if ~(isequal(studyMode,'apstv') || isequal(studyMode,'apstv2_stut_s'))
    figure('Position',[100,100,720,360]);
    subplot(1,2,1);
    plot(a_tShift_iuMidF.accel*1e3,chg_I_UYMidF_Int(:,1)*1e3,'o'); hold on;
    ylabel('Compensation in [i]-[u_1,j_1]_m_i_d_F (ms)','FontSize',fontSize);
    xlabel('Time shift of [i,u_1]_m_i_d_F (ms)','FontSize',fontSize);
    [k,r2,p]=lincorr(a_tShift_iuMidF.accel*1e3,chg_I_UYMidF_Int(:,1)*1e3);
    xs=get(gca,'XLim'); ys=get(gca,'YLim');
    plot(xs,[0,0],'color',[0.5,0.5,0.5]);
    plot(xs,k(1)+k(2)*xs,'k--');
    set(gca,'XLim',xs,'YLim',ys,'FontSize',fontSize);
    text(xs(1)+0.05*range(xs),ys(2)-0.06*range(ys),sprintf('Linear corr.: r^2=%.3f; p=%.3f',r2,p));
    ratioTComp.accel=mean(chg_I_UYMidF_Int(:,1))/mean(a_tShift_iuMidF.accel);
    [r,t,p]=spear(transpose(a_tShift_iuMidF.accel*1e3),chg_I_UYMidF_Int(:,1)*1e3);
    text(xs(1)+0.05*range(xs),ys(2)-0.12*range(ys),sprintf('Spearman corr.: p=%.3f',p));
    text(xs(1)+0.05*range(xs),ys(2)-0.18*range(ys),sprintf('Pct. of compensation: %.1f%%',ratioTComp.accel*1e2));

    subplot(1,2,2);
    plot(a_tShift_iuMidF.decel*1e3,chg_I_UYMidF_Int(:,2)*1e3,'o'); hold on;
    ylabel('Compensation in [i]-[u_1,j_1]_m_i_d_F (ms)','FontSize',fontSize);
    xlabel('Time shift of [i,u_1]_m_i_d_F (ms)','FontSize',fontSize);
    [k,r2,p]=lincorr(a_tShift_iuMidF.decel*1e3,chg_I_UYMidF_Int(:,2)*1e3);
    xs=get(gca,'XLim'); ys=get(gca,'YLim');
    plot(xs,[0,0],'color',[0.5,0.5,0.5]);
    plot(xs,k(1)+k(2)*xs,'k--');
    set(gca,'XLim',xs,'YLim',ys,'FontSize',fontSize);
    text(xs(1)+0.05*range(xs),ys(2)-0.06*range(ys),sprintf('Linear corr.: r^2=%.3f; p=%.3f',r2,p));
    ratioTComp.decel=mean(chg_I_UYMidF_Int(:,2))/mean(a_tShift_iuMidF.decel);
    [r,t,p]=spear(transpose(a_tShift_iuMidF.decel*1e3),chg_I_UYMidF_Int(:,2)*1e3);
    text(xs(1)+0.05*range(xs),ys(2)-0.12*range(ys),sprintf('Spearman corr.: p=%.3f',p));
    text(xs(1)+0.05*range(xs),ys(2)-0.18*range(ys),sprintf('Pct. of compensation: %.1f%%',ratioTComp.decel*1e2));
end

%% Check the consistency between the time intervals extracted from the individual and average
%% trajectories
if isequal(studyMode, 'apstv')
    chg_IUInt_aTraj(:, 1) = ioaTraj_tInts.down(:, 1) - ioaTraj_tInts.none(:, 1);
    chg_IUInt_aTraj(:, 2) = ioaTraj_tInts.up(:, 1) - ioaTraj_tInts.none(:, 1);
    chg_IYInt_aTraj(:, 1) = ioaTraj_tInts.down(:, 2) - ioaTraj_tInts.none(:, 2);
    chg_IYInt_aTraj(:, 2) = ioaTraj_tInts.up(:, 2) - ioaTraj_tInts.none(:, 2);
    chg_IU2Int_aTraj(:, 1) = ioaTraj_tInts.down(:, 3) - ioaTraj_tInts.none(:, 3);
    chg_IU2Int_aTraj(:, 2) = ioaTraj_tInts.up(:, 3) - ioaTraj_tInts.none(:, 3);

    figure('Position', [100, 50, 1000, 700]);
    subplot(2, 3, 1);
    plot_eq(chg_IUInt_aTraj(:, 1), chg_IUInt(:, 1), 'o', 'chg_iuInt aTraj', 'chg_iuInt', 'down');
    subplot(2, 3, 2);
    plot_eq(chg_IYInt_aTraj(:, 1), chg_IYInt(:, 1), 'o', 'chg_iyInt aTraj', 'chg_iyInt', 'down');
    subplot(2, 3, 3);
    plot_eq(chg_IU2Int_aTraj(:, 1), chg_IU2Int(:, 1), 'o', 'chg_iu2Int aTraj', 'chg_iu2Int', 'down');
    subplot(2, 3, 4);
    plot_eq(chg_IUInt_aTraj(:, 2), chg_IUInt(:, 2), 'o', 'chg_iuInt aTraj', 'chg_iuInt', 'up');
    subplot(2, 3, 5);
    plot_eq(chg_IYInt_aTraj(:, 2), chg_IYInt(:, 2), 'o', 'chg_iyInt aTraj', 'chg_iyInt', 'up');
    subplot(2, 3, 6);
    plot_eq(chg_IU2Int_aTraj(:, 2), chg_IU2Int(:, 2), 'o', 'chg_iu2Int aTraj', 'chg_iu2Int', 'up');
end

%% Visualization 1.1: three-level plots
% if ~(isequal(studyMode,'apstv')  || isequal(studyMode,'apstv2_stut_s'))
    t_int_sum(chg_IUInt,chg_IYInt,chg_IU2Int,chg_IY2Int,chg_IU3Int,chg_IY3Int,colors,...
        'windowWidth',500,'windowHeight',260,'fontSize',12);
% end

if isequal(studyMode,'apstv2')
    metaPlot_5(a_IUInt*1e3,chg_IUInt*1e3,'[i]-[u]_1 interval','ms',fontSize,colors);
    metaPlot_5(a_IYInt*1e3,chg_IYInt*1e3,'[i]-[j]_1 interval','ms',fontSize,colors);
    metaPlot_5(a_IU2Int*1e3,chg_IU2Int*1e3,'[i]-[u]_2 interval','ms',fontSize,colors);
    metaPlot_5(a_IY2Int*1e3,chg_IY2Int*1e3,'[i]-[j]_2 interval','ms',fontSize,colors);
    metaPlot_5(a_IU3Int*1e3,chg_IU3Int*1e3,'[i]-[u]_3 interval','ms',fontSize,colors);
    metaPlot_5(a_IY3Int*1e3,chg_IY3Int*1e3,'[i]-[j]_2 interval','ms',fontSize,colors);
    
    metaPlot_5(a_I_UYMid_Int*1e3,chg_I_UYMid_Int*1e3,'[i]-[u,j] interval','ms',fontSize,colors,'YLim',[205,240]);
    metaPlot_5(a_I_YU2Mid_Int*1e3,chg_I_YU2Mid_Int*1e3,'[i]-[j,u_2] interval','ms',fontSize,colors,'YLim',[350,410]);

    metaPlot_5(a_I_UYMidF_Int*1e3,chg_I_UYMidF_Int*1e3,'[i]-[u_1,j_1]_m_i_d_F interval','ms',fontSize,colors,'YLim',[200,250]);
    metaPlot_5(a_I_YU2MidF_Int*1e3,chg_I_YU2MidF_Int*1e3,'[i]-[j_1,u_2]_m_i_d_F interval','ms',fontSize,colors,'YLim',[360,420]);

    metaPlot_5(a_uF2,chg_uF2,'F2 at [u]_1','Hz',fontSize,colors,'YLim',[1075,1225]);
    metaPlot_5(a_yF2,chg_yF2,'F2 at [j]_1','Hz',fontSize,colors,'YLim',[2000,2125]);

    idx_uy_mid=find(config.normT2_bw_uy_ratio==0.5);
    metaPlot_5(normT2_bw_uy_f2{idx_uy_mid},change_bw_uy_f2{idx_uy_mid},'F2 at the [u]_1-[j]_1 midpoint','Hz',fontSize,colors,'YLim',[1625,1700]);

    idx_yu_1q=find(config.normT3_bw_yu_ratio==0.25);
    metaPlot_5(normT3_bw_yu_f2{idx_yu_1q},change_bw_yu_f2{idx_yu_1q},'F2 at the [j]_1-[u]_2 1st quater point','Hz',fontSize,colors,'YLim',[1900,2100]);

    idx_yu_mid=find(config.normT3_bw_yu_ratio==0.5);
    metaPlot_5(normT3_bw_yu_f2{idx_yu_mid},change_bw_yu_f2{idx_yu_mid},'F2 at the [j]_1-[u]_2 midpoint','Hz',fontSize,colors,'YLim',[1675,1775]);
    
else
    metaPlot_3(studyMode,a_IUInt*1e3,chg_IUInt*1e3,'[i]-[u]_1 interval','ms',fontSize,colors,posthocOpt,...
        multiSubOpt,'windowWidth',ana_winWidth,'windowHeight',ana_winHeight*1.35,'horizontal');
    metaPlot_3(studyMode,a_IYInt*1e3,chg_IYInt*1e3,'[i]-[j]_1 interval','ms',fontSize,colors,posthocOpt,...
        multiSubOpt,'windowWidth',ana_winWidth,'windowHeight',ana_winHeight*1.35,'horizontal');
    metaPlot_3(studyMode,a_IU2Int*1e3,chg_IU2Int*1e3,'[i]-[u]_2 interval','ms',fontSize,colors,posthocOpt,...
        multiSubOpt,'windowWidth',ana_winWidth,'windowHeight',ana_winHeight*1.35,'horizontal');
    metaPlot_3(studyMode,a_IY2Int*1e3,chg_IY2Int*1e3,'[i]-[j]_2 interval','ms',fontSize,colors,posthocOpt,...
        multiSubOpt,'windowWidth',ana_winWidth,'windowHeight',ana_winHeight*1.35,'horizontal');
    metaPlot_3(studyMode,a_IU3Int*1e3,chg_IU3Int*1e3,'[i]-[u]_3 interval','ms',fontSize,colors,posthocOpt,...
        multiSubOpt,'windowWidth',ana_winWidth,'windowHeight',ana_winHeight*1.35,'horizontal');
    metaPlot_3(studyMode,a_IY3Int*1e3,chg_IY3Int*1e3,'[i]-[j]_2 interval','ms',fontSize,colors,posthocOpt,...
        multiSubOpt,'windowWidth',ana_winWidth,'windowHeight',ana_winHeight*1.35,'horizontal');

    metaPlot_3(studyMode,a_Y1Y2Int*1e3,chg_Y1Y2Int*1e3,'[j]_1-[j]_2 interval','ms',fontSize,colors,posthocOpt,...
               'YLim',[270,320],multiSubOpt,'windowWidth',ana_winWidth,'windowHeight',ana_winHeight);
    metaPlot_3(studyMode,a_Y2Y3Int*1e3,chg_Y2Y3Int*1e3,'[j]_2-[j]_3 interval','ms',fontSize,colors,posthocOpt,...
               'YLim',[270,320],multiSubOpt,'windowWidth',ana_winWidth,'windowHeight',ana_winHeight);
    
    plot_interval_changes(a_IYInt*1e3,a_Y1Y2Int*1e3,a_Y2Y3Int*1e3,chg_IYInt*1e3,chg_Y1Y2Int*1e3,chg_Y2Y3Int*1e3,studyMode,fontSize,colors,posthocOpt);
    
    
    metaPlot_3(studyMode,a_I_UYMid_Int*1e3,chg_I_UYMid_Int*1e3,'[i]-[u,j] interval','ms',fontSize,colors,posthocOpt,'YLim',[205,240],multiSubOpt,'windowWidth',ana_winWidth,'windowHeight',ana_winHeight);
    metaPlot_3(studyMode,a_I_YU2Mid_Int*1e3,chg_I_YU2Mid_Int*1e3,'[i]-[j,u_2] interval','ms',fontSize,colors,posthocOpt,'YLim',[350,410],multiSubOpt,'windowWidth',ana_winWidth,'windowHeight',ana_winHeight);
    
    metaPlot_3(studyMode,a_I_UYMidF_Int*1e3,chg_I_UYMidF_Int*1e3,'[i]-[u_1,j_1]_m_i_d_F interval','ms',fontSize,colors,posthocOpt,'YLim',[200,250],multiSubOpt,'windowWidth',ana_winWidth,'windowHeight',ana_winHeight);
    metaPlot_3(studyMode,a_I_YU2MidF_Int*1e3,chg_I_YU2MidF_Int*1e3,'[i]-[j_1,u_2]_m_i_d_F interval','ms',fontSize,colors,posthocOpt,'YLim',[360,420],multiSubOpt,'windowWidth',ana_winWidth,'windowHeight',ana_winHeight);

    idx_uy_mid=find(config.normT2_bw_uy_ratio==0.5);
    idx_yu_1q=find(config.normT3_bw_yu_ratio==0.25);
    idx_yu_mid=find(config.normT3_bw_yu_ratio==0.5);
    idx_u2=find(config.normT3_bw_yu_ratio==1.0);
    
    commonYLim = [-27, 48];
    metaPlot_3(studyMode,a_uF2,chg_uF2,'F2 at [u]_1','Hz',fontSize,colors,posthocOpt,...
               'YLim',commonYLim,multiSubOpt,'windowWidth',ana_winWidth,'windowHeight',ana_winHeight);
    metaPlot_3(studyMode,normT2_bw_uy_f2{idx_uy_mid},change_bw_uy_f2{idx_uy_mid},'F2 at [u]_1-[j]_1 midpoint','Hz',...
               fontSize,colors,posthocOpt,'YLim',commonYLim,multiSubOpt,'windowWidth',ana_winWidth,'windowHeight',ana_winHeight);
    metaPlot_3(studyMode,a_yF2,chg_yF2,'F2 at [j]_1','Hz',fontSize,colors,posthocOpt,...
               'YLim',commonYLim,multiSubOpt,'windowWidth',ana_winWidth,'windowHeight',ana_winHeight);
    metaPlot_3(studyMode,normT3_bw_yu_f2{idx_yu_mid},change_bw_yu_f2{idx_yu_mid},'F2 at  [j]_1-[u]_2 midpoint','Hz',...
               fontSize,colors,posthocOpt,'YLim',commonYLim,multiSubOpt,'windowWidth',ana_winWidth,'windowHeight',ana_winHeight);

%     metaPlot_3(studyMode,a_u2F2,chg_u2F2,'F2 at [u]_2','Hz',fontSize,colors,posthocOpt,'YLim',[1000,1150],multiSubOpt,'windowWidth',ana_winWidth,'windowHeigth',ana_winHeight);
    metaPlot_3(studyMode,normT3_bw_yu_f2{idx_yu_1q},change_bw_yu_f2{idx_yu_1q},'F2 at [j]-[u]_2 1st quater point','Hz',...
        fontSize,colors,posthocOpt,'YLim',[1900,2100],multiSubOpt,'windowWidth',ana_winWidth,'windowHeight',ana_winHeight);
    metaPlot_3(studyMode,normT3_bw_yu_f2{idx_u2},change_bw_yu_f2{idx_u2},'F2 at [u]_2','Hz',...
        fontSize,colors,posthocOpt,'YLim',[900,1200],multiSubOpt,'windowWidth',ana_winWidth,'windowHeight',ana_winHeight);

%     metaPlot_3_4items(studyMode, a_uF2, chg_uF2, normT2_bw_uy_f2{idx_uy_mid}, change_bw_uy_f2{idx_uy_mid} , ...
%                       a_yF2, chg_yF2, normT3_bw_yu_f2{idx_yu_mid}, change_bw_yu_f2{idx_yu_mid});
                      
end


%% Visualization 1.2: two-level plots
% metaPlot_2(changeIF2,'Change in F2 of [i] in "I"','Hz',fontSize,colors);
% metaPlot_2(changeYF2,'Change in F2 of [j] in "you"','Hz',fontSize,colors);
% tScore_int(1)=metaPlot_2(changeIUInt,'Change in [i]-[u]_1 interval','ms',fontSize,colors,'YLim',[-3,6]);
% tScore_int_UY=metaPlot_2(changeUYInt,'Change in [u]_1-[j]_1 interval','ms',fontSize,colors,'YLim',[-6,6]);
% tScore_int(2)=metaPlot_2(changeIYInt,'Change in [i]-[j]_1 interval','ms',fontSize,colors,'YLim',[-6,8]);
% tScore_int(3)=metaPlot_2(changeIU2Int,'Change in [i]-[u]_2 interval','ms',fontSize,colors,'YLim',[-6,12]);
% tScore_int(4)=metaPlot_2(changeIY2Int,'Change in [i]-[j]_2 interval','ms',fontSize,colors,'YLim',[-6,12]);

%% Visualization 1.3: the decay of local time perturbation 
if exist('tScore')
    figure('NumberTitle','off','Name','Decay of local time perturbation','Position',[100,100,500,300]);
    set(gca,'FontSize',fontSize,'FontWeight','Bold','LineWidth',2);
    plot([1,2,3,4],-tScore_int,'bo-','LineWidth',2);
    hold on;
    plot([0.5,4.5],repmat(-tinv(0.05/2,2*(length(changeIUInt.accel)-1)),1,2),'k--','LineWidth',1)
    text(4,-tinv(0.05/2,2*(length(changeIUInt.accel)-1)),'p=0.05','FontSize',fontSize,'FontWeight','Bold');
    plot([0.5,4.5],repmat(-tinv(0.01/2,2*(length(changeIUInt.accel)-1)),1,2),'k--','LineWidth',1)
    text(4,-tinv(0.01/2,2*(length(changeIUInt.accel)-1)),'p=0.01','FontSize',fontSize,'FontWeight','Bold');
    box off;
    set(gca,'XLim',[0.5,4.5]);
    set(gca,'YLim',[1.25,3.25]);
    set(gca,'XTick',[1:4],'XTickLabel',[]);
    ylabel('Decel vs. Accel t-score','FontSize',fontSize,'FontWeight','Bold');
    set(gca,'FontSize',fontSize,'FontWeight','Bold','LineWidth',2);
end

%%
t_norm_tab.mean=nan(0,3);
t_norm_tab.sem=nan(0,3);
% 
% for i1=1:length(config.normT1_bw_iu_ratio)
%     [t_means,t_sems]=metaPlot_errorBar(normT1_bw_iu_f2{i1},sprintf('normT1: F2 @ %.2f point between i and u',config.normT1_bw_iu_ratio(i1)),'Hz',fontSize);
%     t_norm_tab.mean=[t_norm_tab.mean;t_means];
%     t_norm_tab.sem=[t_norm_tab.sem;t_sems];
% end
% % for i1=1:length(config.normT1_bw_uy_ratio)
% %     [t_means,t_sems]=metaPlot_errorBar(normT1_bw_uy_f2{i1},sprintf('normT1: F2 @ %.2f point between u and y',config.normT1_bw_uy_ratio(i1)),'Hz',fontSize);
% %     t_norm_tab.mean=[t_norm_tab.mean;t_means];
% %     t_norm_tab.sem=[t_norm_tab.sem;t_sems];    
% % end
% for i1=1:length(config.normT2_bw_uy_ratio)
%     [t_means,t_sems]=metaPlot_errorBar(normT2_bw_uy_f2{i1},sprintf('normT2: 2 @ %.2f point between u and y',config.normT2_bw_uy_ratio(i1)),'Hz',fontSize);
%     t_norm_tab.mean=[t_norm_tab.mean;t_means];
%     t_norm_tab.sem=[t_norm_tab.sem;t_sems];  
% end
% for i1=1:length(config.normT3_bw_yu_ratio)
%     [t_means,t_sems]=metaPlot_errorBar(normT3_bw_yu_f2{i1},sprintf('normT3: F2 @ %.2f point between y and u in "you"',config.normT3_bw_yu_ratio(i1)),'Hz',fontSize);
%     t_norm_tab.mean=[t_norm_tab.mean;t_means];
%     t_norm_tab.sem=[t_norm_tab.sem;t_sems];
% end
% 
% for i1=1:length(config.normT1_bw_iu_ratio)
% %     metaPlot_2(change_bw_iu_f2{i1},sprintf('Change in F2 @ %.2f normT1 between i and u',config.normT1_bw_iu_ratio(i1)),'Hz',fontSize,'YLim',[-10,15]);
%     if i1<length(config.normT1_bw_iu_ratio)
%         metaPlot_2(change_bw_iu_f2{i1},sprintf('Change in F2'),'Hz',fontSize,colors,'YLim',[-10,15]);
%     else
%         metaPlot_2(change_bw_iu_f2{i1},sprintf('Change in F2'),'Hz',fontSize,colors,'YLim',[-50,40]);
%     end
% end
% for i1=1:length(config.normT1_bw_uy_ratio)
% %     metaPlot_2(change_bw_uy_f2{i1},sprintf('Change in F2 @ %.2f normT2 between u and y',config.normT1_bw_uy_ratio(i1)),'Hz',fontSize);
%     metaPlot_2(change_bw_uy_f2{i1},sprintf('Change in F2'),'Hz',fontSize,colors,'YLim',[-50,40]);
% end
% for i1=1:length(config.normT3_bw_yu_ratio)
% %     metaPlot_2(change_bw_yu_f2{i1},sprintf('Change in F2 @ %.2f normT3 between y and u in "you"',config.normT3_bw_yu_ratio(i1)),'Hz',fontSize);
%     metaPlot_2(change_bw_yu_f2{i1},sprintf('Change in F2'),'Hz',fontSize,colors,'YLim',[-30,60]);
% end

%% Visualization 2
% figure('Position',[100,100,900,360])
% set(gca,'FontSize',fontSize,'FontWeight','Bold','LineWidth',2);
% flds=fields(colors);
% for i1=1:size(t_norm_tab.mean,2)
%     plot(1:size(t_norm_tab.mean,1),t_norm_tab.mean(:,i1),'o-','Color',colors.(flds{i1}),'LineWidth',1);
%     hold on;
% end
% % set(gca,'XLim',[1,size(t_norm_tab.mean,1)],'YLim',[1050,2250],'XTick',[1,5,9,13],'XTickLabel',{'[i]','[u]_1','[j]_1','[u]_2'},'FontWeight','Bold');
% set(gca,'XLim',[1,size(t_norm_tab.mean,1)],'YLim',[1050,2270],'XTick',[]);
% ys=get(gca,'YLim');
% text(1,ys(1)-0.05*range(ys),'[i]','FontSize',fontSize,'FontWeight','Bold');
% text(5,ys(1)-0.05*range(ys),'[u]_1','FontSize',fontSize,'FontWeight','Bold');
% text(9,ys(1)-0.05*range(ys),'[j]_1','FontSize',fontSize,'FontWeight','Bold');
% text(13,ys(1)-0.05*range(ys),'[u]_2','FontSize',fontSize,'FontWeight','Bold');
% xlabel('Normalized time','FontSize',fontSize,'FontWeight','Bold');
% ylabel('F2 (Hz)','FontSize',fontSize,'FontWeight','Bold');
% legend({'noPert','accel','decel'},'Location','Northeast');
% 
% plot([5,5],ys,'-','Color',[0.5,0.5,0.5]);
% plot([9,9],ys,'-','Color',[0.5,0.5,0.5]);
% box on;

return


%% Sub-functions
function varargout=metaPlot_errorBar(meas,measName,measUnit,fontSize,varargin)
flds=fields(meas);
x_tick_label=cell(1,0);
meas_mean=[];
meas_sem=[];
rmanova_tab=[];
for i1=1:length(flds)
	t_meas=meas.(flds{i1});
	x_tick_label{length(x_tick_label)+1}=flds{i1};
	meas_mean=[meas_mean,mean(t_meas(~isnan(t_meas)))];
	meas_sem=[meas_sem,ste(t_meas(~isnan(t_meas)))];
    rmanova_tab=[rmanova_tab;...
        [transpose(t_meas),i1*ones(length(t_meas),1),transpose(1:length(t_meas))]];
end
figure('NumberTitle','off','name',measName);
set(gca,'FontSize',fontSize,'LineWidth',2);
errorbar([1:length(flds)],meas_mean,meas_sem,'LineWidth',2);
set(gca,'XLim',[0.5,length(flds)+0.5],'XTick',[1:length(flds)]);
set(gca,'XTickLabel',x_tick_label);
ylabel([measName,' (',measUnit,')'],'FontSize',fontSize);
title(measName,'FontSize',12);

% Perform one-way RM-ANOVA
rmanova_res=rmaov1(rmanova_tab,0.05);
xs=get(gca,'XLim'); ys=get(gca,'YLim');
text(xs(1)+0.05*range(xs),ys(2)-0.075*range(ys),...
    sprintf('RM-ANOVA: F(%d) = %.3f, p = %.4f',rmanova_res.v1,rmanova_res.F1,rmanova_res.P1),...
    'FontSize',fontSize);

% Output
if nargout==2
    varargout{1}=meas_mean;
    varargout{2}=meas_sem;
end
return

%%
function varargout=metaPlot_2(meas,measName,measUnit,fontSize,colors,varargin)
figure('NumberTitle','off','Name',measName,'Position',[100,100,300,200]);
set(gca,'FontSize',fontSize,'LineWidth',1);
% bar([1,2],[mean(meas.accel),mean(meas.decel)]); hold on;
bar(1,mean(meas.accel),'EdgeColor','k','FaceColor',colors.accel,'LineWidth',1); hold on;
bar(2,mean(meas.decel),'EdgeColor','k','FaceColor',colors.decel,'LineWidth',1);
errorbar([1,2],[mean(meas.accel),mean(meas.decel)],[ste(meas.accel),ste(meas.decel)],'ko','LineWidth',1);
set(gca,'XTick',[1,2],'XTickLabel',{'AP','DP'});
% set(gca,'YLim',[-20,20]);
ylabel(sprintf('%s (%s)',measName,measUnit),'FontSize',fontSize-2);

if ~isempty(fsic(varargin,'YLim'))
    idx=fsic(varargin,'YLim');
    set(gca,'YLim',varargin{idx+1});
end
box off;

% Perform RM-ANOVA
rmanova_tab=[];
flds=fields(meas);
for i1=1:length(flds)
	t_meas=meas.(flds{i1});
    rmanova_tab=[rmanova_tab;...
        [transpose(t_meas),i1*ones(length(t_meas),1),transpose(1:length(t_meas))]];
end
rmanova_res=rmaov1(rmanova_tab,0.05);
xs=get(gca,'XLim'); ys=get(gca,'YLim');
% text(xs(1)+0.05*range(xs),ys(2)-0.075*range(ys),...
%     sprintf('RM-ANOVA: F(%d,%d) = %.3f, p = %.4f',rmanova_res.v1,rmanova_res.v3,rmanova_res.F1,rmanova_res.P1),...
%     'FontSize',fontSize);

% Perform paired t-test
[h,p_t,ci_t,stats_t]=ttest2(meas.accel,meas.decel);
fprintf('2-sample t-test: t(%d) = %.3f, p = %.3f\n',stats_t.df,stats_t.tstat,p_t);
% text(xs(1)+0.05*range(xs),ys(2)-0.150*range(ys),...
%     sprintf('Paired t-test: t(%d) = %.3f, p = %.4f',stats_t.df,stats_t.tstat,p_t),...
%     'FontSize',fontSize);

y_bar=ys(2)-0.05*range(ys);
y_bar_low=ys(2)-0.065*range(ys);

plot([1,2],[y_bar,y_bar],'k-','LineWidth',1);
plot([1,1],[y_bar_low,y_bar],'k-','LineWidth',1);
plot([2,2],[y_bar_low,y_bar],'k-','LineWidth',1);
if (rmanova_res.P1<0.05);
    plot(1.5,y_bar+0.025*range(ys),'k*','MarkerSize',9,'LineWidth',2);
end
set(gca,'FontSize',fontSize,'LineWidth',1);

if nargout==1
    varargout{1}=stats_t.tstat;
end
return

%% metaPlot_3_line
function metaPlot_3_line(meas,measName,measUnit,fontSize,varargin)
figure('NumberTitle','off','Name',measName,'Position',[100,100,400,300]);
set(gca,'FontSize',fontSize,'FontWeight','Bold','LineWidth',2);
for i1=1:size(meas,1)
    plot([1,2,3],meas(i1,:),'bo-','LineWidth',2);
    hold on;
end
set(gca,'XLim',[0.5,3.5],'XTick',[1,2,3],'XTickLabel',{'noPert','Accel','Decel'});
xlabel('');
ylabel(sprintf('%s (%s)',measName,measUnit),'FontWeight','Bold');

if ~isempty(fsic(varargin,'YLim'))
    idx=fsic(varargin,'YLim');
    set(gca,'YLim',varargin{idx+1});
end

return

%% metaPlot_5
function metaPlot_5(meas,measChg,measName,measUnit,fontSize,colors,varargin)
figure('NumberTitle','off','Name',measName,'Position',[100,100,800,600]);

subplot('Position',[0.08,0.075,0.4,0.375]);
set(gca,'FontSize',fontSize,'FontWeight','Bold','LineWidth',2);

bar(1,mean(meas(:,1)),'FaceColor','w','EdgeColor','k','LineWidth',2);
hold on;
bar(2,mean(meas(:,2)),'FaceColor',colors.accel,'EdgeColor','k','LineWidth',2);
bar(3,mean(meas(:,3)),'FaceColor',colors.decel,'EdgeColor','k','LineWidth',2);
bar(4,mean(meas(:,4)),'FaceColor',colors.down,'EdgeColor','k','LineWidth',2);
bar(5,mean(meas(:,5)),'FaceColor',colors.up,'EdgeColor','k','LineWidth',2);
errorbar([1,2,3,4,5],[mean(meas(:,1)),mean(meas(:,2)),mean(meas(:,3)),mean(meas(:,4)),mean(meas(:,5))],...
    [ste(meas(:,1)),ste(meas(:,2)),ste(meas(:,3)),ste(meas(:,4)),ste(meas(:,5))],'ko','LineWidth',2);
set(gca,'XLim',[0.5,5.5],'XTick',[1,2,3,4,5],'XTickLabel',{'noPert','Accel','Decel','Down','Up'});
xlabel('');
ylabel(sprintf('%s (%s) (mean\\pmSEM)',measName,measUnit),'FontWeight','Bold','FontSize',fontSize+1);

if ~isempty(fsic(varargin,'YLim'))
    idx=fsic(varargin,'YLim');
    set(gca,'YLim',varargin{idx+1});
end
box off;
set(gca,'FontSize',fontSize,'FontWeight','Bold','LineWidth',2);

ys=get(gca,'YLim');
y_bar=ys(2)-0.06*range(ys);
y_bar_low=ys(2)-0.075*range(ys);

subplot('Position',[0.575,0.5,0.4,0.45]);
bar(1,mean(measChg(:,1)),'FaceColor',colors.accel,'EdgeColor','k','LineWidth',2);
hold on;
bar(2,mean(measChg(:,2)),'FaceColor',colors.decel,'EdgeColor','k','LineWidth',2);
bar(3,mean(measChg(:,3)),'FaceColor',colors.down,'EdgeColor','k','LineWidth',2);
bar(4,mean(measChg(:,4)),'FaceColor',colors.up,'EdgeColor','k','LineWidth',2);
errorbar([1,2,3,4],[mean(measChg(:,1)),mean(measChg(:,2)),mean(measChg(:,3)),mean(measChg(:,4))],...
    [ste(measChg(:,1)),ste(measChg(:,2)),ste(measChg(:,3)),ste(measChg(:,4))],'ko','LineWidth',2);
set(gca,'XLim',[0.5,4.5],'XTick',[1,2,3,4],'XTickLabel',{'Accel','Decel','Down','Up'});
ylabel(sprintf('Change in %s (%s)',measName,measUnit),'FontWeight','Bold','FontSize',fontSize+1);
xlabel('');
box off;
set(gca,'FontSize',fontSize,'FontWeight','Bold','LineWidth',2);

ys=get(gca,'YLim');
y_bar=ys(2)-0.05*range(ys);
y_bar_low=ys(2)-0.075*range(ys);
y_text=ys(2)-0.025*range(ys);

plot([1,2],[y_bar,y_bar],'k-','LineWidth',2);
plot([1,1],[y_bar,y_bar_low],'k-','LineWidth',2);
plot([2,2],[y_bar,y_bar_low],'k-','LineWidth',2);
plot([3,4],[y_bar,y_bar],'k-','LineWidth',2);
plot([3,3],[y_bar,y_bar_low],'k-','LineWidth',2);
plot([4,4],[y_bar,y_bar_low],'k-','LineWidth',2);
[h_chg_accelDecel,p_chg_accelDecel]=ttest2(measChg(:,1),measChg(:,2));
[h_chg_downUp,p_chg_downUp]=ttest2(measChg(:,3),measChg(:,4));
if p_chg_accelDecel<0.05
    text(1.2,y_text,sprintf('p=%.3f',p_chg_accelDecel),'FontWeight','Bold');
else
    text(1.2,y_text,sprintf('p=%.3f',p_chg_accelDecel),'FontWeight','Normal');
end
if p_chg_downUp<0.05
    text(3.2,y_text,sprintf('p=%.3f',p_chg_downUp),'FontWeight','Bold');
else
    text(3.2,y_text,sprintf('p=%.3f',p_chg_downUp),'FontWeight','Normal');
end 

subplot('Position',[0.575,0.075,0.4,0.375]);
plot(transpose(measChg),'o-','LineWidth',2); hold on;
plot([0.5,4.5],[0,0],'k-','LineWidth',1);
set(gca,'XLim',[0.5,4.5]);
set(gca,'FontSize',fontSize,'FontWeight','Bold','LineWidth',2);
set(gca,'XLim',[0.5,4.5],'XTick',[1,2,3,4],'XTickLabel',{'Accel','Decel','Down','Up'});
xlabel('');
ylabel(sprintf('Change in %s (%s)',measName,measUnit),'FontWeight','Bold','FontSize',fontSize+1);
box off; 

return



%% 
function metaPlot_3_scatter(meas,measName,measUnit,fontSize,colors,varargin)
figure('NumberTitle','off','Name',measName,'Position',[100,100,400,300]);
set(gca,'FontSize',fontSize,'FontWeight','Bold','LineWidth',2);
% for i1=1:size(meas,1)
%     plot([1,2,3],meas(i1,:),'bo-','LineWidth',2);
%     hold on;
% end

plot(meas(:,1),meas(:,2),'o','Color',colors.accel,'LineWidth',2);
hold on;
plot(meas(:,1),meas(:,3),'o','Color',colors.decel,'LineWidth',2);

xs=get(gca,'XLim'); ys=get(gca,'YLim');
lims=[min([xs(1),ys(1)]),max([ys(2),ys(2)])];
set(gca,'XLim',lims,'YLim',lims);
plot(lims,lims,'-','Color',[0.5,0.5,0.5]);

xlabel(sprintf('%s (%s) noPert',measName,measUnit),'FontWeight','Bold');
ylabel(sprintf('%s (%s) accel',measName,measUnit),'FontWeight','Bold');

return

%%


function [dacacheDir, systatDir] = get_dacacheDir(subjID, SUBJ_LIST)

flds = fields(SUBJ_LIST);
bFound = 0;
for i1 = 1 : numel(flds)
    fld = flds{i1};
    
    
    
    if ~isempty(fsic(SUBJ_LIST.(fld),subjID))
        bFound = 1;
        break;
    end
end

if bFound == 0
    dacacheDir = '';
    systatDir = '';
    return
end

studyMode = flds{i1};

if isequal(studyMode(end-3:end),'_PWS') || isequal(studyMode(end-3:end),'_PFS')
    studyMode = studyMode(1:end-4);
end

%% Path config
hostName=getHostName;
if isequal(studyMode,'apstv2')
    if isequal(hostName,'smcg-w510') || isequal(hostName,'smcgw510') || isequal(hostName,'smcg_w510')
        dacacheDir='E:/speechres/apstv2/mcode/dacache';
        systatDir='E:/speechres/apstv2/systat_files_apstv2/';
    else
        dacacheDir='Z:/speechres/apstv2/mcode/dacache';
        systatDir='Z:/speechres/apstv2/systat_files_apstv2/';
    end
elseif isequal(studyMode,'apstv2t')
    if isequal(hostName,'smcg-w510') || isequal(hostName,'smcgw510') || isequal(hostName,'smcg_w510')
        dacacheDir='E:/speechres/apstv2/mcode/dacache';
        systatDir='E:/speechres/apstv2/systat_files_apstv2t/';
    else
        dacacheDir='Z:/speechres/apstv2/mcode/dacache';
        systatDir='Z:/speechres/apstv2/systat_files_apstv2t/';
    end
elseif isequal(studyMode,'apstv')
    if isequal(hostName,'smcg-w510') || isequal(hostName,'smcgw510') || isequal(hostName,'smcg_w510')
        dacacheDir='E:/speechres/apstv/mcode/dacache';
        systatDir='E:/speechres/apstv/systat_files/';
    else
        dacacheDir='Z:/speechres/apstv/mcode/dacache';
        systatDir='Z:/speechres/apstv/systat_files/';
    end
    
elseif isequal(studyMode,'apstv2_stut_s')
    if isequal(hostName,'smcg-w510') || isequal(hostName,'smcgw510') || isequal(hostName,'smcg_w510')
        dacacheDir='E:/speechres/apstv2/mcode/dacache_stut';
        systatDir='E:/speechres/apstv2/systat_files_stut_s/';
    else
        dacacheDir='Z:/speechres/apstv/mcode/dacache_stut';
        systatDir='Z:/speechres/apstv/systat_files_stut_s/';
    end

elseif isequal(studyMode,'apstv2_stut_t')
    if isequal(hostName,'smcg-w510') || isequal(hostName,'smcgw510') || isequal(hostName,'smcg_w510')
        dacacheDir='E:/speechres/apstv2/mcode/dacache_stut';
        systatDir='E:/speechres/apstv2/systat_files_stut_t/';
    else
        dacacheDir='Z:/speechres/apstv/mcode/dacache_stut';
        systatDir='Z:/speechres/apstv/systat_files_stut_t/';
    end    
end
return