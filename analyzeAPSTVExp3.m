function varargout = analyzeAPSTVExp3(subjID,varargin)
%% Config
fontSize=12;
frameDur=16/12e3;
nInterp=250; % 250
mvaWinWidth=21;

config.normT1_bw_iu_ratio=[0,1/4,2/4,3/4,1];
config.normT1_bw_uy_ratio=[1/4,2/4,3/4,1];
config.normT2_bw_uy_ratio=[1/4,2/4,3/4,1];
config.normT3_bw_yu_ratio=[1/4,2/4,3/4,1];

bSave=0;
if ~isempty(fsic(varargin,'save')) || ~isempty(fsic(varargin, 'saveOnly'))
    bSave=1;
end

bInterpCompare = 0;
if ~isempty(fsic(varargin, 'interpCompare'))
    bInterpCompare = 1;
end

colors.noPert='k';
colors.up=[1,0,0];
colors.down=[0,0.4,1];
colors.accel=[1,0.5,0];
colors.decel=[0,0.5,0];

% colors_pert.up=0.5+0.5*colors.up;
% colors_pert.down=0.5+0.5*colors.down;
colors_pert.up='m';
colors_pert.down=[0,0.75,0.5];
colors_pert.accel=0.5*colors.accel+0.5;
colors_pert.decel=0.5*colors.decel+0.5;

lw=1;
bDiscard=1;

%%

hostName=getHostName();
if isequal(hostName,'smcg-w510') || isequal(hostName,'smcgw510') || isequal(hostName,'smcg_w510')
    rawDataDir.apstv='E:\DATA\APSTV';
    rawDataDir.apstv2='E:\DATA\APSTV2';
    rawDataDir.apstv2t='E:\DATA\APSTV2';
    rawDataDir.apstv_STUT='E:\STUT_DATA';
    rawDataDir.apstv2t_STUT='E:\STUT_DATA';
    dacacheDir.apstv='E:\speechres\apstv\mcode\dacache';
    dacacheDir.apstv2='E:\speechres\apstv2\mcode\dacache';
    dacacheDir.apstv2t='E:\speechres\apstv2\mcode\dacache';
    dacacheDir.apstv_STUT='E:\speechres\apstv2\mcode\dacache_STUT';
    dacacheDir.apstv2t_STUT='E:\speechres\apstv2\mcode\dacache_STUT';
elseif isequal(hostName,'glossa') 
    rawDataDir.apstv='G:\DATA\APSTV';
    rawDataDir.apstv2='G:\DATA\APSTV2';
    rawDataDir.apstv2t='G:\DATA\APSTV2';
    rawDataDir.apstv_STUT='G:\STUT_DATA';
    rawDataDir.apstv2t_STUT='G:\STUT_DATA';
    dacacheDir.apstv='Z:\speechres\apstv\mcode\dacache';
    dacacheDir.apstv2='Z:\speechres\apstv2\mcode\dacache';
    dacacheDir.apstv2t='Z:\speechres\apstv2\mcode\dacache';
    dacacheDir.apstv_STUT='Z:\speechres\apstv2\mcode\dacache_STUT';
    dacacheDir.apstv2t_STUT='Z:\speechres\apstv2\mcode\dacache_STUT';
elseif isequal(hostName,'Killick0')
    rawDataDir.apstv='Z:\DATA\APSTV';
    rawDataDir.apstv2='Z:\DATA\APSTV2';
    rawDataDir.apstv2t='Z:\DATA\APSTV2';
    rawDataDir.apstv_STUT='Z:\STUT_DATA';
    rawDataDir.apstv2t_STUT='Z:\STUT_DATA';
    dacacheDir.apstv='Z:\speechres\apstv\mcode\dacache';
    dacacheDir.apstv2='Z:\speechres\apstv2\mcode\dacache';
    dacacheDir.apstv2t='Z:\speechres\apstv2\mcode\dacache';
    dacacheDir.apstv_STUT='Z:\speechres\apstv2\mcode\dacache_STUT';
    dacacheDir.apstv2t_STUT='Z:\speechres\apstv2\mcode\dacache_STUT';
elseif isequal(lower(hostName),'northroom_745')
    rawDataDir.apstv='Z:\DATA\APSTV';    % TODO
    rawDataDir.apstv2='Z:\DATA\APSTV2';   % TODO
    rawDataDir.apstv2t='Z:\DATA\APSTV2';   % TODO
    rawDataDir.apstv_STUT='Z:\STUT_DATA';
    rawDataDir.apstv2t_STUT='Z:\STUT_DATA';
    dacacheDir.apstv='Z:\speechres\apstv\mcode\dacache';
    dacacheDir.apstv2='Z:\speechres\apstv2\mcode\dacache';
    dacacheDir.apstv2t='Z:\speechres\apstv2\mcode\dacache';
    dacacheDir.apstv_STUT='Z:\speechres\apstv2\mcode\dacache_STUT';
    dacacheDir.apstv2t_STUT='Z:\speechres\apstv2\mcode\dacache_STUT';
end


%% Decide whether the subject is a APSTV or APSTV2 subject
if isequal(subjID(1:5),'AS_PS') || isequal(subjID(1:6),'APSTV_')
    et='apstv'; % et stands for "experiment type"
    expDir=fullfile(rawDataDir.(et),subjID);
    dacacheFN=fullfile(dacacheDir.(et),[subjID,'.mat']);
elseif isequal(subjID(1:7),'APSTV2_') || isequal(subjID(1:3),'IOU')
    et='apstv2';
    expDir=fullfile(rawDataDir.(et),subjID);
    dacacheFN=fullfile(dacacheDir.(et),[subjID,'.mat']);
elseif length(subjID)>=8 && isequal(subjID(1:8),'APSTV2T_') 
    et='apstv2t';
    expDir=fullfile(rawDataDir.(et),subjID);
    dacacheFN=fullfile(dacacheDir.(et),[subjID,'.mat']);
else
    if ~isempty(fsic(varargin,'STUT_S'))
        et='apstv_STUT';
        expDir=fullfile(rawDataDir.(et),subjID,'APSTV2_STUT_S');
        dacacheFN=fullfile(dacacheDir.(et),[subjID,'_S.mat']);
    elseif ~isempty(fsic(varargin,'STUT_T'))
        et='apstv2t_STUT';
        expDir=fullfile(rawDataDir.(et),subjID,'APSTV2_STUT_T');
        dacacheFN=fullfile(dacacheDir.(et),[subjID,'_T.mat']);
    else
        et='unknown';
        fprintf('ERROR: unknown experiment type (et).\n');
        return
    end        
end


if ~isdir(expDir)
    fprintf('ERROR: Cannot find directory %s. Terminated.\n',expDir);
    return
end
if ~isfile(fullfile(expDir,'expt.mat'))
    fprintf('ERROR: Cannot find expt.mat in directory %s. Terminated.\n',expDir);
    return
end


if ~isfile(dacacheFN)
    fprintf('ERROR: Cannot find pdata file %s.\n',dacacheFN);
    return
end

load(fullfile(expDir,'expt.mat'));  % gives expt
load(dacacheFN); % gives pdata
fprintf('Subject ID: \t%s\n',expt.subject.name);
fprintf('Subject sex: \t%s\n',expt.subject.sex);
fprintf('Experiment type (et): %s\n',et);

if isequal(et,'apstv') || isequal(et,'apstv_STUT')
    pertTypes={'noPert','down','up'};
    idx_noPert=fsic(pertTypes,'noPert');
    fprintf('\tUtter perts: \t[Down: %.2f, Up: %.2f]\n',abs(min(expt.utterPert.main)),max(expt.utterPert.main));
elseif isequal(et,'apstv2') || isequal(et,'apstv2t') || isequal(et,'apstv2t_STUT')
    pertTypes=unique(expt.utterPertType.main);
    fprintf('\tUtter perts: \n');
    for i1=1:numel(pertTypes)
        if isequal(pertTypes{i1},'none')            
            pertTypes{i1}='noPert';
            continue;
        else
            idx1=fsic(expt.utterPertType.main,pertTypes{i1});
            fprintf('\t%s: \t%.2f\n',pertTypes{i1},expt.utterPert.main(idx1));
        end
            
    end
    idx_noPert=fsic(pertTypes,'noPert');
    pertTypes=[pertTypes{idx_noPert},pertTypes(setxor(1:numel(pertTypes),idx_noPert))];
end

if ~isequal(expt.subject.name,subjID)
    fprintf('WARNING: subject ID mismatch between input argument and the expt.subject field "name"\n')
end
if ~isequal(pdata.subject.name,subjID)
    fprintf('WARNING: subject ID mismatch between input argument and the pdata.subject field "name"\n')
end

%% Obtain information about the perturbation field
d0=dir(fullfile(expDir,'main','rep1','trial-*-1.mat'));
if isequal(et,'apstv') || isequal(et,'apstv_STUT')
    for i1=1:length(d0)
        load(fullfile(expDir,'main','rep1',d0(i1).name));    % gives data
        if data.params.bShift==1 && (~isfield(data.params,'bTimeWarp') || data.params.bTimeWarp==0)
            pdata.pertField=data.params.pertField;
            break;
        end
    end
    if (~isfield(pdata,'pertField') || ~isfield(pdata.pertField,'F2UB'))
        fprintf('ERROR: failed to obatain the necessary perturbation-field information. Terminated.\n');
        return;
    end
else
    pdata.pertField=expt.pertInfo;
end


%%
protoStruct_cell=struct;
for i1=1:numel(pertTypes)
    protoStruct_cell.(pertTypes{i1})={};
end
protoStruct_array=struct;
for i1=1:numel(pertTypes)
    protoStruct_array.(pertTypes{i1})=[];
end

protoStruct_array_AT=struct;
flds0=fields(protoStruct_cell);
for i1=1:numel(flds0)
    protoStruct_array_AT.(['aft_',flds0{i1}])=[];
end

if ~isdir(fullfile(expDir,'main'))
	fprintf('"Main" sub-directory doesn''t exist in %s. Terminated.\n',expDir);
	return
end

d1=dir(fullfile(expDir,'main','rep*'));
nReps=length(d1);

iuInt=protoStruct_array;
iyInt=protoStruct_array;
uyInt=protoStruct_array;
iu2Int=protoStruct_array;

i_uyMid_int=protoStruct_array;
i_yu2Mid_int=protoStruct_array;
i_uyMidF_int=protoStruct_array;
i_yu2MidF_int=protoStruct_array;

a_rawFNs=protoStruct_cell;

iF2=protoStruct_array;
yF2=protoStruct_array;
uF2=protoStruct_array;
yF2_you=protoStruct_array;
uF2_you=protoStruct_array;
iuF2Dec=protoStruct_array;
uyF2Inc=protoStruct_array;
youF2Dec=protoStruct_array;
iuF2Vel=protoStruct_array;
uyF2Vel=protoStruct_array;
f1Trajs=protoStruct_cell;
f2Trajs=protoStruct_cell;
f2Trajs_pert=protoStruct_cell;
f1Trajs_pert=protoStruct_cell;
f2_pertShift=protoStruct_cell;
f1_pertShift=protoStruct_cell;


rmsTrajs=protoStruct_cell;

% --- New initialization --- %
fmtFields = {'F1', 'F2'};
segTraj = struct; % Linear inerpolation (default)
if bInterpCompare
    segTraj_nn = struct; % Nearest neighbor interpolation
    segTraj_cs = struct; % Cubic spline interpolation
end

N_segs = 6;
% 1 - tNorm1a: between [i] and [u]_1
% 2 - tNorm2: between [u]_1 and [j]_1
% 3 - tNorm3: between [j]_1 and [u]_2
% 4 - tNorm4: between [u]_2 and [j]_2
% 5 - tNorm5: between [j]_2 and [u]_3
% 6 - tNorm6: between [u]_3 and [j]_3

for i1 = 1 : numel(fmtFields)
    ff = fmtFields{i1};
    
    segTraj.(ff) = cell(1, N_segs);
    if bInterpCompare
        segTraj_nn.(ff) = cell(1, N_segs);
        segTraj_cs.(ff) = cell(1, N_segs);
    end

    for i2 = 1 : N_segs
        segTraj.(ff){i2} = protoStruct_cell;
        
        if bInterpCompare
            segTraj_nn.(ff){i2} = protoStruct_cell;
            segTraj_cs.(ff){i2} = protoStruct_cell;
        end
    end
end 

% --- Old initializations --- %
segTrajF2_tNorm1=protoStruct_cell;   % The time normalized F2 trajectories between iouOnset and youOnset. Anchorage between iouOnset and youOnset
% segTrajF2_tNorm2=protoStruct_cell;   % The time normalized F2 trajectories between iouOnset and youOnset. Anchorage between uTime and youOnset
segTrajF2_tNorm2_pert=protoStruct_cell;

% segTrajF2_tNorm2_nn=protoStruct_cell;
% segTrajF2_tNorm2_cs=protoStruct_cell;

segTrajF2_tNorm1a=protoStruct_cell;
segTrajF2_tNorm1a_pert=protoStruct_cell;

segTrajF2_tNorm1a_nn = protoStruct_cell;    % Nearest neighbor
segTrajF2_tNorm1a_cs = protoStruct_cell;    % Cubic spline

segVelF2_tNorm1=protoStruct_cell;

% segTrajF2_tNorm3=protoStruct_cell;	% Time normalized F2 trajectories between youOnset and uayOnset
% segTrajF2_tNorm4=protoStruct_cell;	% Time normalized F2 trajectories between uayOnset and yo1Onset
% segTrajF2_tNorm5=protoStruct_cell;  % Time normalized F2 trajectories between yo1Onset and u2Time
% segTrajF2_tNorm6=protoStruct_cell;  % Time normalized F2 trajectories between u2Time and yo2Onset

% segTrajF2_tNorm3_nn=protoStruct_cell;
% segTrajF2_tNorm4_nn=protoStruct_cell;
% segTrajF2_tNorm5_nn=protoStruct_cell;
% segTrajF2_tNorm6_nn=protoStruct_cell;

% segTrajF2_tNorm3_cs=protoStruct_cell;
% segTrajF2_tNorm4_cs=protoStruct_cell;
% segTrajF2_tNorm5_cs=protoStruct_cell;
% segTrajF2_tNorm6_cs=protoStruct_cell;

% --- ~Old initializations --- %

segTrajF2_tNorm3_pert=protoStruct_cell;	% Time normalized F2 trajectories between youOnset and uayOnset
segTrajF2_tNorm4_pert=protoStruct_cell;	% Time normalized F2 trajectories between uayOnset and yo1Onset
segTrajF2_tNorm5_pert=protoStruct_cell;  % Time normalized F2 trajectories between yo1Onset and u2Time
segTrajF2_tNorm6_pert=protoStruct_cell;  % Time normalized F2 trajectories between u2Time and yo2Onset

segPertShiftF2_tNorm1=protoStruct_cell;
segPertShiftF2_tNorm2=protoStruct_cell;
segPertShiftF2_tNorm3=protoStruct_cell;
segPertShiftF2_tNorm4=protoStruct_cell;
segPertShiftF2_tNorm5=protoStruct_cell;
segPertShiftF2_tNorm6=protoStruct_cell;

spect_iu2=protoStruct_cell;

trajF1_POA=protoStruct_cell;          % POA stands for perturbation-onset aligned.
trajF2_POA=protoStruct_cell;
traj_POA_fn=protoStruct_cell;

normT1_bw_iu_f2=cell(1,0);
for i1=1:length(config.normT1_bw_iu_ratio)
    normT1_bw_iu_f2{i1}=protoStruct_array;
end
normT1_bw_uy_f2=cell(1,0);
for i1=1:length(config.normT1_bw_uy_ratio)
    normT1_bw_uy_f2{i1}=protoStruct_array;
end
normT2_bw_uy_f2=cell(1,0);
for i1=1:length(config.normT2_bw_uy_ratio)
    normT2_bw_uy_f2{i1}=protoStruct_array;
end
normT3_bw_yu_f2=cell(1,0);
for i1=1:length(config.normT3_bw_yu_ratio)
    normT3_bw_yu_f2{i1}=protoStruct_array;
end

timeShifts.minF2.accel=[];
timeShifts.minF2.decel=[];
timeShifts.iuMidF.accel=[];
timeShifts.iuMidF.decel=[];

prevPertTypes=protoStruct_cell;

%% Data extraction
for i1=1:numel(pertTypes)
    fld=pertTypes{i1};
    fld0=fld;
    if isequal(et,'apstv2') || isequal(et,'apstv2t') || isequal(et,'apstv_STUT') || isequal(et,'apstv2t_STUT')
        if isequal(fld,'noPert')
            fld='none';
        end
    else
       if isequal(fld,'down')
           fld='accel';
       elseif isequal(fld,'up')
           fld='decel';
       end 
    end
    for i2=1:numel(pdata.utters.(fld))        
        this_utter=pdata.utters.(fld){i2};
        
        if isempty(this_utter)
            fprintf('WARNING: trial %s #%d is empty. Discarded.\n',fld0,i2);
            continue;
        end
        if ~isfield(this_utter,'bDiscard')
            fprintf('WARNING: trial %s #%d has no field "bDiscard". Discarded.\n',fld0,i2);
            continue;
        end
        
        if bDiscard==1
            if this_utter.bDiscard==1
                fprintf('%s : %s discarded\n',this_utter.pertStr,this_utter.rawDataFN);
                continue;
            end
            if (isfield(this_utter,'rating') && this_utter.rating==0)
                fprintf('%s : %s discard due to rating==0\n',this_utter.pertStr,this_utter.rawDataFN);
                continue;
            end
        end
        
        if isfield(this_utter,'bPertOkay')
            if this_utter.bPertOkay==0 || isequal(this_utter.bPertOkay,'0')
                fprintf('%s : %s discarded due to bPertOkay == 0\n',this_utter.pertStr,this_utter.rawDataFN);
                continue;
            end
        end
        
        iouOnset=this_utter.iouOnset;
        youOnset=this_utter.youOnset;
        uTime=this_utter.uTime;
        u2Time=this_utter.uayOnset;
        
        if iouOnset>=uTime || uTime>=youOnset || iouOnset>=youOnset
            fprintf('WARNING: trial %s contains wrong order of the events {iouOnset, uTime and youOnset}. Discarded.\n',this_utter.rawDataFN);
            continue;
        end
        
        t_iF2=this_utter.iF2;
        t_uF2=this_utter.uF2;
        t_yF2=this_utter.yF2;
        t_uF2_you=this_utter.uF2_you;
        
        if t_uF2>t_iF2 || t_uF2>t_yF2
            fprintf('WARNING: trial %s contains uF2 greater than iF2 or yF2. Discarded.\n',this_utter.rawDataFN);
            continue;
        end
        
        pertStr=this_utter.pertStr;
        if ~isequal(pertStr,fld)
            fprintf('WARNING: trial %s contains pertStr that does not match its field location in pdata.utters.\n',this_utter.rawDataFN);
            continue;
        end
        
        if ~isfield(this_utter,'traj_F2') 
            fprintf('WARNING: trial %s does not contain the "traj_F2" field.\n',this_utter.rawDataFN);
            continue;
        end
        
        if isempty(this_utter.traj_F2)
            fprintf('WARNING: trial %s contains an empty traj_F2 field.\n',this_utter.rawDataFN);
            continue;
        end
        
        raw_data_fn=getRawFN_(expDir,this_utter.rawDataFN);                        
        
        if ~(this_utter.uTime>this_utter.iouOnset && this_utter.youOnset>this_utter.uTime && ...
             this_utter.uayOnset>this_utter.youOnset && this_utter.yo1Onset>this_utter.uayOnset && ...
             this_utter.yo1End>this_utter.yo1Onset && this_utter.yo2Onset>this_utter.yo1End)
            fprintf('Temporal landmarks order problematic: %s\n',this_utter.rawDataFN);           
        end
        
        if (~isempty(iouOnset) && ~isempty(youOnset) && ~isempty(uTime) && ~isnan(iouOnset) && ~isnan(youOnset) && ~isempty(uTime))
			iuInt.(fld0)(end+1)=uTime-iouOnset;
			iyInt.(fld0)(end+1)=youOnset-iouOnset;
			uyInt.(fld0)(end+1)=youOnset-uTime;
            iu2Int.(fld0)(end+1)=u2Time-iouOnset;
            
%             if iuInt.(fld0)(end) < 0.15
%                 pause(0);   % DEBUG
%             end
            
            
            i_uyMid_int.(fld0)(end+1)=(uTime+youOnset)/2-iouOnset;
            i_yu2Mid_int.(fld0)(end+1)=(youOnset+u2Time)/2-iouOnset;
            
            taxis1=iouOnset:frameDur:iouOnset+frameDur*(length(this_utter.traj_F2)-1);
            f2Seg=this_utter.traj_F2(taxis1>=uTime & taxis1<youOnset);                        
                                    
            if ~isempty(find(this_utter.traj_F2==0))
                pause(0);
            end
            
            f2Mid=mean([t_uF2,t_yF2]);
            i_uyMidF_int.(fld0)(end+1)=get_midF(f2Seg,f2Mid,uTime,frameDur,'up')-iouOnset;
            
            f2Seg=this_utter.traj_F2(taxis1>=youOnset & taxis1<u2Time);
            f2Mid=mean([t_yF2,t_uF2_you]);
            i_yu2MidF_int.(fld0)(end+1)=get_midF(f2Seg,f2Mid,youOnset,frameDur,'down')-iouOnset;
            
			iF2.(fld0)(end+1)=t_iF2;
			yF2.(fld0)(end+1)=t_yF2;
			uF2.(fld0)(end+1)=t_uF2;            
            uF2_you.(fld0)(end+1)=t_uF2_you;
			iuF2Dec.(fld0)(end+1)=t_iF2-t_uF2;
			uyF2Inc.(fld0)(end+1)=t_yF2-t_uF2;
            youF2Dec.(fld0)(end+1)=t_yF2-t_uF2_you;
            
			iuF2Vel.(fld0)(end+1)=(t_uF2-t_iF2)/(uTime-iouOnset);
			uyF2Vel.(fld0)(end+1)=(t_yF2-t_uF2)/(youOnset-uTime);
            
            a_rawFNs.(fld0){end+1}=raw_data_fn;
            
            
            
        end
        
        if ~isempty(find(this_utter.traj_F2 == 0))
            error('Zero values in F2 trajectory: rawDataFN = %s', this_utter.rawDataFN);
        end
        
        f1Trajs.(fld0){end+1}=this_utter.traj_F1;
        f2Trajs.(fld0){end+1}=this_utter.traj_F2;
        
        if ~isempty(find(this_utter.traj_F2==0))
%             pause(0);
        end
         
        load(raw_data_fn);  % gives data
        taxis0=0:(1/data.params.sr):(1/data.params.sr*(length(data.signalIn)-1));
        sigSeg=data.signalIn(taxis0>this_utter.iouOnset & taxis0<this_utter.uayOnset);
        [s,f,t]=spectrogram(sigSeg,128,96,1024,data.params.sr);
        s=s(f<3000,:);
        
        spect_iu2.(fld0){end+1}=10*log10(abs(s));
%         plot(t); hold on;

        if ~isequal(fld0,'noPert')
            pertType_data=getPertType_data(data);
            if isequal(et,'apstv')
                if isequal(pertType_data,'down')
                    pertType_data='accel';
                elseif isequal(pertType_data,'up')
                    pertType_data='decel';
                end
            end            
            
            if ~isequal(fld,pertType_data)
                fprintf('WARNING: mismatch between fld0 and pertType_data in %s.\n',this_utter.rawDataFN);
            end
            
            frameDur=data.params.frameLen/data.params.sr;
            taxis1p=0:frameDur:(frameDur*(size(data.fmts(:,2)-1)));
            of1=data.fmts(:,1);
            of2=data.fmts(:,2);
            sf1=data.sfmts(:,1);
            sf2=data.sfmts(:,2);
            
            ssf1=sf1;
            ssf1(ssf1==0)=of1(ssf1==0);
            ssf2=sf2;
            ssf2(ssf2==0)=of2(ssf2==0);
            
            of1=mva_nz(of1,mvaWinWidth,'Hamming');
            of2=mva_nz(of2,mvaWinWidth,'Hamming');
            sf1=mva_nz(sf1,mvaWinWidth,'Hamming');
            sf2=mva_nz(sf2,mvaWinWidth,'Hamming');
            ssf1=mva_nz(ssf1,mvaWinWidth,'Hamming');
            ssf2=mva_nz(ssf2,mvaWinWidth,'Hamming');
            
%             orig_f2_a=data.fmts(:,2);
%             orig_f2_a=mva_nz(orig_f2_a,mvaWinWidth,'Hamming');
%             orig_f2_a=orig_f2_a(taxis1p>=this_utter.iouOnset & taxis1p<=this_utter.yo2Onset);
%             orig_f1_a=data.fmts(:,1);
%             orig_f1_a=mva_nz(orig_f1_a,mvaWinWidth,'Hamming');
%             orig_f1_a=orig_f1_a(taxis1p>=this_utter.iouOnset & taxis1p<=this_utter.yo2Onset);
            orig_f2=this_utter.traj_F2;
            orig_f1=this_utter.traj_F1;
            
            f2Trajs_pert.(fld0){end+1}=sf2(taxis1p>=this_utter.iouOnset & taxis1p<=this_utter.yo2Onset);
            f2Trajs_pert.(fld0){end}(f2Trajs_pert.(fld0){end}==0)=orig_f2(f2Trajs_pert.(fld0){end}==0);
            f1Trajs_pert.(fld0){end+1}=sf1(taxis1p>=this_utter.iouOnset & taxis1p<=this_utter.yo2Onset);
            f1Trajs_pert.(fld0){end}(f1Trajs_pert.(fld0){end}==0)=orig_f1(f1Trajs_pert.(fld0){end}==0);
            
            minLen=min([length(f2Trajs_pert.(fld0){end}),length(orig_f2)]);
            f2_pertShift.(fld0){end+1}=f2Trajs_pert.(fld0){end}(1:minLen)-orig_f2(1:minLen);
            minLen=min([length(f1Trajs_pert.(fld0){end}),length(orig_f1)]);
            f1_pertShift.(fld0){end+1}=f1Trajs_pert.(fld0){end}(1:minLen)-orig_f1(1:minLen);
        end
        
        if isequal(fld0,'noPert')
            % For AT (adaptation analysis): figure out the pert type of the previous trial
            prevPertTypes.(fld0){end+1}=get_prev_pert_type(this_utter.repNum,this_utter.trialNum,expt); 
            if isequal(prevPertTypes.(fld0){end},'none')
                prevPertTypes.(fld0){end}='noPert';
            end
        end
        
        taxis1 = 0 : frameDur : frameDur * (length(this_utter.traj_F2) - 1);
        
        % Time-normalized F1 and F2 trajectory: anchorage 1: between [i] and [j]_1 (not used in generating
        % the complete time-normalized trajectory
        t_segTrajF2 = this_utter.traj_F2(taxis1 <= this_utter.youOnset - this_utter.iouOnset);
        t_segTrajF2 = interp1(1 : length(t_segTrajF2), t_segTrajF2, linspace(1, length(t_segTrajF2), nInterp));
        segTrajF2_tNorm1.(fld0){end + 1}=t_segTrajF2;
        
        assert(isequal(size(this_utter.traj_F1), size(this_utter.traj_F2)));
        for i3 = 1 : 6
            if i3 == 1  % tNorm1a: between [i] and [u]_1
                t_beg = this_utter.iouOnset; 
                t_end = this_utter.uTime;
            elseif i3 == 2 % tNorm2: between [u]_1 and [j]_1
                t_beg = this_utter.uTime; 
                t_end = this_utter.youOnset;
            elseif i3 == 3 % tNorm3: between [j]_1 and [u]_2
                t_beg = this_utter.youOnset;
                t_end = this_utter.uayOnset;
            elseif i3 == 4 % tNorm4: between [u]_2 and [j]_2
                t_beg = this_utter.uayOnset;
                t_end = this_utter.yo1Onset;
            elseif i3 == 5 % tNorm5: between [j]_2 and [u]_3
                t_beg = this_utter.yo1Onset;
                t_end = this_utter.yo1End;
            elseif i3 == 6 % tNorm6: between [u]_3 and [j]_3
                t_beg = this_utter.yo1End;
                t_end = this_utter.yo2Onset;
            end
            assert(t_end > t_beg);
                        
            for i4 = 1 : numel(fmtFields)
                ff = fmtFields{i4};
                tff = sprintf('traj_%s', ff);
            
                t_segTraj = this_utter.(tff)(taxis1 >= t_beg - this_utter.iouOnset & ...
                                               taxis1 <= t_end - this_utter.iouOnset);
                if ~isempty(t_segTraj)
                    segTraj.(ff){i3}.(fld0){end + 1} = ...
                        interp1(1 : length(t_segTraj), t_segTraj, linspace(1,length(t_segTraj), nInterp));
                    
                    if bInterpCompare
                        segTraj_nn.(ff){i3}.(fld0){end + 1} = ...
                            interp1(1 : length(t_segTraj), t_segTraj, linspace(1,length(t_segTraj), nInterp), 'nearest');
                        segTraj_cs.(ff){i3}.(fld0){end + 1} = ...
                            interp1(1 : length(t_segTraj), t_segTraj, linspace(1,length(t_segTraj), nInterp), 'spline');
                    end
                else
                    segTraj.(ff){i3}.(fld0){end + 1} = nan(1, nInterp);
                    
                    if bInterpCompare
                        segTraj_nn.(ff){i3}.(fld0){end + 1} = nan(1, nInterp);
                        segTraj_cs.(ff){i3}.(fld0){end + 1} = nan(1, nInterp);
                    end
                end
            end
        end
        
        % Anchorage 1a: between [i] and [u]_1 (used in genereating complete time-normalized trajectory)
%         t_segTrajF2=this_utter.traj_F2(taxis1<=this_utter.uTime-this_utter.iouOnset);
%         if ~isempty(t_segTrajF2)
%             segTrajF2_tNorm1a.(fld0){end+1}= ... % t_segTrajF2;
%                 interp1(1:length(t_segTrajF2),t_segTrajF2,linspace(1,length(t_segTrajF2),nInterp));            
%             segTrajF2_tNorm1a_nn.(fld0){end+1}= ... % t_segTrajF2;
%                 interp1(1:length(t_segTrajF2),t_segTrajF2,linspace(1,length(t_segTrajF2),nInterp), 'nearest');
%             segTrajF2_tNorm1a_cs.(fld0){end+1}= ... % t_segTrajF2;
%                 interp1(1:length(t_segTrajF2),t_segTrajF2,linspace(1,length(t_segTrajF2),nInterp), 'spline');            
%         else
%             segTrajF2_tNorm1a.(fld0){end+1}=nan(1,nInterp);
%             segTrajF2_tNorm1a_nn.(fld0){end+1}=nan(1,nInterp);
%             segTrajF2_tNorm1a_cs.(fld0){end+1}=nan(1,nInterp);
%         end
        if ~isequal(fld0,'noPert')
            t_segTrajF2_pert=f2Trajs_pert.(fld0){end}(taxis1<=this_utter.uTime-this_utter.iouOnset);
            t_segTrajF2_pert=interp1(1:length(t_segTrajF2_pert),t_segTrajF2_pert,linspace(1,length(t_segTrajF2_pert),nInterp));
            segTrajF2_tNorm1a_pert.(fld0){end+1}=t_segTrajF2_pert;
        end

        % Time-normalized F2 trajectory: anchorage 2: between [u]_1 and [j]_1
%         t_segTrajF2=this_utter.traj_F2(taxis1>=this_utter.uTime-this_utter.iouOnset & taxis1<=this_utter.youOnset-this_utter.iouOnset);
%         if ~isempty(t_segTrajF2)
%             segTrajF2_tNorm2.(fld0){end+1} = ...
%                 interp1(1:length(t_segTrajF2),t_segTrajF2,linspace(1,length(t_segTrajF2),nInterp));
%             segTrajF2_tNorm2_nn.(fld0){end+1} = ...
%                 interp1(1:length(t_segTrajF2),t_segTrajF2,linspace(1,length(t_segTrajF2),nInterp), 'nearest');
%             segTrajF2_tNorm2_cs.(fld0){end+1} = ...
%                 interp1(1:length(t_segTrajF2),t_segTrajF2,linspace(1,length(t_segTrajF2),nInterp), 'spline');
%         else
%             segTrajF2_tNorm2.(fld0){end+1}=nan(1,nInterp);
%             segTrajF2_tNorm2_nn.(fld0){end+1}=nan(1,nInterp);
%             segTrajF2_tNorm2_cs.(fld0){end+1}=nan(1,nInterp);
%         end
        if ~isequal(fld0,'noPert')
            t_segTrajF2_pert=f2Trajs_pert.(fld0){end}(taxis1>=this_utter.uTime-this_utter.iouOnset & taxis1<=this_utter.youOnset-this_utter.iouOnset);
            t_segTrajF2_pert=interp1(1:length(t_segTrajF2_pert),t_segTrajF2_pert,linspace(1,length(t_segTrajF2_pert),nInterp));
            segTrajF2_tNorm2_pert.(fld0){end+1}=t_segTrajF2_pert;
        end
        
        % Time-normalized F1 trajectory: anchorage 3: between [j]_1 and [u]_2
%         t_segTrajF2=this_utter.traj_F2(taxis1>=this_utter.youOnset-this_utter.iouOnset & taxis1<=this_utter.uayOnset-this_utter.iouOnset);
%         if ~isempty(t_segTrajF2)            
%             segTrajF2_tNorm3.(fld0){end+1} = ...
%                 interp1(1:length(t_segTrajF2),t_segTrajF2,linspace(1,length(t_segTrajF2),nInterp));
%             segTrajF2_tNorm3_nn.(fld0){end+1} = ...
%                 interp1(1:length(t_segTrajF2),t_segTrajF2,linspace(1,length(t_segTrajF2),nInterp), 'nearest');
%             segTrajF2_tNorm3_cs.(fld0){end+1} = ...
%                 interp1(1:length(t_segTrajF2),t_segTrajF2,linspace(1,length(t_segTrajF2),nInterp), 'spline');
%         else
%             segTrajF2_tNorm3.(fld0){end+1}=nan(1,nInterp);
%             segTrajF2_tNorm3_nn.(fld0){end+1}=nan(1,nInterp);
%             segTrajF2_tNorm3_cs.(fld0){end+1}=nan(1,nInterp);
%         end   
        if ~isequal(fld0,'noPert')
            t_segTrajF2_pert=f2Trajs_pert.(fld0){end}(taxis1>=this_utter.youOnset-this_utter.iouOnset & taxis1<=this_utter.uayOnset-this_utter.iouOnset);
            t_segTrajF2_pert=interp1(1:length(t_segTrajF2_pert),t_segTrajF2_pert,linspace(1,length(t_segTrajF2_pert),nInterp));
            segTrajF2_tNorm3_pert.(fld0){end+1}=t_segTrajF2_pert;
        end
        
        % Time normalied F2 trajectory: anchorage 4: between [u]_2 and [j]_2
%         t_segTrajF2=this_utter.traj_F2(taxis1>=this_utter.uayOnset-this_utter.iouOnset & taxis1<=this_utter.yo1Onset-this_utter.iouOnset);
%         if ~isempty(t_segTrajF2)            
%             segTrajF2_tNorm4.(fld0){end+1} = ...
%                 interp1(1:length(t_segTrajF2),t_segTrajF2,linspace(1,length(t_segTrajF2),nInterp));
%             segTrajF2_tNorm4_nn.(fld0){end+1} = ...
%                 interp1(1:length(t_segTrajF2),t_segTrajF2,linspace(1,length(t_segTrajF2),nInterp), 'nearest');
%             segTrajF2_tNorm4_cs.(fld0){end+1} = ...
%                 interp1(1:length(t_segTrajF2),t_segTrajF2,linspace(1,length(t_segTrajF2),nInterp), 'spline');
%         else            
%             segTrajF2_tNorm4.(fld0){end+1}=nan(1,nInterp);
%             segTrajF2_tNorm4_nn.(fld0){end+1}=nan(1,nInterp);
%             segTrajF2_tNorm4_cs.(fld0){end+1}=nan(1,nInterp);
%         end
        if ~isequal(fld0,'noPert')
            t_segTrajF2_pert=f2Trajs_pert.(fld0){end}(taxis1>=this_utter.uayOnset-this_utter.iouOnset & taxis1<=this_utter.yo1Onset-this_utter.iouOnset);
            t_segTrajF2_pert=interp1(1:length(t_segTrajF2_pert),t_segTrajF2_pert,linspace(1,length(t_segTrajF2_pert),nInterp));
            segTrajF2_tNorm4_pert.(fld0){end+1}=t_segTrajF2_pert;
        end
        
        % Time normalied F2 trajectory: anchorage 5: between [j]_2 and [u]_3
        t_segTrajF2=this_utter.traj_F2(taxis1>=this_utter.yo1Onset-this_utter.iouOnset & taxis1<=this_utter.yo1End-this_utter.iouOnset);
%         if ~isempty(t_segTrajF2)            
%             segTrajF2_tNorm5.(fld0){end+1} = ...
%                 interp1(1:length(t_segTrajF2),t_segTrajF2,linspace(1,length(t_segTrajF2),nInterp));
%             segTrajF2_tNorm5_nn.(fld0){end+1} = ...
%                 interp1(1:length(t_segTrajF2),t_segTrajF2,linspace(1,length(t_segTrajF2),nInterp), 'nearest');
%             segTrajF2_tNorm5_cs.(fld0){end+1} = ...
%                 interp1(1:length(t_segTrajF2),t_segTrajF2,linspace(1,length(t_segTrajF2),nInterp), 'spline');
%         else
%             segTrajF2_tNorm5.(fld0){end+1}=nan(1,nInterp);
%             segTrajF2_tNorm5_nn.(fld0){end+1}=nan(1,nInterp);
%             segTrajF2_tNorm5_cs.(fld0){end+1}=nan(1,nInterp);
%         end
        if ~isequal(fld0,'noPert')
            t_segTrajF2_pert=f2Trajs_pert.(fld0){end}(taxis1>=this_utter.yo1Onset-this_utter.iouOnset & taxis1<=this_utter.yo1End-this_utter.iouOnset);
            t_segTrajF2_pert=interp1(1:length(t_segTrajF2_pert),t_segTrajF2_pert,linspace(1,length(t_segTrajF2_pert),nInterp));
            segTrajF2_tNorm5_pert.(fld0){end+1}=t_segTrajF2_pert;
        end
        
        % Time normalied F2 trajectory: anchorage 6: between [u]_3 and [j]_3
        t_segTrajF2=this_utter.traj_F2(taxis1>=this_utter.yo1End-this_utter.iouOnset & taxis1<=this_utter.yo2Onset-this_utter.iouOnset);
%         if ~isempty(t_segTrajF2)            
%             segTrajF2_tNorm6.(fld0){end+1} = ...
%                 interp1(1:length(t_segTrajF2),t_segTrajF2,linspace(1,length(t_segTrajF2),nInterp));
%             segTrajF2_tNorm6_nn.(fld0){end+1} = ...
%                 interp1(1:length(t_segTrajF2),t_segTrajF2,linspace(1,length(t_segTrajF2),nInterp), 'nearest');
%             segTrajF2_tNorm6_cs.(fld0){end+1} = ...
%                 interp1(1:length(t_segTrajF2),t_segTrajF2,linspace(1,length(t_segTrajF2),nInterp), 'spline');
%         else
%             segTrajF2_tNorm6.(fld0){end+1}=nan(1,nInterp);
%             segTrajF2_tNorm6_nn.(fld0){end+1}=nan(1,nInterp);
%             segTrajF2_tNorm6_cs.(fld0){end+1}=nan(1,nInterp);
%         end
        if ~isequal(fld0,'noPert')
            idx_o=find(taxis1>=this_utter.yo1End-this_utter.iouOnset & taxis1<=this_utter.yo2Onset-this_utter.iouOnset);
            t_segTrajF2_pert=f2Trajs_pert.(fld0){end}(idx_o(1):min([idx_o(end),length(f2Trajs_pert.(fld0){end})]));
            t_segTrajF2_pert=interp1(1:length(t_segTrajF2_pert),t_segTrajF2_pert,linspace(1,length(t_segTrajF2_pert),nInterp));
            segTrajF2_tNorm6_pert.(fld0){end+1}=t_segTrajF2_pert;
        end
        
        if ~(isequal(fld0,'noPert') || isequal(fld,'none'))
            t_segPertShiftF2=ssf2(taxis1p>=this_utter.iouOnset & taxis1p<=this_utter.uTime)-of2(taxis1p>=this_utter.iouOnset & taxis1p<=this_utter.uTime);
            if~isempty(t_segPertShiftF2)
                t_segPertShiftF2=interp1(1:length(t_segPertShiftF2),t_segPertShiftF2,linspace(1,length(t_segPertShiftF2),nInterp));
                segPertShiftF2_tNorm1.(fld0){end+1}=t_segPertShiftF2;
            else
                segPertShiftF2_tNorm1.(fld0){end+1}=nan(1,nInterp);
            end
            
            t_segPertShiftF2=ssf2(taxis1p>=this_utter.uTime & taxis1p<=this_utter.youOnset)-of2(taxis1p>=this_utter.uTime & taxis1p<=this_utter.youOnset);
            if~isempty(t_segPertShiftF2)
                t_segPertShiftF2=interp1(1:length(t_segPertShiftF2),t_segPertShiftF2,linspace(1,length(t_segPertShiftF2),nInterp));
                segPertShiftF2_tNorm2.(fld0){end+1}=t_segPertShiftF2;
            else
                segPertShiftF2_tNorm2.(fld0){end+1}=nan(1,nInterp);
            end
            
            t_segPertShiftF2=ssf2(taxis1p>=this_utter.youOnset & taxis1p<=this_utter.uayOnset)-of2(taxis1p>=this_utter.youOnset & taxis1p<=this_utter.uayOnset);
            if~isempty(t_segPertShiftF2)
                t_segPertShiftF2=interp1(1:length(t_segPertShiftF2),t_segPertShiftF2,linspace(1,length(t_segPertShiftF2),nInterp));
                segPertShiftF2_tNorm3.(fld0){end+1}=t_segPertShiftF2;
            else
                segPertShiftF2_tNorm3.(fld0){end+1}=nan(1,nInterp);
            end
            
            t_segPertShiftF2=ssf2(taxis1p>=this_utter.uayOnset & taxis1p<=this_utter.yo1Onset)-of2(taxis1p>=this_utter.uayOnset & taxis1p<=this_utter.yo1Onset);
            if~isempty(t_segPertShiftF2)
                t_segPertShiftF2=interp1(1:length(t_segPertShiftF2),t_segPertShiftF2,linspace(1,length(t_segPertShiftF2),nInterp));
                segPertShiftF2_tNorm4.(fld0){end+1}=t_segPertShiftF2;
            else
                segPertShiftF2_tNorm4.(fld0){end+1}=nan(1,nInterp);
            end
            
            t_segPertShiftF2=ssf2(taxis1p>=this_utter.yo1Onset & taxis1p<=this_utter.yo1End)-of2(taxis1p>=this_utter.yo1Onset & taxis1p<=this_utter.yo1End);
            if~isempty(t_segPertShiftF2)
                t_segPertShiftF2=interp1(1:length(t_segPertShiftF2),t_segPertShiftF2,linspace(1,length(t_segPertShiftF2),nInterp));
                segPertShiftF2_tNorm5.(fld0){end+1}=t_segPertShiftF2;
            else
                segPertShiftF2_tNorm5.(fld0){end+1}=nan(1,nInterp);
            end
            
            t_segPertShiftF2=ssf2(taxis1p>=this_utter.yo1End & taxis1p<=this_utter.yo2Onset)-of2(taxis1p>=this_utter.yo1End & taxis1p<=this_utter.yo2Onset);
            if~isempty(t_segPertShiftF2)
                t_segPertShiftF2=interp1(1:length(t_segPertShiftF2),t_segPertShiftF2,linspace(1,length(t_segPertShiftF2),nInterp));
                segPertShiftF2_tNorm6.(fld0){end+1}=t_segPertShiftF2;
            else
                segPertShiftF2_tNorm6.(fld0){end+1}=nan(1,nInterp);
            end
        end

        % Time-normalized F2 velocity: anchorage 1
        t_velF2=diff(this_utter.traj_F2) / frameDur;    % Hz/sec
        t_segVelF2=t_velF2(taxis1<=this_utter.youOnset-this_utter.iouOnset);
        t_segVelF2=interp1(1:length(t_segVelF2),t_segVelF2,linspace(1,length(t_segVelF2),nInterp));
        segVelF2_tNorm1.(fld0){end+1}=t_segVelF2;

        % Perturbation-onset aligned (POA) formant trajectories 
        idx1=find(this_utter.traj_F2<pdata.pertField.F2UB,1);
        if ~isempty(idx1) 
            if idx1*frameDur>=this_utter.uTime-this_utter.iouOnset;
                fprintf('Warning: %s: POA trajectory extraction failed: idx1 larger than expected.\n',this_utter.rawDataFN);
            else                        
                trajF1_POA.(fld0){end+1}=this_utter.traj_F1(idx1:end);                        
                trajF2_POA.(fld0){end+1}=this_utter.traj_F2(idx1:end);
                traj_POA_fn.(fld0){end+1}=this_utter.rawDataFN;
            end
        else
            fprintf('Warning: %s: POA trajectory extraction failed.\n',this_utter.rawDataFN);
        end
        
        % Calculate the amount of timing perturbation, in the case of a
        % accel/decel perturbation
        if ~isequal(et,'apstv') && ~isequal(et,'apstv_STUT')
            if isequal(this_utter.pertStr,'accel') || isequal(this_utter.pertStr,'decel')
                taxis1=0:frameDur:frameDur*(length(f2Trajs.(fld0){end})-1);
                seg_f2=f2Trajs.(fld0){end}(taxis1<=this_utter.youOnset-this_utter.iouOnset);
                seg_sf2=f2Trajs_pert.(fld0){end}(taxis1<=this_utter.youOnset-this_utter.iouOnset);
                seg_f2(seg_f2==0)=NaN;
                seg_sf2(seg_sf2==0)=NaN;
                [foo,idx_minF2]=nanmin(seg_f2);
                [foo,idx_s_minF2]=nanmin(seg_sf2);
                this_utter.tPert_minF2=(idx_s_minF2-idx_minF2)*frameDur;
                timeShifts.minF2.(fld0)(end+1)=(idx_s_minF2-idx_minF2)*frameDur;            

                seg_f2=f2Trajs.(fld0){end}(taxis1<=this_utter.uTime-this_utter.iouOnset);
                seg_sf2=f2Trajs_pert.(fld0){end}(taxis1<=this_utter.uTime-this_utter.iouOnset);
                iuMidF=mean([this_utter.iF2,this_utter.uF2]);
                seg_f2(seg_f2==0)=NaN;
                seg_sf2(seg_sf2==0)=NaN;
                idx_iuMidF=find(seg_f2(1:end-1)>=iuMidF & seg_f2(2:end)<iuMidF);                        
                idx_s_iuMidF=find(seg_sf2(1:end-1)>=iuMidF & seg_sf2(2:end)<iuMidF);

                if isempty(idx_iuMidF) || isempty(idx_s_iuMidF) || length(idx_s_iuMidF) > 1
                    this_utter.tPert_iuMidF=NaN;
                else
                    if length(idx_iuMidF)>1 
                        idx_iuMidF=mean(idx_iuMidF);
                        this_utter.tPert_iuMidF=(idx_s_iuMidF-idx_iuMidF)*frameDur;
                    elseif length(idx_s_iuMidF)>1
                        idx_s_iuMidF=mean(idx_s_iuMidF);
                        this_utter.tPert_iuMidF=(idx_s_iuMidF-idx_iuMidF)*frameDur;
                    else
                        idx_iuMidF=idx_iuMidF+frameDur*(seg_f2(idx_iuMidF)-iuMidF)/(seg_f2(idx_iuMidF)-seg_f2(idx_iuMidF+1));                    
                        idx_s_iuMidF=idx_s_iuMidF+frameDur*(seg_sf2(idx_s_iuMidF)-iuMidF)/(seg_sf2(idx_s_iuMidF)-seg_sf2(idx_s_iuMidF+1));
                        this_utter.tPert_iuMidF=(idx_s_iuMidF-idx_iuMidF)*frameDur;
                    end
                end
                timeShifts.iuMidF.(fld0)(end+1)=this_utter.tPert_iuMidF;

                pdata.utters.(fld){i2}=this_utter;
            end
        end
        
        % Obtain the rms trajectory (from the rawDataFile)
        taxis_1=0 : (frameDur) : (frameDur*(length(data.fmts(:, 2))-1));
        [foo, idx_0] = min(abs(taxis_1 - this_utter.iouOnset));
        [foo, idx_1] = min(abs(taxis_1 - this_utter.yo2Onset));
        t_rmsTraj = data.rms(idx_0 : idx_1, 1);
        t_rmsTraj = log10(t_rmsTraj);
        t_rmsTraj = (t_rmsTraj - mean(t_rmsTraj)) / std(t_rmsTraj);
        rmsTrajs.(fld0){end + 1} = t_rmsTraj;
        
        clear('data');
    end
    
end

%% Calculate average spectrograms
avgSpect_iu2=protoStruct_array;
avgSpect_iu2_df=f(2)-f(1);
avgSpect_iu2_dt=t(2)-t(1);
pertTLens=[];
for i1=1:numel(pertTypes)
    fld=pertTypes{i1};
    tlens=nan(1,numel(spect_iu2.(fld)));
    flens=nan(1,numel(spect_iu2.(fld)));
    for i2=1:numel(spect_iu2.(fld))
        flens(i2)=size(spect_iu2.(fld){i2},1);
        tlens(i2)=size(spect_iu2.(fld){i2},2);
    end
    min_tLen=min(tlens);
    min_fLen=min(flens);
    
    mat_spect=nan(min_fLen,min_tLen,numel(spect_iu2.(fld)));
    for i2=1:numel(spect_iu2.(fld))
        mat_spect(:,:,i2)=spect_iu2.(fld){i2}(1:min_fLen,1:min_tLen);
    end
    pertTLens(end+1)=min_tLen;
    
    avgSpect_iu2.(fld)=mean(mat_spect,3);
end

for i1=1:numel(pertTypes)
    fld=pertTypes{i1};
    figure('Position',[100,100,800,600]);
    t_axis=0:avgSpect_iu2_dt:avgSpect_iu2_dt*(size(avgSpect_iu2.(fld),2)-1);
    f_axis=0:avgSpect_iu2_df:avgSpect_iu2_df*(size(avgSpect_iu2.(fld),1)-1);
    imagesc(t_axis,f_axis,avgSpect_iu2.(fld));
    axis xy;
    hold on;
    title(fld);
    xlabel('Time (s)'); ylabel('Frequency (Hz)');
    set(gca,'XLim',[0,avgSpect_iu2_dt*(max(pertTLens)-1)])
end

%% Calculate normT1 measures (between i and y)
a=avgTrace1(segTrajF2_tNorm1.noPert);
a=avgTrace1(segTrajF2_tNorm1.noPert);
[foo,idx_min]=min(a(:,1));
normT1_uTime=idx_min;
flds=fields(segTrajF2_tNorm1);
for i0=1:numel(config.normT1_bw_iu_ratio)
    for i1=1:length(flds)
        fld=flds{i1};
        for i2=1:numel(segTrajF2_tNorm1.(fld))
            normT1_bw_iu_f2{i0}.(fld)(end+1)=...
                segTrajF2_tNorm1.(fld){i2}(round(1+config.normT1_bw_iu_ratio(i0)*(normT1_uTime-1)));        
%             normT1_bw_uy_f2.(fld)=[normT1_bw_uy_f2.(fld),...
%                 segTrajF2_tNorm1.(fld){i2}(round(normT1_uTime+normT1_bw_uy_ratio*(length(segTrajF2_tNorm1.(fld){i2})-normT1_uTime)))];        
        end
    end
end
for i0=1:length(config.normT1_bw_uy_ratio)
    for i1=1:numel(flds)
        fld=flds{i1};
        for i2=1:numel(segTrajF2_tNorm1.(fld)) 
            normT1_bw_uy_f2{i0}.(fld)(end+1)=...
                segTrajF2_tNorm1.(fld){i2}(round(normT1_uTime+config.normT1_bw_uy_ratio(i0)*(length(segTrajF2_tNorm1.(fld){i2})-normT1_uTime)));
        end
    end
end

%% Calculate normT2 measures (between u and y)
% flds=fields(segTrajF2_tNorm2);
flds=fields(segTraj.F2{2});
for i0=1:length(config.normT2_bw_uy_ratio)
    for i1=1:length(flds)
        fld=flds{i1};
        for i2=1:length(segTraj.F2{2}.(fld)) 
            normT2_bw_uy_f2{i0}.(fld)(end+1)=...
                segTraj.F2{2}.(fld){i2}(round(1+config.normT2_bw_uy_ratio(i0)*(length(segTraj.F2{2}.(fld){i2})-1)));
        end
    end
end

%% Calculate normT3 measures (between y and u in "you")
% flds=fields(segTrajF2_tNorm3);
flds=fields(segTraj.F2{3});
for i0=1:length(config.normT3_bw_yu_ratio)
    for i1=1:length(flds)
        fld=flds{i1};
        for i2=1:length(segTraj.F2{3}.(fld)) 
            normT3_bw_yu_f2{i0}.(fld)(end+1)=...
                segTraj.F2{3}.(fld){i2}(round(1+config.normT3_bw_yu_ratio(i0)*(length(segTraj.F2{3}.(fld){i2})-1)));
        end
    end
end

%% Construct the normalized trajectory between [i] and [j]3
% On the normalized time axis
%   0 - [i]
%   1 - [u]_1
%   2 - [j]_1
%   3 - [u]_2
%   4 - [j]_2
%   5 - [u]_3
%   6 - [j]_3
fullTNorm_tAxis=0 : 1 / (nInterp-1) : 6;

fullTNormF1=protoStruct_cell;
fullTNormF2=protoStruct_cell;

avg_fullTNormF1=protoStruct_cell;
avg_fullTNormF1_nn=protoStruct_cell;
avg_fullTNormF1_cs=protoStruct_cell;

fullTNormF1_nn = protoStruct_cell;  % Nearest neighbor
fullTNormF1_cs = protoStruct_cell;  % Cubic spline

avg_fullTNormF2=protoStruct_cell;
avg_fullTNormF2_nn=protoStruct_cell;
avg_fullTNormF2_cs=protoStruct_cell;

fullTNormF2_nn = protoStruct_cell;  % Nearest neighbor
fullTNormF2_cs = protoStruct_cell;  % Cubic spline

fullTNormF2_pert=protoStruct_cell;
avg_fullTNormF2_pert=protoStruct_cell;

fullTNormF2_FN=protoStruct_cell;
avg_fullTNormF2_FN=protoStruct_cell;

fullTNormF2_AT=protoStruct_array_AT;     % AT stands for "adaptation test"
avg_fullTNormF2_AT=protoStruct_array_AT;

for i1=1:numel(flds)
    fld=flds{i1};
    fullTNormF2.(fld)=nan(0,length(fullTNorm_tAxis));
    fullTNormF2_nn.(fld)=nan(0,length(fullTNorm_tAxis));
    fullTNormF2_cs.(fld)=nan(0,length(fullTNorm_tAxis));
    fullTNormF2_pert.(fld)=nan(0,length(fullTNorm_tAxis));
    for i2=1:numel(segTrajF2_tNorm1a.(fld))
        fullTNormF2.(fld)=[fullTNormF2.(fld); [segTrajF2_tNorm1a.(fld){i2}(1:end), segTrajF2_tNorm2.(fld){i2}(2:end), ...
                               segTrajF2_tNorm3.(fld){i2}(2:end), segTrajF2_tNorm4.(fld){i2}(2:end), ...
                               segTrajF2_tNorm5.(fld){i2}(2:end), segTrajF2_tNorm6.(fld){i2}(2:end)]];
        fullTNormF2_nn.(fld)=[fullTNormF2_nn.(fld); [segTrajF2_tNorm1a_nn.(fld){i2}(1:end), segTrajF2_tNorm2_nn.(fld){i2}(2:end), ...
                               segTrajF2_tNorm3_nn.(fld){i2}(2:end), segTrajF2_tNorm4_nn.(fld){i2}(2:end), ...
                               segTrajF2_tNorm5_nn.(fld){i2}(2:end), segTrajF2_tNorm6_nn.(fld){i2}(2:end)]];
        fullTNormF2_cs.(fld)=[fullTNormF2_cs.(fld); [segTrajF2_tNorm1a_cs.(fld){i2}(1:end), segTrajF2_tNorm2_cs.(fld){i2}(2:end), ...
                               segTrajF2_tNorm3_cs.(fld){i2}(2:end), segTrajF2_tNorm4_cs.(fld){i2}(2:end), ...
                               segTrajF2_tNorm5_cs.(fld){i2}(2:end), segTrajF2_tNorm6_cs.(fld){i2}(2:end)]];
        
        if isequal(fld,'noPert') || isequal(fld,'none') % AT: adaptation test
            fldAT=['aft_',prevPertTypes.noPert{i2}];
            fullTNormF2_AT.(fldAT)=[fullTNormF2_AT.(fldAT);fullTNormF2.(fld)(end,:)];
        else
            fullTNormF2_pert.(fld)=[fullTNormF2_pert.(fld); [segTrajF2_tNorm1a_pert.(fld){i2}(1:end), segTrajF2_tNorm2_pert.(fld){i2}(2:end), ...
                   segTrajF2_tNorm3_pert.(fld){i2}(2:end), segTrajF2_tNorm4_pert.(fld){i2}(2:end), ...
                   segTrajF2_tNorm5_pert.(fld){i2}(2:end), segTrajF2_tNorm6_pert.(fld){i2}(2:end)]];
        end
    end
    
    t_sum=sum(fullTNormF2.(fld)');
    idx_ok=find(~isnan(t_sum));
    fullTNormF2.(fld)=fullTNormF2.(fld)(idx_ok,:);
    
    t_avg=nanmean(fullTNormF2.(fld));
    t_std=nanstd(fullTNormF2.(fld));
    t_n=size(fullTNormF2.(fld),1)*ones(1,size(fullTNormF2.(fld),2));    
    avg_fullTNormF2.(fld)=[t_avg',t_std',t_n'];
    
    t_avg=nanmean(fullTNormF2_nn.(fld));
    t_std=nanstd(fullTNormF2_nn.(fld));
    t_n=size(fullTNormF2_nn.(fld),1)*ones(1,size(fullTNormF2_nn.(fld),2));    
    avg_fullTNormF2_nn.(fld)=[t_avg',t_std',t_n'];
    
    t_avg=nanmean(fullTNormF2_cs.(fld));
    t_std=nanstd(fullTNormF2_cs.(fld));
    t_n=size(fullTNormF2_cs.(fld),1)*ones(1,size(fullTNormF2_cs.(fld),2));    
    avg_fullTNormF2_cs.(fld)=[t_avg',t_std',t_n'];
    
    if ~(isequal(fld,'noPert') || isequal(fld,'none'))
        t_avg=nanmean(fullTNormF2_pert.(fld));
        t_std=nanstd(fullTNormF2_pert.(fld));
        t_n=size(fullTNormF2_pert.(fld),1)*ones(1,size(fullTNormF2.(fld),2));
    
        avg_fullTNormF2_pert.(fld)=[t_avg',t_std',t_n'];
    end
end

fldsAT=fields(fullTNormF2_AT);
for i1=1:numel(fldsAT)
    fld=fldsAT{i1};
    
    t_avg=nanmean(fullTNormF2_AT.(fld));
    t_std=nanstd(fullTNormF2_AT.(fld));
    t_n=size(fullTNormF2_AT.(fld),1)*ones(1,size(fullTNormF2_AT.(fld),2));
    
    avg_fullTNormF2_AT.(fld)=[t_avg',t_std',t_n'];
end

fullTNormPertShiftF2=protoStruct_cell;
avg_fullTNormPertShiftF2=protoStruct_cell;
for i1=1:numel(flds)
    fld=flds{i1};
    fullTNormPertShift.(fld)=nan(0,length(fullTNorm_tAxis));
    for i2=1:numel(segPertShiftF2_tNorm1.(fld))
        fullTNormPertShift.(fld)=[fullTNormPertShift.(fld); [segPertShiftF2_tNorm1.(fld){i2}(1:end), segPertShiftF2_tNorm2.(fld){i2}(2:end), ...
                               segPertShiftF2_tNorm3.(fld){i2}(2:end), segPertShiftF2_tNorm4.(fld){i2}(2:end), ...
                               segPertShiftF2_tNorm5.(fld){i2}(2:end), segPertShiftF2_tNorm6.(fld){i2}(2:end)]];
    end
    
    t_sum=sum(fullTNormPertShift.(fld)');
    idx_ok=find(~isnan(t_sum));
    fullTNormPertShift.(fld)=fullTNormPertShift.(fld)(idx_ok,:);
    
    t_avg=nanmean(fullTNormPertShift.(fld));
    t_std=nanstd(fullTNormPertShift.(fld));
    t_n=size(fullTNormPertShift.(fld),1)*ones(1,size(fullTNormPertShift.(fld),2));
    
    avg_fullTNormPertShiftF2.(fld)=[t_avg',t_std',t_n'];
end

%%
frameDur=16/12e3;
    
for i1=1:length(fields(f2Trajs))
    flds=fields(f2Trajs);
    fld=flds{i1};

    avg_f1Traj.(fld)=avgTrace1(f1Trajs.(fld));
    avg_f2Traj.(fld)=avgTrace1(f2Trajs.(fld));
    avg_rmsTraj.(fld)=avgTrace1(rmsTrajs.(fld));

    if ~isempty(f2Trajs_pert.(fld))
        avg_f2Traj_pert.(fld)=avgTrace1(f2Trajs_pert.(fld));
        avg_f1Traj_pert.(fld)=avgTrace1(f1Trajs_pert.(fld));
        avg_f2_pertShift.(fld)=avgTrace1(f2_pertShift.(fld));
        avg_f1_pertShift.(fld)=avgTrace1(f1_pertShift.(fld));
    end

end


%% Compute IOA_TN (IOA, time-normalized) trajectories
if isfield(avg_f2Traj,'none')
    fld1='none';
else
    fld1='noPert';
end
t_traj_f2=avg_f2Traj.(fld1)(:,1);
diff_t_traj_f2=diff(t_traj_f2);
idx_turn=find(diff_t_traj_f2(1:end-1).*diff_t_traj_f2(2:end)<0);
idx_turn = idx_turn(idx_turn > 5);
f2LB_IOA=t_traj_f2(idx_turn(1)+1);
f2UB_IOA=t_traj_f2(idx_turn(2)+1);

d_t_traj_f2=diff(t_traj_f2);
idx_peak_1=find(d_t_traj_f2(1:end-1)>0 & d_t_traj_f2(2:end)<0);
idx_peak_1 = idx_peak_1(idx_peak_1 > 5);

idxNormInt_IOA=idx_peak_1(1)+1;

t_traj_f1=avg_f1Traj.(fld1)(:,1);
diff_t_traj_f1=diff(t_traj_f1);
idx_maxima=find(diff_t_traj_f1(1:end-1).*diff_t_traj_f1(2:end)<0 & diff_t_traj_f1(1:end-1)>0);
idx_maxima = idx_maxima(idx_maxima > 5);
idx_maxima=idx_maxima(1);            
idx_minima=find(diff_t_traj_f1(1:end-1).*diff_t_traj_f1(2:end)<0 & diff_t_traj_f1(1:end-1)<0);
idx_minima = idx_minima(idx_minima > 5);
idx_minima=idx_minima(idx_minima>idx_maxima);            
f1LB_IOA=t_traj_f1(idx_minima(1)+1);
f1UB_IOA=t_traj_f1(idx_maxima+1);

trajF2_IOA_TN=protoStruct_cell;
trajF1_IOA_TN=protoStruct_cell;
trajF2_IOA_TN_FN=protoStruct_cell;
trajF1_IOA_TN_FN=protoStruct_cell;
trajF2_IOA_FN=protoStruct_cell;
trajF1_IOA_FN=protoStruct_cell;

% trajF2_pert_IOA=protoStruct_cell;
% trajF1_pert_IOA=protoStruct_cell;
trajF2_pert_IOA_TN=protoStruct_cell;
trajF1_pert_IOA_TN=protoStruct_cell;
trajF2_pert_IOA_FN=protoStruct_cell;
trajF1_pert_IOA_FN=protoStruct_cell;
trajF2_pert_IOA_TN_FN=protoStruct_cell;
trajF1_pert_IOA_TN_FN=protoStruct_cell;
f2_pertShift_IOA_TN=protoStruct_cell;
f1_pertShift_IOA_TN=protoStruct_cell;
f2_pertShift_IOA_FN=protoStruct_cell;
f1_pertShift_IOA_FN=protoStruct_cell;
f2_pertShift_IOA_TN_FN=protoStruct_cell;
f1_pertShift_IOA_TN_FN=protoStruct_cell;

flds=fields(f2Trajs);
for i1=1:length(flds)
    fld=flds{i1};
    for i2=1:length(f2Trajs.(fld))
        t_ind_traj_f2=f2Trajs.(fld){i2};
        t_ind_traj_f1=f1Trajs.(fld){i2};

        fNormF2=(t_ind_traj_f2-f2LB_IOA)/(f2UB_IOA-f2LB_IOA);
        fNormF1=(t_ind_traj_f1-f1LB_IOA)/(f1UB_IOA-f1LB_IOA);

        oldTAxis=1:length(t_ind_traj_f2);
        newTAxis=linspace(1,idxNormInt_IOA,400);
        dt=newTAxis(2)-newTAxis(1);
        newTAxis=[newTAxis,newTAxis(end)+dt:dt:oldTAxis(end)];
        trajF2_IOA_TN.(fld){end+1}=interp1(oldTAxis,t_ind_traj_f2,newTAxis);
        trajF1_IOA_TN.(fld){end+1}=interp1(oldTAxis,t_ind_traj_f1,newTAxis);
        
        trajF2_IOA_FN.(fld){end+1}=fNormF2;
        trajF1_IOA_FN.(fld){end+1}=fNormF1;
        
        trajF2_IOA_TN_FN.(fld){end+1}=interp1(oldTAxis,fNormF2,newTAxis);
        trajF1_IOA_TN_FN.(fld){end+1}=interp1(oldTAxis,fNormF1,newTAxis);
        
    end
    
    if ~(isequal(fld,'none') || isequal(fld,'noPert'))
        for i2=1:length(f2Trajs_pert.(fld))
            t_ind_traj_f2p=f2Trajs_pert.(fld){i2};
            t_ind_traj_f1p=f1Trajs_pert.(fld){i2};
            t_f2_pertShift=f2_pertShift.(fld){i2};
            t_f1_pertShift=f1_pertShift.(fld){i2};

            fNormF2=(t_ind_traj_f2p-f2LB_IOA)/(f2UB_IOA-f2LB_IOA);
            fNormF1=(t_ind_traj_f1p-f1LB_IOA)/(f1UB_IOA-f1LB_IOA);
            pertShiftNormF2=t_f2_pertShift/(f2UB_IOA-f2LB_IOA);
            pertShiftNormF1=t_f1_pertShift/(f1UB_IOA-f1LB_IOA);

            oldTAxis=1:length(t_ind_traj_f2p);
            newTAxis=linspace(1,idxNormInt_IOA,400);
            dt=newTAxis(2)-newTAxis(1);
            newTAxis=[newTAxis,newTAxis(end)+dt:dt:oldTAxis(end)];
            
%             trajF2_pert_IOA.(fld){end+1}=t_ind_traj_f2p;
%             trajF1_pert_IOA.(fld){end+1}=t_ind_traj_f1p;
            
            trajF2_pert_IOA_TN.(fld){end+1}=interp1(oldTAxis,t_ind_traj_f2p,newTAxis);
            trajF1_pert_IOA_TN.(fld){end+1}=interp1(oldTAxis,t_ind_traj_f1p,newTAxis);

            trajF2_pert_IOA_FN.(fld){end+1}=fNormF2;
            trajF1_pert_IOA_FN.(fld){end+1}=fNormF1;
            
            trajF2_pert_IOA_TN_FN.(fld){end+1}=interp1(oldTAxis,fNormF2,newTAxis);
            trajF1_pert_IOA_TN_FN.(fld){end+1}=interp1(oldTAxis,fNormF1,newTAxis);
            
            f2_pertShift_IOA_TN.(fld){end+1}=interp1(oldTAxis,t_f2_pertShift,newTAxis);
            f1_pertShift_IOA_TN.(fld){end+1}=interp1(oldTAxis,t_f1_pertShift,newTAxis);
            
            f2_pertShift_IOA_FN.(fld){end+1}=pertShiftNormF2;
            f1_pertShift_IOA_FN.(fld){end+1}=pertShiftNormF1;
            
            f2_pertShift_IOA_TN_FN.(fld){end+1}=interp1(oldTAxis,pertShiftNormF2,newTAxis);
            f1_pertShift_IOA_TN_FN.(fld){end+1}=interp1(oldTAxis,pertShiftNormF1,newTAxis);
        end
    end
    
    avg_f2Traj_IOA_TN.(fld)=avgTrace1(trajF2_IOA_TN.(fld));
    avg_f1Traj_IOA_TN.(fld)=avgTrace1(trajF1_IOA_TN.(fld));
    avg_f2Traj_IOA_FN.(fld)=avgTrace1(trajF2_IOA_FN.(fld));
    avg_f1Traj_IOA_FN.(fld)=avgTrace1(trajF1_IOA_FN.(fld));
    avg_f2Traj_IOA_TN_FN.(fld)=avgTrace1(trajF2_IOA_TN_FN.(fld));
    avg_f1Traj_IOA_TN_FN.(fld)=avgTrace1(trajF1_IOA_TN_FN.(fld));

    avg_f2Traj_pert_IOA_TN.(fld)=avgTrace1(trajF2_pert_IOA_TN.(fld));
    avg_f1Traj_pert_IOA_TN.(fld)=avgTrace1(trajF1_pert_IOA_TN.(fld));
    avg_f2Traj_pert_IOA_FN.(fld)=avgTrace1(trajF2_pert_IOA_FN.(fld));
    avg_f1Traj_pert_IOA_FN.(fld)=avgTrace1(trajF1_pert_IOA_FN.(fld));
    avg_f2Traj_pert_IOA_TN_FN.(fld)=avgTrace1(trajF2_pert_IOA_TN_FN.(fld));
    avg_f1Traj_pert_IOA_TN_FN.(fld)=avgTrace1(trajF1_pert_IOA_TN_FN.(fld));
    
    avg_f2_pertShift_IOA_TN.(fld)=avgTrace1(f2_pertShift_IOA_TN.(fld));
    avg_f2_pertShift_IOA_FN.(fld)=avgTrace1(f2_pertShift_IOA_FN.(fld));
    avg_f2_pertShift_IOA_TN_FN.(fld)=avgTrace1(f2_pertShift_IOA_TN_FN.(fld));
    
    avg_f1_pertShift_IOA_TN.(fld)=avgTrace1(f1_pertShift_IOA_TN.(fld));
    avg_f1_pertShift_IOA_FN.(fld)=avgTrace1(f1_pertShift_IOA_FN.(fld));
    avg_f1_pertShift_IOA_TN_FN.(fld)=avgTrace1(f1_pertShift_IOA_TN_FN.(fld));
end

%% Compute the FTN time-normalized & frequency-normalized (FN) trajectories

f2LB_FTN=avg_fullTNormF2.noPert(nInterp);
f2UB_FTN=avg_fullTNormF2.noPert(nInterp*2-1);

flds=fields(fullTNormF2);
for i1=1:numel(flds)
    fld=flds{i1};
    fullTNormF2_FN.(fld)=(fullTNormF2.(fld)-f2LB_FTN)/(f2UB_FTN-f2LB_FTN);
    
    t_avg=nanmean(fullTNormF2_FN.(fld));
    t_std=nanstd(fullTNormF2_FN.(fld));
    t_n=size(fullTNormF2_FN.(fld),1)*ones(1,size(fullTNormF2_FN.(fld),2));
    
    avg_fullTNormF2_FN.(fld)=[t_avg',t_std',t_n'];
    
    if ~(isequal(fld,'noPert') || isequal(fld,'none'))
        fullTNormF2_pert_FN.(fld)=(fullTNormF2_pert.(fld)-f2LB_FTN)/(f2UB_FTN-f2LB_FTN);
    
        t_avg=nanmean(fullTNormF2_pert_FN.(fld));
        t_std=nanstd(fullTNormF2_pert_FN.(fld));
        t_n=size(fullTNormF2_FN.(fld),1)*ones(1,size(fullTNormF2_pert_FN.(fld),2));
    
        avg_fullTNormF2_pert_FN.(fld)=[t_avg',t_std',t_n'];
    end
end

%%

pertTypes_pertOnly=pertTypes(setxor(1:numel(pertTypes),idx_noPert));
for j1=1:length(pertTypes)
    fld=pertTypes{j1};
    f2Trajs_pert_mat.(fld)=[];
    
    nTrials=length(f2Trajs_pert.(fld));
    minLen=Inf;
    for j2=1:length(f2Trajs_pert.(fld))
        if length(f2Trajs_pert.(fld){j2})<minLen
            minLen=length(f2Trajs_pert.(fld){j2});
        end
    end
    for j2=1:length(f2Trajs_pert.(fld))
        f2Trajs_pert_mat.(fld)=[f2Trajs_pert_mat.(fld),f2Trajs_pert.(fld){j2}(1:minLen)];
    end
    
end

% 
% for i1=1:length(fields(f2Trajs))
%     flds=fields(f2Trajs);
%     fld=flds{i1};
% 
%     if ~isempty(f2Trajs_pert.(fld))
%         avg_f2Traj_pert.(fld)=avgTrace1(f2Trajs_pert.(fld));
%         avg_f1Traj_pert.(fld)=avgTrace1(f1Trajs_pert.(fld));
%         avg_f2_pertShift.(fld)=avgTrace1(f2_pertShift.(fld));
%         avg_f1_pertShift.(fld)=avgTrace1(f1_pertShift.(fld));
%     end
% 
% end

%% DEBUG: The relation between yF2 (manual) and fullTNormF2
yF2_FTN=struct;
flds1=fields(fullTNormF2);
figure;
axs=[];
ays=[];
afns={};


if ~isempty(fsic(varargin,'pick'))
    for i1=1:numel(flds1)
        fld=flds1{i1};

        yF2_FTN.(fld)=nan(1,size(fullTNormF2.(fld),1));

        for i2=1:size(fullTNormF2.(fld),1)
            [foo,idx_2]=min(abs(fullTNorm_tAxis-2));
            yF2_FTN.(fld)(i2)=fullTNormF2.(fld)(i2,idx_2);
        end

        plot(yF2.(fld),yF2_FTN.(fld),'o','color',colors.(fld));
        axs=[axs,yF2.(fld)];
        ays=[ays,yF2_FTN.(fld)];
        afns=[afns,a_rawFNs.(fld)];
        hold on;
    end
    
    bContinue=0;
    while bContinue==0
        coord=ginput(1);

        xs=get(gca,'XLim'); ys=get(gca,'YLim');
        if coord(1)<xs(1) || coord(1)>xs(2) || coord(2)<ys(1) || coord(2)>ys(2)
            bContinue=1;
        else
            dist=(axs-coord(1)).^2+(ays-coord(2)).^2;
            [foo,idx_min]=min(dist);
            fprintf('rawFN: %s\n',afns{idx_min});
        end
    end
end

%% Average and plot F1 and F2 trajectories (original and perturbed)
figure('Position',[50,50,1000,700]);

frameDur=16/12e3;
for i0=1:2
    subplot('Position',[0.1+(i0-1)*0.45,0.65+0.015,0.4,0.3]);
    set(gca,'FontSize',fontSize);
    
    for i1=1:length(fields(f2Trajs))
        flds=fields(f2Trajs);
        fld=flds{i1};

        for i2=1:length(f2Trajs.(fld))
            t_f1=f1Trajs.(fld){i2};
            t_f2=f2Trajs.(fld){i2};
            taxis1=0:(frameDur):(frameDur*(length(t_f2)-1));
            if ~isempty(find(t_f2==0))
%                 pause(0);
            end

            if i0==1
                if ~isempty(fsic(flds,'up'))
                    if isequal(fld,'up') continue; end
                else 
                    if isequal(fld,'decel') continue; end
                end
            else
                if ~isempty(fsic(flds,'up'))
                    if isequal(fld,'down') continue; end
                else
                    if isequal(fld,'accel') continue; end
                end
            end
            plot(taxis1*1e3,t_f1,'Color',colors.(fld),'LineWidth',0.5);
            plot(taxis1*1e3,t_f2,'Color',colors.(fld),'LineWidth',0.5);
            hold on;            
        end
    end
    
    if isfield(f2Trajs_pert,'down')
        if i0==1
            pertFld='down';
        else
            pertFld='up';
        end
    else
        if i0==1
            pertFld='accel';
        else
            pertFld='decel';
        end
    end
    for i2=1:length(f2Trajs_pert.(pertFld))
        t_f2p=f2Trajs_pert.(pertFld){i2};
        taxis1=0:(frameDur):(frameDur*(length(t_f2p)-1));
        plot(taxis1*1e3,t_f2p,'--','Color',colors_pert.(pertFld));
        hold on;
    end
    
    if ~isempty(fsic(varargin,'trajXLim'))
        trajXLim=varargin{fsic(varargin,'trajXLim')+1};
        set(gca,'XLim',trajXLim);
    end
    if ~isempty(fsic(varargin,'trajYLim'))
        trajYLim=varargin{fsic(varargin,'trajYLim')+1};
        set(gca,'YLim',trajYLim);
    end
    
    set(gca,'XTickLabel',[]);
    
    if i0==1
        ylabel('F2 (Hz)','FontSize',fontSize);
%         xlabel('Time (ms)','FontSize',fontSize);
    end
    
    xs=get(gca,'XLim'); ys=get(gca,'YLim');
    if i0==1
        if ~isempty(fsic(flds,'up'))
            title('Down - noPert');
            ezlegend([xs(1)+0.025*range(xs),ys(1)+0.75*range(ys),0.375*range(xs),0.225*range(ys)],0.4,...
                {'noPert','Down (prod.)','Down (AF)'},{colors.noPert,colors.down,colors_pert.down},...
                [fontSize,fontSize,fontSize]-1,{'-','-','--'},{colors.noPert,colors.down,colors_pert.down},[0.5,0.5,0.5],[0,0,0]);
        else
            title('Accel - noPert');
            ezlegend([xs(1)+0.025*range(xs),ys(1)+0.75*range(ys),0.375*range(xs),0.225*range(ys)],0.4,...
                {'noPert','Accel (prod.)','Accel (AF)'},{colors.noPert,colors.accel,colors_pert.accel},...
                [fontSize,fontSize,fontSize]-1,{'-','-','--'},{colors.noPert,colors.accel,colors_pert.accel},[0.5,0.5,0.5],[0,0,0]);
        end
    else
        if ~isempty(fsic(flds,'up'))
            title('Up - noPert');
            ezlegend([xs(1)+0.025*range(xs),ys(1)+0.75*range(ys),0.375*range(xs),0.225*range(ys)],0.4,...
                {'noPert','Up (prod.)','Up (AF)'},{colors.noPert,colors.up,colors_pert.up},...
                [fontSize,fontSize,fontSize]-1,{'-','-','--'},{colors.noPert,colors.up,colors_pert.up},[0.5,0.5,0.5],[0,0,0]);
        else
            title('Decel - noPert');
            ezlegend([xs(1)+0.025*range(xs),ys(1)+0.75*range(ys),0.375*range(xs),0.225*range(ys)],0.4,...
                {'noPert','Decel (prod.)','Decel (AF)'},{colors.noPert,colors.decel,colors_pert.decel},...
                [fontSize,fontSize,fontSize]-1,{'-','-','--'},{colors.noPert,colors.decel,colors_pert.decel},[0.5,0.5,0.5],[0,0,0]);
        end
    end
    set(gca,'XTickLabel',[]);
end

%%
% figure; 
for i0=1:2
    subplot('Position',[0.1+(i0-1)*0.45,0.35+0.015,0.4,0.3]);
    set(gca,'FontSize',fontSize);
    for i1=1:length(pertTypes)
        fld=pertTypes{i1};
        taxis1=0:(frameDur):(frameDur*(length(avg_f2Traj.(fld)(:,1))-1));
        if i0==1
            if ~isempty(fsic(flds,'up'))
                if isequal(fld,'up') continue; end
            else 
                if isequal(fld,'decel') continue; end
            end           
        else
            if ~isempty(fsic(flds,'up'))
                if isequal(fld,'down') continue; end
            else
                if isequal(fld,'accel') continue; end
            end
        end
        
        if ~(isequal(fld,'noPert') || isequal(fld,'none'))
            taxis1=0:(frameDur):(frameDur*(length(avg_f2Traj_pert.(fld)(:,1))-1));
            plot(taxis1*1e3,avg_f2Traj_pert.(fld)(:,1),'-','Color',colors_pert.(fld)); hold on;
            plot(taxis1*1e3,avg_f2Traj_pert.(fld)(:,1)-avg_f2Traj_pert.(fld)(:,2)./sqrt(avg_f2Traj_pert.(fld)(:,3)),'--','Color',colors_pert.(fld));
            plot(taxis1*1e3,avg_f2Traj_pert.(fld)(:,1)+avg_f2Traj_pert.(fld)(:,2)./sqrt(avg_f2Traj_pert.(fld)(:,3)),'--','Color',colors_pert.(fld));
        end
        
        taxis1=0:(frameDur):(frameDur*(length(avg_f1Traj.(fld)(:,1))-1));
        plot(taxis1*1e3,avg_f1Traj.(fld)(:,1),'-','Color',colors.(fld)); hold on;
        plot(taxis1*1e3,avg_f1Traj.(fld)(:,1)-avg_f1Traj.(fld)(:,2)./sqrt(avg_f1Traj.(fld)(:,3)),'--','Color',colors.(fld));
        plot(taxis1*1e3,avg_f1Traj.(fld)(:,1)+avg_f1Traj.(fld)(:,2)./sqrt(avg_f1Traj.(fld)(:,3)),'--','Color',colors.(fld));
        plot(taxis1*1e3,avg_f2Traj.(fld)(:,1),'-','Color',colors.(fld)); hold on;
        plot(taxis1*1e3,avg_f2Traj.(fld)(:,1)-avg_f2Traj.(fld)(:,2)./sqrt(avg_f2Traj.(fld)(:,3)),'--','Color',colors.(fld));
        plot(taxis1*1e3,avg_f2Traj.(fld)(:,1)+avg_f2Traj.(fld)(:,2)./sqrt(avg_f2Traj.(fld)(:,3)),'--','Color',colors.(fld));
        
        
    end
    
    if i0==1
        ylabel('F2 (Hz, mean\pm 1SEM)');
        if ~isempty(fsic(flds,'up'))
            title('Down - noPert');
            ezlegend([xs(1)+0.025*range(xs),ys(1)+0.75*range(ys),0.375*range(xs),0.225*range(ys)],0.4,...
                {'noPert','Down (prod.)','Down (AF)'},{colors.noPert,colors.down,colors_pert.down},...
                [fontSize,fontSize,fontSize]-1,{'-','-','--'},{colors.noPert,colors.down,colors_pert.down},[0.5,0.5,0.5],[0,0,0]);
        else
            title('Accel - noPert');
            ezlegend([xs(1)+0.025*range(xs),ys(1)+0.75*range(ys),0.375*range(xs),0.225*range(ys)],0.4,...
                {'noPert','Accel (prod.)','Accel (AF)'},{colors.noPert,colors.accel,colors_pert.accel},...
                [fontSize,fontSize,fontSize]-1,{'-','-','--'},{colors.noPert,colors.accel,colors_pert.accel},[0.5,0.5,0.5],[0,0,0]);
        end
    else
        if ~isempty(fsic(flds,'up'))
            title('Up - noPert');
            ezlegend([xs(1)+0.025*range(xs),ys(1)+0.75*range(ys),0.375*range(xs),0.225*range(ys)],0.4,...
                {'noPert','Up (prod.)','Up (AF)'},{colors.noPert,colors.up,colors_pert.up},...
                [fontSize,fontSize,fontSize]-1,{'-','-','--'},{colors.noPert,colors.up,colors_pert.up},[0.5,0.5,0.5],[0,0,0]);
        else
            title('Decel - noPert');
            ezlegend([xs(1)+0.025*range(xs),ys(1)+0.75*range(ys),0.375*range(xs),0.225*range(ys)],0.4,...
                {'noPert','Decel (prod.)','Decel (AF)'},{colors.noPert,colors.decel,colors_pert.decel},...
                [fontSize,fontSize,fontSize]-1,{'-','-','--'},{colors.noPert,colors.decel,colors_pert.decel},[0.5,0.5,0.5],[0,0,0]);
        end
    end

    if ~isempty(fsic(varargin,'trajXLim'))
        trajXLim=varargin{fsic(varargin,'trajXLim')+1};
        set(gca,'XLim',trajXLim);
    end
    if ~isempty(fsic(varargin,'trajYLim'))
        trajYLim=varargin{fsic(varargin,'trajYLim')+1};
        set(gca,'YLim',trajYLim);
    end

    
    set(gca,'XTickLabel',[]);
    
end

%% Plot the deviation from noPert (IOA, un-normalized)
if isfield(avg_f2Traj,'down')
    pertFld1='down';
    pertFld2='up';
else
    pertFld1='accel';
    pertFld2='decel';
end
for i0=1:2
    subplot('Position',[0.1+(i0-1)*0.45,0.05+0.015,0.4,0.3]);
    set(gca,'FontSize',fontSize);
    commonLen=min([size(avg_f2Traj.noPert,1),size(avg_f2Traj.(pertFld1),1),size(avg_f2Traj.(pertFld2),1),...
        size(avg_f2_pertShift.(pertFld1),1),size(avg_f2_pertShift.(pertFld2),1)]);
    taxis2=0:frameDur:frameDur*(commonLen-1);
    if i0==1
        plot(taxis2*1e3,avg_f2_pertShift.(pertFld1)(1:commonLen,1),'Color',colors_pert.(pertFld1));       
        hold on;
        plot(taxis2*1e3,avg_f2_pertShift.(pertFld1)(1:commonLen,1)-avg_f2_pertShift.(pertFld1)(1:commonLen,2)./sqrt(avg_f2_pertShift.(pertFld1)(1:commonLen,2)),...
            '--','Color',colors_pert.(pertFld1));
        plot(taxis2*1e3,avg_f2_pertShift.(pertFld1)(1:commonLen,1)+avg_f2_pertShift.(pertFld1)(1:commonLen,2)./sqrt(avg_f2_pertShift.(pertFld1)(1:commonLen,2)),...
            '--','Color',colors_pert.(pertFld1));
        
        plot(taxis2*1e3,avg_f2Traj.(pertFld1)(1:commonLen,1)-avg_f2Traj.noPert(1:commonLen,1),'Color',colors.(pertFld1));
        hold on;
        plot(taxis2*1e3,avg_f2Traj.(pertFld1)(1:commonLen,1)-avg_f2Traj.noPert(1:commonLen,1)-avg_f2Traj.(pertFld1)(1:commonLen,2)./sqrt(avg_f2Traj.(pertFld1)(1:commonLen,2)),...
            '--','Color',colors.(pertFld1));
        plot(taxis2*1e3,avg_f2Traj.(pertFld1)(1:commonLen,1)-avg_f2Traj.noPert(1:commonLen,1)+avg_f2Traj.(pertFld1)(1:commonLen,2)./sqrt(avg_f2Traj.(pertFld1)(1:commonLen,2)),...
            '--','Color',colors.(pertFld1));
    else
        plot(taxis2*1e3,avg_f2_pertShift.(pertFld2)(1:commonLen,1),'Color',colors_pert.(pertFld2));
        hold on;
        plot(taxis2*1e3,avg_f2_pertShift.(pertFld2)(1:commonLen,1)-avg_f2_pertShift.(pertFld2)(1:commonLen,2)./sqrt(avg_f2_pertShift.(pertFld2)(1:commonLen,2)),...
            '--','Color',colors_pert.(pertFld2));
        plot(taxis2*1e3,avg_f2_pertShift.(pertFld2)(1:commonLen,1)+avg_f2_pertShift.(pertFld2)(1:commonLen,2)./sqrt(avg_f2_pertShift.(pertFld2)(1:commonLen,2)),...
            '--','Color',colors_pert.(pertFld2));        
        
        plot(taxis2*1e3,avg_f2Traj.(pertFld2)(1:commonLen,1)-avg_f2Traj.noPert(1:commonLen,1),'Color',colors.(pertFld2));
        hold on;
        plot(taxis2*1e3,avg_f2Traj.(pertFld2)(1:commonLen,1)-avg_f2Traj.noPert(1:commonLen,1)-avg_f2Traj.(pertFld2)(1:commonLen,2)./sqrt(avg_f2Traj.(pertFld2)(1:commonLen,2)),...
            '--','Color',colors.(pertFld2));
        plot(taxis2*1e3,avg_f2Traj.(pertFld2)(1:commonLen,1)-avg_f2Traj.noPert(1:commonLen,1)+avg_f2Traj.(pertFld2)(1:commonLen,2)./sqrt(avg_f2Traj.(pertFld2)(1:commonLen,2)),...
            '--','Color',colors.(pertFld2));
    end
        

    if ~isempty(fsic(varargin,'trajXLim'))
        trajXLim=varargin{fsic(varargin,'trajXLim')+1};
        set(gca,'XLim',trajXLim);
    end
    % if ~isempty(fsic(varargin,'trajYLim'))
    %     trajYLim=varargin{fsic(varargin,'trajYLim')+1};
    %     set(gca,'YLim',trajYLim);
    % end
    xs=get(gca,'XLim');
    plot(xs,[0,0],'-','Color',[0.5,0.5,0.5]);
    
    xs=get(gca,'XLim'); ys=get(gca,'YLim');
    
    if i0==1
        ylabel('F2 perturbation & compensation (Hz)');
        if ~isempty(fsic(flds,'up'))
            title('Down - noPert');
            ezlegend([xs(1)+0.025*range(xs),ys(1)+0.75*range(ys),0.375*range(xs),0.225*range(ys)],0.4,...
                {'noPert','Down (prod.)','Down (AF)'},{colors.noPert,colors.down,colors_pert.down},...
                [fontSize,fontSize,fontSize]-1,{'-','-','--'},{colors.noPert,colors.down,colors_pert.down},[0.5,0.5,0.5],[0,0,0]);
        else
            title('Accel - noPert');
            ezlegend([xs(1)+0.025*range(xs),ys(1)+0.75*range(ys),0.375*range(xs),0.225*range(ys)],0.4,...
                {'noPert','Accel (prod.)','Accel (AF)'},{colors.noPert,colors.accel,colors_pert.accel},...
                [fontSize,fontSize,fontSize]-1,{'-','-','--'},{colors.noPert,colors.accel,colors_pert.accel},[0.5,0.5,0.5],[0,0,0]);
        end
    else
        if ~isempty(fsic(flds,'up'))
            title('Up - noPert');
            ezlegend([xs(1)+0.025*range(xs),ys(1)+0.75*range(ys),0.375*range(xs),0.225*range(ys)],0.4,...
                {'noPert','Up (prod.)','Up (AF)'},{colors.noPert,colors.up,colors_pert.up},...
                [fontSize,fontSize,fontSize]-1,{'-','-','--'},{colors.noPert,colors.up,colors_pert.up},[0.5,0.5,0.5],[0,0,0]);
        else
            title('Decel - noPert');
            ezlegend([xs(1)+0.025*range(xs),ys(1)+0.75*range(ys),0.375*range(xs),0.225*range(ys)],0.4,...
                {'noPert','Decel (prod.)','Decel (AF)'},{colors.noPert,colors.decel,colors_pert.decel},...
                [fontSize,fontSize,fontSize]-1,{'-','-','--'},{colors.noPert,colors.decel,colors_pert.decel},[0.5,0.5,0.5],[0,0,0]);
        end
    end
    
    xlabel('Time re [i] (ms)');
end

%% F2 trajectory plots: POA (perturbation-onset aligned)
figure; set(gca,'FontSize',fontSize);
frameDur=16/12e3;
for i1=1:length(pertTypes)
	fld=pertTypes{i1};
    avg_f1Traj_POA.(fld)=avgTrace1(trajF1_POA.(fld));    
    avg_f2Traj_POA.(fld)=avgTrace1(trajF2_POA.(fld));
	for i2=1:length(trajF2_POA.(fld))
		t_f2=trajF2_POA.(fld){i2};       
		taxis1=0:(frameDur):(frameDur*(length(t_f2)-1));
		plot(taxis1,t_f2,'Color',colors.(fld)); 
		hold on;
	end
end
xlabel('Time (sec)'); ylabel('F2 (Hz)');

figure; set(gca,'FontSize',fontSize);
for i1=1:length(pertTypes)
	fld=pertTypes{i1};        
    taxis1=0:(frameDur):(frameDur*(length(avg_f2Traj_POA.(fld)(:,1))-1));
    plot(taxis1,avg_f2Traj_POA.(fld)(:,1),'-','Color',colors.(fld)); hold on;
    plot(taxis1,avg_f2Traj_POA.(fld)(:,1)-avg_f2Traj_POA.(fld)(:,2)./sqrt(avg_f2Traj_POA.(fld)(:,3)),'--','Color',colors.(fld));
    plot(taxis1,avg_f2Traj_POA.(fld)(:,1)+avg_f2Traj_POA.(fld)(:,2)./sqrt(avg_f2Traj_POA.(fld)(:,3)),'--','Color',colors.(fld));    
end
xlabel('Time (sec)'); ylabel('F2 (Hz)');


%% Compute POA_TN (POA, time-normalized) trajectories
if isfield(avg_f2Traj_POA,'none')
    fld1='none';
else
    fld1='noPert';
end
t_traj_f2=avg_f2Traj_POA.(fld1)(:,1);
diff_t_traj_f2=diff(t_traj_f2);
idx_turn=find(diff_t_traj_f2(1:end-1).*diff_t_traj_f2(2:end)<0);
f2LB_POA=t_traj_f2(idx_turn(1)+1);
f2UB_POA=t_traj_f2(idx_turn(2)+1);

d_t_traj_f2=diff(t_traj_f2);
idx_peak_1=find(d_t_traj_f2(1:end-1)>0 & d_t_traj_f2(2:end)<0);

idxNormInt_POA=idx_peak_1(1)+1;

t_traj_f1=avg_f1Traj_POA.(fld1)(:,1);
diff_t_traj_f1=diff(t_traj_f1);
idx_maxima=find(diff_t_traj_f1(1:end-1).*diff_t_traj_f1(2:end)<0 & diff_t_traj_f1(1:end-1)>0);
idx_maxima=idx_maxima(1);            
idx_minima=find(diff_t_traj_f1(1:end-1).*diff_t_traj_f1(2:end)<0 & diff_t_traj_f1(1:end-1)<0);
idx_minima=idx_minima(idx_minima>idx_maxima);            
f1LB_POA=t_traj_f1(idx_minima(1)+1);
f1UB_POA=t_traj_f1(idx_maxima+1);

trajF2_POA_TN=protoStruct_cell;
trajF1_POA_TN=protoStruct_cell;
trajF2_POA_TN_FN=protoStruct_cell;
trajF1_POA_TN_FN=protoStruct_cell;

flds=fields(trajF1_POA);
for i1=1:length(flds)
    fld=flds{i1};
    for i2=1:length(trajF2_POA.(fld))
        t_ind_traj_f2=trajF2_POA.(fld){i2};
        t_ind_traj_f1=trajF1_POA.(fld){i2};

        fNormF2=(t_ind_traj_f2-f2LB_POA)/(f2UB_POA-f2LB_POA);
        fNormF1=(t_ind_traj_f1-f1LB_POA)/(f1UB_POA-f1LB_POA);

        oldTAxis=1:length(t_ind_traj_f2);
        newTAxis=linspace(1,idxNormInt_POA,400);
        dt=newTAxis(2)-newTAxis(1);
        newTAxis=[newTAxis,newTAxis(end)+dt:dt:oldTAxis(end)];
        trajF2_POA_TN.(fld){end+1}=interp1(oldTAxis,t_ind_traj_f2,newTAxis);
        trajF1_POA_TN.(fld){end+1}=interp1(oldTAxis,t_ind_traj_f1,newTAxis);
        
        trajF2_POA_TN_FN.(fld){end+1}=interp1(oldTAxis,fNormF2,newTAxis);
        trajF1_POA_TN_FN.(fld){end+1}=interp1(oldTAxis,fNormF1,newTAxis);
    end
    
    avg_f2Traj_POA_TN.(fld)=avgTrace1(trajF2_POA_TN.(fld));
    avg_f1Traj_POA_TN.(fld)=avgTrace1(trajF1_POA_TN.(fld));
    avg_f2Traj_POA_TN_FN.(fld)=avgTrace1(trajF2_POA_TN_FN.(fld));
    avg_f1Traj_POA_TN_FN.(fld)=avgTrace1(trajF1_POA_TN_FN.(fld));
end



%% Add stage2 data to pdata
pdata.stage2=struct;
% pdata.isEMMA=isEMMA;

pdata.stage2.timeStamp  = clock;
pdata.stage2.nInterp    = nInterp;
pdata.stage2.normT1_uTime           = normT1_uTime;
pdata.stage2.normT1_bw_iu_ratio     = config.normT1_bw_iu_ratio;
pdata.stage2.normT1_bw_uy_ratio     = config.normT1_bw_uy_ratio;
pdata.stage2.normT1_bw_iu_f2        = normT1_bw_iu_f2;
pdata.stage2.normT1_bw_uy_f2        = normT1_bw_uy_f2;
pdata.stage2.normT2_bw_uy_ratio     = config.normT2_bw_uy_ratio;
pdata.stage2.normT2_bw_uy_f2        = normT2_bw_uy_f2;
pdata.stage2.normT3_bw_yu_ratio     = config.normT3_bw_yu_ratio;
pdata.stage2.normT3_bw_yu_f2        = normT3_bw_yu_f2;

pdata.stage2.avg_f1Traj=avg_f1Traj;
pdata.stage2.avg_f2Traj=avg_f2Traj;
pdata.stage2.avg_f1Traj_POA=avg_f1Traj_POA;
pdata.stage2.avg_f2Traj_POA=avg_f2Traj_POA;

pdata.stage2.avg_rmsTraj = avg_rmsTraj;

pdata.stage2.avg_f1Traj_pert=avg_f1Traj_pert;
pdata.stage2.avg_f2Traj_pert=avg_f2Traj_pert;
pdata.stage2.avg_f1_pertShift=avg_f1_pertShift;
pdata.stage2.avg_f2_pertShift=avg_f2_pertShift;
pdata.stage2.avg_f1Traj_IOA_TN=avg_f1Traj_IOA_TN;
pdata.stage2.avg_f2Traj_IOA_TN=avg_f2Traj_IOA_TN;
pdata.stage2.avg_f1Traj_IOA_FN=avg_f1Traj_IOA_FN;
pdata.stage2.avg_f2Traj_IOA_FN=avg_f2Traj_IOA_FN;
pdata.stage2.avg_f1Traj_IOA_TN_FN=avg_f1Traj_IOA_TN_FN;
pdata.stage2.avg_f2Traj_IOA_TN_FN=avg_f2Traj_IOA_TN_FN;

pdata.stage2.avg_f1Traj_POA_TN=avg_f1Traj_POA_TN;
pdata.stage2.avg_f2Traj_POA_TN=avg_f2Traj_POA_TN;
pdata.stage2.avg_f1Traj_POA_TN_FN=avg_f1Traj_POA_TN_FN;
pdata.stage2.avg_f2Traj_POA_TN_FN=avg_f2Traj_POA_TN_FN;

pdata.stage2.avg_f1Traj_pert_IOA_TN=avg_f1Traj_pert_IOA_TN;
pdata.stage2.avg_f2Traj_pert_IOA_TN=avg_f2Traj_pert_IOA_TN;
pdata.stage2.avg_f1Traj_pert_IOA_FN=avg_f1Traj_pert_IOA_FN;
pdata.stage2.avg_f2Traj_pert_IOA_FN=avg_f2Traj_pert_IOA_FN;
pdata.stage2.avg_f1Traj_pert_IOA_TN_FN=avg_f1Traj_pert_IOA_TN_FN;
pdata.stage2.avg_f2Traj_pert_IOA_TN_FN=avg_f2Traj_pert_IOA_TN_FN;

pdata.stage2.trajF1_POA=trajF1_POA;
pdata.stage2.trajF2_POA=trajF2_POA;

pdata.stage2.avg_f1_pertShift_IOA_TN=avg_f1_pertShift_IOA_TN;
pdata.stage2.avg_f1_pertShift_IOA_FN=avg_f1_pertShift_IOA_FN;
pdata.stage2.avg_f1_pertShift_IOA_TN_FN=avg_f1_pertShift_IOA_TN_FN;
pdata.stage2.avg_f2_pertShift_IOA_TN=avg_f2_pertShift_IOA_TN;
pdata.stage2.avg_f2_pertShift_IOA_FN=avg_f2_pertShift_IOA_FN;
pdata.stage2.avg_f2_pertShift_IOA_TN_FN=avg_f2_pertShift_IOA_TN_FN;

pdata.stage2.fullTNormF2_tAxis=fullTNorm_tAxis;  % For backward compatibility
pdata.stage2.fullTNorm_tAxis=fullTNorm_tAxis;
pdata.stage2.fullTNormF2=fullTNormF2;

pdata.stage2.avg_fullTNormF2=avg_fullTNormF2;
pdata.stage2.avg_fullTNormF2_nn=avg_fullTNormF2_nn;
pdata.stage2.avg_fullTNormF2_cs=avg_fullTNormF2_cs;

pdata.stage2.fullTNormPertShiftF2=fullTNormPertShiftF2;
pdata.stage2.avg_fullTNormPertShiftF2=avg_fullTNormPertShiftF2;

pdata.stage2.fullTNormF2_pert_FN=fullTNormF2_pert_FN;
pdata.stage2.avg_fullTNormF2_pert_FN=avg_fullTNormF2_pert_FN;

pdata.stage2.fullTNormF2_FN=fullTNormF2_FN;
pdata.stage2.avg_fullTNormF2_FN=avg_fullTNormF2_FN;

pdata.stage2.fullTNormF2_AT=fullTNormF2_AT;
pdata.stage2.avg_fullTNormF2_AT=avg_fullTNormF2_AT;

pdata.stage2.timeShifts=timeShifts;

% Some clean-up
fields_to_clean={'avg_f1Traj','avg_f2Traj','avg_f1Traj_POA','avg_f2Traj_POA'};
for i1=1:length(fields_to_clean)
    if isfield(pdata,fields_to_clean{i1})
        pdata=rmfield(pdata,fields_to_clean{i1});
        fprintf('Clean-up: removing pdata.%s\n',fields_to_clean{i1});
    end
end

if bSave
    save(dacacheFN,'pdata');
    fprintf('pdata file: %s updated.\n',dacacheFN);
else
    fprintf('*** bSave==0. Skipping the pdata update. ***\n');
end

if ~isempty(fsic(varargin, 'saveOnly'))
    return
end

%% Time normalized F2 trajectories
metaTracePlot_(segTrajF2_tNorm1,colors,'segTrajF2_tNorm1: anchored between iouOnset and youOnset','F2 (Hz)','Time (normalized)');
metaTracePlot_(segTrajF2_tNorm2,colors,'segTrajF2_tNorm2: anchored between uTime and youOnset','F2 (Hz)','Time (normalized)');
metaTracePlot_(segTrajF2_tNorm3,colors,'segTrajF2_tNorm3: anchored between youOnset and uayOnset','F2 (Hz)','Time (normalized)');

%% Time normalized F2 velocities
% metaTracePlot_(segVelF2_tNorm1,colors,'segVelF2_tNorm1: anchored between uTime and youOnset','F2 velocity (Hz/sec)','Time (normalized)');

%% Adaptation test
% metaTracePlot_(adaptTestData.segTrajF2_tNorm1,fColors,'adaptTestData.segTrajF2_tNorm1: anchored between iouOnset and youOnset','F2 (Hz)','Time (normalized)');

%%
metaPlot_errorBar(iuInt,'i-u interval','sec',fontSize);
metaPlot_errorBar(iyInt,'i-y interval','sec',fontSize);
metaPlot_errorBar(uyInt,'u-y interval','sec',fontSize);
metaPlot_errorBar(iu2Int,'i-u2 interval','sec',fontSize);

metaPlot_errorBar(i_uyMid_int,'i-[u,y] interval','sec',fontSize);
metaPlot_errorBar(i_yu2Mid_int,'i-[y,u2] interval','sec',fontSize);

metaPlot_errorBar(i_uyMidF_int,'i-[u,y]F interval','sec',fontSize);
metaPlot_errorBar(i_yu2MidF_int,'i-[y,u2]F interval','sec',fontSize);

metaPlot_errorBar(iF2,'F2 of i','Hz',fontSize);
metaPlot_errorBar(uF2,'F2 of [u]','Hz',fontSize);
metaPlot_errorBar(yF2,'F2 of [j]_1','Hz',fontSize);
metaPlot_errorBar(uF2_you,'F2 u in "you"','Hz',fontSize);
metaPlot_errorBar(iuF2Dec,'i-u F2 decrease','Hz',fontSize);
metaPlot_errorBar(uyF2Inc,'u-y F2 increase','Hz',fontSize);
metaPlot_errorBar(youF2Dec,'y-u F2 decrease in "you"','Hz',fontSize);
metaPlot_errorBar(iuF2Vel,'i-u F2 velocity','Hz/s',fontSize);
metaPlot_errorBar(uyF2Vel,'u-y F2 velocity','Hz/s',fontSize);

for i1=1:length(config.normT1_bw_iu_ratio)
    metaPlot_errorBar(normT1_bw_iu_f2{i1},sprintf('normT1: F2 @ %.2f point between i and u',config.normT1_bw_iu_ratio(i1)),'Hz',fontSize);
end
for i1=1:length(config.normT1_bw_uy_ratio)
    metaPlot_errorBar(normT1_bw_uy_f2{i1},sprintf('normT1: F2 @ %.2f point between u and y',config.normT1_bw_uy_ratio(i1)),'Hz',fontSize);
end
for i1=1:length(config.normT2_bw_uy_ratio)
    metaPlot_errorBar(normT2_bw_uy_f2{i1},sprintf('normT2: F2 @ %.2f point between u and y',config.normT2_bw_uy_ratio(i1)),'Hz',fontSize);
end
for i1=1:length(config.normT3_bw_yu_ratio)
    metaPlot_errorBar(normT3_bw_yu_f2{i1},sprintf('normT3: F2 @ %.2f point between y and u',config.normT3_bw_yu_ratio(i1)),'Hz',fontSize);
end

%% Output
if ~isempty(fsic(varargin,'output'))
    
    output_str = varargin{fsic(varargin,'output') + 1};
    for i1 = 1 : numel(output_str);
        eval(sprintf('varargout{%d} = %s;', i1, output_str{i1}))
    end
end
return

%% Sub-routines

function raw_fn=getRawFN_(expDir,fn)
[path1,fn1]=fileparts(fn);
[path2,fn2]=fileparts(path1);
[path3,fn3]=fileparts(path2);
[path4,fn4]=fileparts(path3);

raw_fn=fullfile(expDir,fn3,fn2,fn1);
if ~isequal(raw_fn(end-3:end),'.mat')
    raw_fn=[raw_fn,'.mat'];
end

%%
function pertType=getPertType_data(data)
if data.params.bShift==0
    pertType='none';
else
    if isfield(data.params,'bTimeWarp')
        if data.params.bTimeWarp==0
            pertType=lower(data.params.pertField.pertDirection);
        else
            if data.params.bDecelWarp==1
                pertType='decel';
            else
                pertType='accel';
            end
        end
    else
        pertType=lower(data.params.pertField.pertDirection);
    end
end

return

%%
function metaTracePlot_(trajs,colors,name,y_label,x_label,varargin)
fontSize=14;

if isempty(fsic(varargin,'Position'))
    figure('Position',[0,-400,1000,800],'NumberTitle','off','Name',name);
else    
    figure('Position',varargin{fsic(varargin,'Position')+1},'NumberTitle','off','Name',name);
end
subplot('Position',[0.075,0.5,0.9,0.45]);
set(gca,'FontSize',fontSize);
frameDur=16/12e3;
tMin=Inf; tMax=-Inf;
for i1=1:length(fields(trajs))
	flds=fields(trajs);
	fld=flds{i1};
    avg_traj.(fld)=avgTrace1(trajs.(fld));
	for i2=1:length(trajs.(fld))
		t_f2=trajs.(fld){i2};       
        if isempty(strfind(lower(x_label),'normalized'))
    		taxis1=0:(frameDur):(frameDur*(length(t_f2)-1));
        else
            taxis1=linspace(0,1,length(t_f2));
        end
		plot(taxis1,t_f2,'Color',colors.(fld)); 
		hold on;
        if min(taxis1)<tMin
            tMin=min(taxis1);
        end
        if max(taxis1)>tMax
            tMax=max(taxis1);
        end
	end
end
set(gca,'XTick',[]);
set(gca,'XLim',[tMin,tMax]);
ylim=get(gca,'YLim');

ylabel(y_label,'FontSize',fontSize);

subplot('Position',[0.075,0.05,0.9,0.45]);
set(gca,'FontSize',fontSize);
for i1=1:length(fields(avg_traj))
	flds=fields(avg_traj);
	fld=flds{i1};    
    if isempty(strfind(lower(x_label),'normalized'))        
        taxis1=0:(frameDur):(frameDur*(length(avg_traj.(fld)(:,1))-1));
    else
        taxis1=linspace(0,1,length(avg_traj.(fld)(:,1)));
    end
    plot(taxis1,avg_traj.(fld)(:,1),'-','Color',colors.(fld)); hold on;
    plot(taxis1,avg_traj.(fld)(:,1)-avg_traj.(fld)(:,2)./sqrt(avg_traj.(fld)(:,3)),'--','Color',colors.(fld));
    plot(taxis1,avg_traj.(fld)(:,1)+avg_traj.(fld)(:,2)./sqrt(avg_traj.(fld)(:,3)),'--','Color',colors.(fld));    
end
set(gca,'XLim',[tMin,tMax],'YLim',ylim);

xlabel(x_label,'FontSize',fontSize);
ylabel(y_label,'FontSize',fontSize);





return

%%
function metaPlot_errorBar(meas,measName,measUnit,fontSize)
flds=fields(meas);
x_tick_label=cell(1,0);
meas_mean=[];
meas_sem=[];
for i1=1:length(flds)
	t_meas=meas.(flds{i1});
	x_tick_label{length(x_tick_label)+1}=flds{i1};
	meas_mean=[meas_mean,mean(t_meas(~isnan(t_meas)))];
	meas_sem=[meas_sem,ste(t_meas(~isnan(t_meas)))];
end
figure('NumberTitle','off','name',measName,'Position',[200,200,320,240]);
set(gca,'FontSize',fontSize);
errorbar([1:length(flds)],meas_mean,meas_sem,'LineWidth',2);
set(gca,'XLim',[0.5,length(flds)+0.5],'XTick',[1:length(flds)]);
set(gca,'XTickLabel',x_tick_label);
ylabel([measName,' (',measUnit,')']);
title(measName,'FontSize',12);

if isfield(meas,'noPert')
    fldCmp0='noPert';
else
    fldCmp0='none';
end
if isfield(meas,'down') && isfield(meas,'up')
    fldCmp1='down';
    fldCmp2='up';
else
    fldCmp1='accel';
    fldCmp2='decel';
end
[p_rs_12,h]=ranksum(meas.(fldCmp1),meas.(fldCmp2));
[p_rs_01,h]=ranksum(meas.(fldCmp0),meas.(fldCmp1));
[p_rs_02,h]=ranksum(meas.(fldCmp0),meas.(fldCmp2));


fprintf('%s: \n',measName)
fprintf('\t%s vs. %s: p = %.5f\n',fldCmp1,fldCmp2,p_rs_12);
fprintf('\t%s vs. %s: p = %.5f\n',fldCmp0,fldCmp1,p_rs_01);
fprintf('\t%s vs. %s: p = %.5f\n',fldCmp0,fldCmp2,p_rs_02);
return