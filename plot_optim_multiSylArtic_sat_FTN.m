function plot_optim_multiSylArtic_sat_FTN(simDataSets, optimFld, minMode, varargin)
%% Config
frameDur=16/12e3;   % s
colors.none=[0,0,0];
colors.down=[0,0,1];
colors.up=[1,0,0];
colors.accel=[1,0.5,0];
colors.decel=[0,0.5,0];

colorsShade=struct;
flds=fields(colors);
for i1=1:numel(flds)
    fld=flds{i1};
    colorsShade.(fld)=1-0.5*(1-(colors.(fld)));
end

fs=14;
resLW=2;
barW=0.04;

%%
ds_prefix = 'sat_multiSylArtic_velSimp_FTN_simDataSet';

if isnumeric(simDataSets)
    nums = simDataSets;
    
    simDataSets = {};
    for i1 = 1 : numel(nums)
        simDataSets{i1} = [ds_prefix, sprintf('_%.3d', nums(i1))];
    end    
end

if iscell(simDataSets)
    nSets=numel(simDataSets);
else
    nSets=1;
    simDataSets={simDataSets};
end


a.TIMING_UPDATE_COEFF_WS=[];
a.BS_WS_TEMPORAL_RATIO=[];
a.IM_AUD_DELAY=[];
a.UPDATE_INTERVAL=[];
a.FB_GAIN=[];
a.REL_ERR_THRESH=[];

res.F2_err=[];
res.F2_err_up=[];
res.F2_err_down=[];
res.F2_err_upDown=[];
res.F2_err_accel=[];
res.F2_err_decel=[];
res.F2_err_accelDecel=[];
res.tInt_err=[];
res.tInt_err_upDown=[];
res.tInt_err_accelDecel=[];

for i1=1:nSets
    data=load(simDataSets{i1});    
    
    flds=fields(a);
    for i2=1:numel(flds)
        fld=flds{i2};
        if isfield(data,['a_',fld])
            eval(['a.',fld,'=[a.',fld,',data.a_',fld,'];']);
        else
            if isequal(fld,'FB_GAIN')
                a.FB_GAIN=[a.FB_GAIN,zeros(size(data.a_TIMING_UPDATE_COEFF_WS))];
            end
        end
    end
    
    res.F2_err = [res.F2_err, data.a_F2_err];
    res.F2_err_up = [res.F2_err_up, data.a_F2_err_up];
    res.F2_err_down = [res.F2_err_down, data.a_F2_err_down];
    res.F2_err_upDown = [res.F2_err_upDown, data.a_F2_err_upDown];
    res.F2_err_accel = [res.F2_err_accel, data.a_F2_err_accel];
    res.F2_err_decel = [res.F2_err_decel, data.a_F2_err_decel];
    res.F2_err_accelDecel = [res.F2_err_accelDecel, data.a_F2_err_accelDecel];
    res.tInt_err = [res.tInt_err, data.a_tInt_err];
    res.tInt_err_upDown = [res.tInt_err_upDown, data.a_tInt_err_upDown];
    res.tInt_err_accelDecel = [res.tInt_err_accelDecel, data.a_tInt_err_accelDecel];
    
    if ~exist('targF2_FTN')
        targF2_FTN=data.targF2_FTN;
        if size(targF2_FTN.upDown,1)>size(targF2_FTN.upDown,2)
            targF2_FTN.upDown=targF2_FTN.upDown';
            targF2_FTN.accelDecel=targF2_FTN.accelDecel';            
        end
        F2UB=data.F2UB;
        T_STEP=data.T_STEP;
    end
end



%% Load positional fb control data sets, if specified
if ~isempty(fsic(varargin,'posFB')) || ~isempty(fsic(varargin,'velFB'))
    bPos=1;
    if ~isempty(fsic(varargin,'posFB'))
        posDataSets=varargin{fsic(varargin,'posFB')+1};
    else
        posDataSets=varargin{fsic(varargin,'velFB')+1};
    end
    if iscell(posDataSets)
        nSetsP=numel(posDataSets);
    else
        nSetsP=1;
        simDataSets={posDataSets};
    end

    ap.fbDelay=[];
    ap.fbGain=[];
    
    resP.F2_err=[];
    resP.F2_err_upDown=[];
    resP.F2_err_accelDecel=[];
    resP.tInt_err=[];
    resP.tInt_err_upDown=[];
    resP.tInt_err_accelDecel=[];
    
    for i1=1:nSetsP
        dataP=load(posDataSets{i1});    

        flds=fields(ap);
        for i2=1:numel(flds)
            fld=flds{i2};
            if isfield(dataP,['a_',fld])
                eval(['ap.',fld,'=[ap.',fld,',dataP.a_',fld,'];']);
            else
                pause(0);
            end
        end

        resP.F2_err = [resP.F2_err, dataP.a_F2_err];
        resP.F2_err_upDown = [resP.F2_err_upDown, dataP.a_F2_err_upDown];
        resP.F2_err_accelDecel = [resP.F2_err_accelDecel, dataP.a_F2_err_accelDecel];
        resP.tInt_err = [resP.tInt_err, dataP.a_tInt_err];
        resP.tInt_err_upDown = [resP.tInt_err_upDown, dataP.a_tInt_err_upDown];
        resP.tInt_err_accelDecel = [resP.tInt_err_accelDecel, dataP.a_tInt_err_accelDecel];
    end
else
    bPos = 0;
end

% [foo, idx] = sort(res.F2_err_upDown);
% res.order_F2_err_upDown = nan(size(res.F2_err_upDown));
% for i1 = 1 : numel(idx)
%     res.order_F2_err_upDown(i1) = find(res.F2_err_upDown == res.F2_err_upDown(idx(i1)));
% end
res.norm_F2_err_upDown = (res.F2_err_upDown - min(res.F2_err_upDown)) / range(res.F2_err_upDown);
res.norm_tInt_err_upDown = (res.tInt_err_upDown - min(res.tInt_err_upDown)) / range(res.tInt_err_upDown);
res.compNorm_err_upDown = sqrt(res.norm_F2_err_upDown.^2 + res.norm_tInt_err_upDown.^2);

res.norm_F2_err_accelDecel = (res.F2_err_accelDecel - min(res.F2_err_accelDecel)) / range(res.F2_err_accelDecel);
res.norm_tInt_err_accelDecel = (res.tInt_err_accelDecel - min(res.tInt_err_accelDecel)) / range(res.tInt_err_accelDecel);
res.compNorm_err_accelDecel = sqrt(res.norm_F2_err_accelDecel.^2 + res.norm_tInt_err_accelDecel.^2);

res.compNorm_err_total = sqrt(res.norm_F2_err_upDown.^2 + res.norm_tInt_err_upDown.^2 + res.norm_F2_err_accelDecel.^2 + res.norm_tInt_err_accelDecel.^2);

%% Parametric constraints
PARAMS = fields(a);
resFlds = fields(res);

for i1=1:numel(PARAMS)
    par = PARAMS{i1};
    if ~isempty(fsic(varargin,par))
        parConstr = varargin{fsic(varargin,par)+1};
        if numel(parConstr) == 1
            idx_keep = find(a.(par)==parConstr);
            for j1=1:numel(PARAMS)
                a.(PARAMS{j1}) = a.(PARAMS{j1})(idx_keep);
            end
            for j1=1:numel(resFlds)
                if ~isempty(res.(resFlds{j1}))
                    res.(resFlds{j1}) = res.(resFlds{j1})(idx_keep);
                end
            end
        elseif numel(parConstr) ==2
            idx_keep = find(a.(par)>=parConstr(1) && a.(par)<=parConstr(2));
            for j1=1:numel(PARAMS)
                a.(PARAMS{j1}) = a.(PARAMS{j1})(idx_keep);
            end
            for j1=1:numel(resFlds)
                if ~isempty(res.(resFlds{j1}))
                    res.(resFlds{j1}) = res.(resFlds{j1})(idx_keep);
                end
            end
        end
    end
end


%%

% optim_flds={'F2Err','tIntErr'};
% optim_flds={'F2Err'};
optim_flds = {optimFld};

for i1=1:numel(optim_flds)
    optimFld=optim_flds{i1};

    if isequal(optimFld,'F2Err')
        if isequal(minMode,'total')
            [min_err,idx_min]=min(res.F2_err);
            if bPos
                [min_err_p,idx_min_p]=min(resP.F2_err);
            end
        elseif isequal(minMode,'upDown')
            [min_err,idx_min]=min(res.F2_err_upDown);
%             min_err = quantile(res.F2_err_upDown,0.2);
%             idx_min = find(res.F2_err_upDown == min_err);
%             idx_min = idx_min(1);
            if bPos
                [min_err_p,idx_min_p]=min(resP.F2_err_upDown);
            end
        elseif isequal(minMode, 'up')
            [min_err,idx_min]=min(res.F2_err_up);
        elseif isequal(minMode, 'down')
            [min_err,idx_min]=min(res.F2_err_down);
        elseif isequal(minMode,'accelDecel')
            [min_err,idx_min]=min(res.F2_err_accelDecel);
%             min_err = quantile(res.F2_err_upDown,0.2);
%             idx_min = find(res.F2_err_upDown == min_err);
%             idx_min = idx_min(1);
            if bPos
                [min_err_p,idx_min_p]=min(resP.F2_err_accelDecel);
            end
        elseif isequal(minMode, 'accel')
            [min_err, idx_min] = min(res.F2_err_accel);
        elseif isequal(minMode, 'decel')
            [min_err, idx_min] = min(res.F2_err_decel);
        end
    elseif isequal(optimFld,'tInt')
        if isequal(minMode,'total')
            [min_err,idx_min]=min(res.tInt_err);
            if bPos
                [min_err_p,idx_min_p]=min(resP.tInt_err);
            end
        elseif isequal(minMode,'upDown')
            [min_err,idx_min]=min(res.tInt_err_upDown);
            if bPos
                [min_err_p,idx_min_p]=min(resP.tInt_err_upDown);
            end
        elseif isequal(minMode,'accelDecel')
            [min_err,idx_min]=min(res.tInt_err_accelDecel);
            if bPos
                [min_err_p,idx_min_p]=min(resP.tInt_err_accelDecel);
            end
        end 
    elseif isequal(optimFld, 'compNorm')
        if isequal(minMode, 'upDown')
            [min_err,idx_min]=min(res.compNorm_err_upDown);
        elseif isequal(minMode, 'accelDecel')
            [min_err,idx_min]=min(res.compNorm_err_accelDecel);
        end
    end
    
    F2_optim.(optimFld).minErr=min_err;
    F2_optim.(optimFld).F2_err=res.F2_err(idx_min);
    F2_optim.(optimFld).tInt_err=res.tInt_err(idx_min);
%     F2_optim.(optimFld).FM_COMMAND_UPDATE_COEFF_WS=a.FM_COMMAND_UPDATE_COEFF_WS(idx_min);
    F2_optim.(optimFld).TIMING_UPDATE_COEFF_WS=a.TIMING_UPDATE_COEFF_WS(idx_min);
    F2_optim.(optimFld).BS_WS_TEMPORAL_RATIO=a.BS_WS_TEMPORAL_RATIO(idx_min);
%     F2_optim.(optimFld).SPATIAL_UPDATE_LEAD_N=a.SPATIAL_UPDATE_LEAD_N(idx_min);
    F2_optim.(optimFld).IM_AUD_DELAY=a.IM_AUD_DELAY(idx_min);
    F2_optim.(optimFld).UPDATE_INTERVAL=a.UPDATE_INTERVAL(idx_min);
    F2_optim.(optimFld).FB_GAIN=a.FB_GAIN(idx_min);
    F2_optim.(optimFld).REL_ERR_THRESH = a.REL_ERR_THRESH(idx_min);
    
    if bPos
        F2_optim_p.(optimFld).minErr=min_err_p;
        F2_optim_p.(optimFld).F2_err=resP.F2_err(idx_min_p);
        F2_optim_p.(optimFld).tInt_err=resP.tInt_err(idx_min_p);
        F2_optim_p.(optimFld).fbDelay=ap.fbDelay(idx_min_p);
        F2_optim_p.(optimFld).fbGain=ap.fbGain(idx_min_p);
    end

    t_text{1} = sprintf('(%s) Optim parameter setting found: F2Err = %.2f Hz; tIntErr = %.6f ms',optimFld,res.F2_err(idx_min),res.tInt_err(idx_min)*1e3);
%     t_text{2} = sprintf('\tFM_COMMAND_UPDATE_COEFF_WS = %f',F2_optim.(optimFld).FM_COMMAND_UPDATE_COEFF_WS);
    t_text{2} = sprintf('\tTIMING_UPDATE_COEFF_WS = %f',F2_optim.(optimFld).TIMING_UPDATE_COEFF_WS);
    t_text{3} = sprintf('\tBS_WS_TEMPORAL_RATIO = %f',F2_optim.(optimFld).BS_WS_TEMPORAL_RATIO);
%     t_text{5} = sprintf('\tSPATIAL_UPDATE_LEAD_N = %f',F2_optim.(optimFld).SPATIAL_UPDATE_LEAD_N);
    t_text{4} = sprintf('\tIM_AUD_DELAY = %f',F2_optim.(optimFld).IM_AUD_DELAY);
    t_text{5} = sprintf('\tUPDATE_INTERVAL = %f',F2_optim.(optimFld).UPDATE_INTERVAL);
    t_text{6} = sprintf('\tFB_GAIN = %f',F2_optim.(optimFld).FB_GAIN);
    t_text{7} = sprintf('\tREL_ERR_THRESH = %f', F2_optim.(optimFld).REL_ERR_THRESH);
    
    for j1=1:numel(t_text)
        fprintf('%s\n',t_text{j1});
    end
    
    if bPos
        p_text{1} = sprintf('(%s) Optim parameter setting (** FOR POSITIONAL FB CONTROL **) found: F2Err = %.2f Hz; tIntErr = %.6f ms',...
            optimFld,resP.F2_err(idx_min_p),resP.tInt_err(idx_min_p)*1e3);
        p_text{2} = sprintf('\tfbDelay = %f',F2_optim_p.(optimFld).fbDelay);
        p_text{3} = sprintf('\tfbGain = %f',F2_optim_p.(optimFld).fbGain);
        
        for j1=1:numel(p_text)
            fprintf('%s\n',p_text{j1});
        end
    end
    
    for i2=1:2
        if i2==1
            fld1='down';
            fld2='up';
            mode='upDown';
        else
            fld1='accel';
            fld2='decel';
            mode='accelDecel';
        end
        
        if isequal(mode,'upDown')
%             [s_prodF2,s_audF2,OVERSAMPLE_FACT,T_STEP,chg_iuInt,chg_iyInt]=multiSylArtic_velSimp_model1('PERT_TYPES',{fld1,fld2},...
%                 'targF2',targF2_FTN.(mode),T_STEP,'F2UB',F2UB.(mode),...                
%                 'TIMING_UPDATE_COEFF.WS',F2_optim.(optimFld).TIMING_UPDATE_COEFF_WS,...
%                 'BS_WS_TEMPORAL_RATIO',F2_optim.(optimFld).BS_WS_TEMPORAL_RATIO,...                
%                 'IM_AUD_DELAY',F2_optim.(optimFld).IM_AUD_DELAY,...
%                 'UPDATE_INTERVAL',F2_optim.(optimFld).UPDATE_INTERVAL,...
%                 'FB_GAIN',F2_optim.(optimFld).FB_GAIN,...
%                 'UP_PERT_COEFF',data.optimPertCoeff_FTN.up,...
%                 'DOWN_PERT_COEFF',data.optimPertCoeff_FTN.down,...
%                 'OVERSAMPLE_FACT',1,'noPrompt','noPlot');
            
            [s_prodF2, s_audF2, OVERSAMPLE_FACT, T_STEP, rec_idx, idx_0]=sat_velSimp_model1('PERT_TYPES',{fld1, fld2},...
                'kfMap',data.kfMap.(mode),...
                'targF2',data.targF2_FTN.(mode),data.T_STEP,'F2UB',data.F2UB.(mode),...
                'TIMING_UPDATE_COEFF.WS',F2_optim.(optimFld).TIMING_UPDATE_COEFF_WS,...
                'BS_WS_TEMPORAL_RATIO',F2_optim.(optimFld).BS_WS_TEMPORAL_RATIO,...
                'IM_AUD_DELAY',F2_optim.(optimFld).IM_AUD_DELAY,...
                'UPDATE_INTERVAL',F2_optim.(optimFld).UPDATE_INTERVAL,...
                'FB_GAIN',F2_optim.(optimFld).FB_GAIN,...
                'REL_ERR_THRESH',F2_optim.(optimFld).REL_ERR_THRESH,...
                'OVERSAMPLE_FACT',1,...
                'UP_PERT_COEFF',data.optimPertCoeff_FTN.up,...
                'DOWN_PERT_COEFF',data.optimPertCoeff_FTN.down,...
                'noPlot','noPrompt','B_MEX',1);
        elseif isequal(mode,'accelDecel')
            [s_prodF2, s_audF2, OVERSAMPLE_FACT, T_STEP, rec_idx, idx_0]=sat_velSimp_model1('PERT_TYPES',{fld1, fld2},...
                'kfMap',data.kfMap.(mode),...
                'targF2',data.targF2_FTN.(mode),data.T_STEP,'F2UB',data.F2UB.(mode),...
                'TIMING_UPDATE_COEFF.WS',F2_optim.(optimFld).TIMING_UPDATE_COEFF_WS,...
                'BS_WS_TEMPORAL_RATIO',F2_optim.(optimFld).BS_WS_TEMPORAL_RATIO,...
                'IM_AUD_DELAY',F2_optim.(optimFld).IM_AUD_DELAY,...
                'UPDATE_INTERVAL',F2_optim.(optimFld).UPDATE_INTERVAL,...
                'FB_GAIN',F2_optim.(optimFld).FB_GAIN,...
                'REL_ERR_THRESH',F2_optim.(optimFld).REL_ERR_THRESH,...
                'OVERSAMPLE_FACT',1,...
                'ACCEL_PERT_COEFF',data.optimPertCoeff_FTN.accel,...
                'DECEL_PERT_COEFF',data.optimPertCoeff_FTN.decel,...
                'noPlot','noPrompt','B_MEX',1);
        end
        
        if bPos
            d_prodF2_FTN=dataP.exptData.d_prodF2_FTN;
            if isequal(mode,'upDown')
                if ~isempty(fsic(varargin,'posFB'))
                    [err_p.(mode), sChg_iuInt_p.(mode), sChg_iyInt_p.(mode), sChg_prodF2_p.(mode)]=simErr_posFBControl(...
                            mode,[F2_optim_p.(optimFld).fbDelay; F2_optim_p.(optimFld).fbGain],...
                            targF2_FTN.(mode)',F2UB.(mode),T_STEP,d_prodF2_FTN.(mode),'total',[data.optimPertCoeff_FTN.down,data.optimPertCoeff_FTN.up],0);
                else
                    [err_p.(mode), sChg_iuInt_p.(mode), sChg_iyInt_p.(mode), sChg_prodF2_p.(mode)]=simErr_velFBControl(...
                            mode,[F2_optim_p.(optimFld).fbDelay; F2_optim_p.(optimFld).fbGain],...
                            targF2_FTN.(mode)',F2UB.(mode),T_STEP,d_prodF2_FTN.(mode),'total',[data.optimPertCoeff_FTN.down,data.optimPertCoeff_FTN.up],0);
                end
            elseif isequal(mode,'accelDecel')
                if ~isempty(fsic(varargin,'posFB'))
                    [err_p.(mode), sChg_iuInt_p.(mode), sChg_iyInt_p.(mode), sChg_prodF2_p.(mode)]=simErr_posFBControl(...
                            mode,[F2_optim_p.(optimFld).fbDelay; F2_optim_p.(optimFld).fbGain],...
                            targF2_FTN.(mode)',F2UB.(mode),T_STEP,d_prodF2_FTN.(mode),'total',[data.optimPertCoeff_FTN.accel,data.optimPertCoeff_FTN.decel],0);
                else
                    [err_p.(mode), sChg_iuInt_p.(mode), sChg_iyInt_p.(mode), sChg_prodF2_p.(mode)]=simErr_velFBControl(...
                            mode,[F2_optim_p.(optimFld).fbDelay; F2_optim_p.(optimFld).fbGain],...
                            targF2_FTN.(mode)',F2UB.(mode),T_STEP,d_prodF2_FTN.(mode),'total',[data.optimPertCoeff_FTN.accel,data.optimPertCoeff_FTN.decel],0);
                end
            end
        end
        
        %-- Visualization --%
        %- Plot anb compare the perturbations -%
        if i1==1    % taxis problematic!
%             figure;
%             set(gca,'FontSize',fs);
%             N_SKIP=2;
%             taxis_p=0 : frameDur : data.exptData.upDown.frameDur_FTN*(length(data.exptData.(mode).avgPertShiftF2_FTN.(fld1)(:,1))-1);
%             plot(taxis_p(1:N_SKIP:end),data.exptData.(mode).avgPertShiftF2_FTN.(fld1)(1:N_SKIP:end,1),'.','Color',colors.(fld1),'MarkerSize',4);
%             hold on;
%             taxis_p=0 : data.exptData.upDown.frameDur_FTN : data.exptData.upDown.frameDur_FTN*(length(data.exptData.(mode).avgPertShiftF2_FTN.(fld2)(:,1))-1);
%             plot(taxis_p(1:N_SKIP:end),data.exptData.(mode).avgPertShiftF2_FTN.(fld2)(1:N_SKIP:end,1),'.','Color',colors.(fld2),'MarkerSize',4);
% 
%             N_SKIP=15;
%             s_pertShift.(fld1)=s_audF2.(fld1)-s_prodF2.(fld1);
%             s_pertShift.(fld2)=s_audF2.(fld2)-s_prodF2.(fld2);
%             taxis_p = 0 : T_STEP : T_STEP*(length(s_pertShift.(fld1))-1);
%             plot(taxis_p(1:N_SKIP:end), s_pertShift.(fld1)(1:N_SKIP:end), '.', 'Color', colors.(fld1),'MarkerSize',5);
%             plot(taxis_p(1:N_SKIP:end), s_pertShift.(fld2)(1:N_SKIP:end), '.', 'Color', colors.(fld2),'MarkerSize',5);          
%             
%             set(gca,'XLim',[0,0.6]);
%             xlabel('Time re. [i] (s)');
%             ylabel('F2 perturbation (Hz)');
%             legend({['Expt. pert. (',fld1,')'], ['Expt. pert. (',fld2,')'], ['Sim. pert. (',fld1,')'], ['Sim. pert. (',fld2,')']});
           
        end
        
        %- Plot the production data -%
        figure('Position',[50,100,1200,720],'Color','w');
        h1=subplot('Position',[0.1,0.15,0.5,0.8]);
%         set(gca,'FontSize',fs);
%         [idx_peaks,idx_valleys]=find_peaks_valleys(data.exptData.(mode).avg_trajF2_IOA.none);
%         idx_j2=idx_peaks(2);

        taxis0=0 : 1 / 250 : 1 / 250 * (length(data.exptData.(mode).avgChgTrajF2_FTN.(fld1)) - 1);
%         taxis0=0:data.exptData.upDown.frameDur_FTN:data.exptData.upDown.frameDur_FTN*(idx_j2-1);
        plot([taxis0(1),taxis0(end)],[0,0],'-','Color',[0.5,0.5,0.5]);
        hold on;
        for k1 = 1 : 2
            if k1 == 1;  fld = fld1; else; fld = fld2; end
            plot_sd_t(taxis0, data.exptData.(mode).avgChgTrajF2_FTN.(fld)(:, 1),...
                data.exptData.(mode).avgChgTrajF2_FTN.(fld)(: ,2) ./ sqrt(data.exptData.(mode).avgChgTrajF2_FTN.(fld)(: ,3)),...
                colorsShade.(fld),'patch','FaceAlpha',0.2);
            plot(taxis0, data.exptData.(mode).avgChgTrajF2_FTN.(fld)(: ,1), '--', 'color', colors.(fld));
        end

        
        
        N = 250;
        [targF2_FTN, idxj_targ] = calc_FTN(data.targF2_FTN.(mode), N);
        [s_prodF2_FTN.(fld1), s_idxj.(fld1)] = calc_FTN(s_prodF2.(fld1), N);
        [s_prodF2_FTN.(fld2), s_idxj.(fld2)] = calc_FTN(s_prodF2.(fld2), N);        
        
        taxis1=linspace(taxis0(1), taxis0(end), length(s_prodF2_FTN.(fld1)));
        plot(taxis1, s_prodF2_FTN.(fld1) - targF2_FTN, '-', 'Color', colors.(fld1), 'LineWidth', resLW);
        plot(taxis1, s_prodF2_FTN.(fld2) - targF2_FTN, '-', 'Color', colors.(fld2), 'LineWidth', resLW);
        
        if bPos
            taxis1=0:T_STEP:T_STEP*(length(sChg_prodF2_p.(mode).(fld1))-1);
            plot(taxis1, sChg_prodF2_p.(mode).(fld1), '-', 'Color', colors.(fld1), 'LineWidth', resLW/2);
            plot(taxis1, sChg_prodF2_p.(mode).(fld2), '-', 'Color', colors.(fld2), 'LineWidth', resLW/2);
                        
            base_chg.(fld1) = interp1(taxis0,data.exptData.(mode).avgChgTrajF2_FTN.(fld1)(1:idx_j2,1),taxis1);
            base_chg.(fld2) = interp1(taxis0,data.exptData.(mode).avgChgTrajF2_FTN.(fld2)(1:idx_j2,1),taxis1);
            t_err = [s_prodF2.(fld1)-targF2_FTN.(mode)-base_chg.(fld1),s_prodF2.(fld2)-targF2_FTN.(mode)-base_chg.(fld2)];
            t_err = t_err(~isnan(t_err));
            t_F2_err_sqDIVA = rms(t_err);
            t_err = [sChg_prodF2_p.(mode).(fld1)'-base_chg.(fld1),sChg_prodF2_p.(mode).(fld2)'-base_chg.(fld2)];
            t_err = t_err(~isnan(t_err));
            t_F2_err_oldModel = rms(t_err);
            fprintf('\t * %s: \tsqDIVA F2 err = %.2f Hz\n',mode,t_F2_err_sqDIVA);
            fprintf('\t * %s: \toldModel F2 err = %.2f Hz\n',mode,t_F2_err_oldModel);
            
            axes('pos',[0.475,0.80,0.125,0.15]);
            bar(-0.2,t_F2_err_sqDIVA,0.4,'LineWidth',resLW,'FaceColor','none','EdgeColor','k');
            hold on;
            bar(0.2,t_F2_err_oldModel,0.4,'LineWidth',resLW/2,'FaceColor','none','EdgeColor','k');
            set(gca,'XLim',[-0.5,0.5]);
            set(gca,'XTick',[-0.2,0.2],'XTickLabel',{'sqDIVA','old model'});
            ylabel('Fitting error');
            set(gcf,'CurrentAxes',h1);
            
        end

        set(gca,'XLim',[taxis0(1),taxis0(end)]);        
        box on;
        xlabel('Time (FTN)','FontSize',fs);
        ylabel('F2 change from noPert','FontSize',fs);
        
        % Show text
        xs=get(gca,'XLim'); ys=get(gca,'YLim');
        for j1=1:numel(t_text)
            text(xs(1)+0.02*range(xs), ys(2)-j1*0.025*range(ys), strrep(t_text{j1},'_','\_'));
        end
        rectangle('Position',[xs(1),ys(1),range(xs),range(ys)],'FaceColor','none','EdgeColor','k');
        box off;
        
        if bPos==0
            ezlegend([xs(1)+0.5*range(xs),ys(1)+0.05*range(ys),0.475*range(xs),0.2*range(ys)],...
                0.3,{['Expt. data (',fld1,', \pm1 SEM)'], ['Simulation (',fld1,')'], ['Expt. data (',fld2,', \pm1 SEM)'], ['Simulation (',fld2,')']},...
                {colors.(fld1),colors.(fld1),colors.(fld2),colors.(fld2)},...
                [12,12,12,12],{'--','-','--','-'},{colors.(fld1),colors.(fld1),colors.(fld2),colors.(fld2)},[0.5,resLW,0.5,resLW],[0,0,0,0]);
        else
            ezlegend([xs(1)+0.5*range(xs),ys(1)+0.05*range(ys),0.475*range(xs),0.3*range(ys)],...
                0.3,{['Expt. data (',fld1,', \pm1 SEM)'], ['Simulation: sqDIVA (',fld1,')'], ['Simulation: old model (',fld1,')'], ...
                     ['Expt. data (',fld2,', \pm1 SEM)'], ['Simulation: sqDIVA (',fld2,')'], ['Simulation: old model (',fld2,')']},...
                {colors.(fld1),colors.(fld1),colors.(fld1),colors.(fld2),colors.(fld2),colors.(fld2)},...
                [12,12,12,12,12,12],{'--','-','-','--','-','-'},...
                {colors.(fld1),colors.(fld1),colors.(fld1),colors.(fld2),colors.(fld2),colors.(fld2)},...
                [0.5,resLW,resLW/2,0.5,resLW,resLW/2],[0,0,0,0,0,0]);
        end
        set(gca, 'XLim', [0, 5]);
        

        subplot('Position',[0.675,0.58,0.3,0.37]);
        plot(data.dChg_tInt.(mode)(:, 1) * 1e3, 's--', 'Color', colors.(fld1)); 
        hold on;
        plot(data.dChg_tInt.(mode)(:, 2) * 1e3, 's--', 'Color', colors.(fld2)); 
        set(gca, 'XLim', [0, size(data.dChg_tInt.(mode), 1) + 1]);
        xs = get(gca, 'XLim');
        plot(xs, [0, 0], '-', 'Color', [0.5, 0.5, 0.5]);
        
        t_idx_targ = [idxj_targ(2 : end - 1), idx_0(end)];
        t_idx.(fld1) = [s_idxj.(fld1)(2 : end - 1), rec_idx.(fld1)(end, end)];
        t_idx.(fld2) = [s_idxj.(fld2)(2 : end - 1), rec_idx.(fld2)(end, end)];
        
        plot((t_idx.(fld1) - t_idx_targ) * data.T_STEP * 1e3, 'o-', 'color', colors.(fld1), 'LineWidth', resLW);
        plot((t_idx.(fld2) - t_idx_targ) * data.T_STEP * 1e3, 'o-', 'color', colors.(fld2), 'LineWidth', resLW);
        ylabel('Time interval change from noPert (ms)');
%         set(gca,'FontSize',fs);
%         if bPos==0
%             bar(-.2, chg_iuInt(1), 0.4, 'EdgeColor', colors.(fld1), 'FaceColor', 'none', 'LineWidth', resLW);
%             hold on;
%             bar(0.2, chg_iuInt(2), 0.4, 'EdgeColor', colors.(fld2), 'FaceColor', 'none', 'LineWidth', resLW);
%             y_s = [chg_iuInt,  mean(data.exptData.(mode).chg_IUInt(:,1))+[-1,1]*ste(data.exptData.(mode).chg_IUInt(:,1)),...
%                                mean(data.exptData.(mode).chg_IUInt(:,2))+[-1,1]*ste(data.exptData.(mode).chg_IUInt(:,2))];            
%         else
%             bar(-.3, chg_iuInt(1), 0.2, 'EdgeColor', colors.(fld1), 'FaceColor', 'none', 'LineWidth', resLW);
%             hold on;
%             bar(-.1, sChg_iuInt_p.(mode).(fld1), 0.2, 'EdgeColor', colors.(fld1), 'FaceColor', 'none', 'LineWidth', resLW/2);
%             bar(0.1, chg_iuInt(2), 0.2, 'EdgeColor', colors.(fld2), 'FaceColor', 'none', 'LineWidth', resLW);
%             bar(0.3, sChg_iuInt_p.(mode).(fld2), 0.2, 'EdgeColor', colors.(fld2), 'FaceColor', 'none', 'LineWidth', resLW/2);
%             y_s = [chg_iuInt, sChg_iuInt_p.(mode).(fld1), sChg_iuInt_p.(mode).(fld2), ...
%                    mean(data.exptData.(mode).chg_IUInt(:,1))+[-1,1]*ste(data.exptData.(mode).chg_IUInt(:,1)),...
%                    mean(data.exptData.(mode).chg_IUInt(:,2))+[-1,1]*ste(data.exptData.(mode).chg_IUInt(:,2))];
%         end
%         y_lim=[min(y_s) - 0.1 * range(y_s), max(y_s) + 0.1 * range(y_s)];
%         plot(-.2, mean(data.exptData.(mode).chg_IUInt(:,1)), 'o', 'Color', colors.(fld1), 'LineWidth', 1);
%         plot([-.2,-.2], mean(data.exptData.(mode).chg_IUInt(:,1))+ste(data.exptData.(mode).chg_IUInt(:,1))*[-1,1], '-', 'Color', colors.(fld1), 'LineWidth', 1);
%         plot([-.2+barW/2,-.2-barW/2], mean(data.exptData.(mode).chg_IUInt(:,1))+ste(data.exptData.(mode).chg_IUInt(:,1))*[-1,-1], '-', 'Color', colors.(fld1), 'LineWidth', 1);
%         plot([-.2+barW/2,-.2-barW/2], mean(data.exptData.(mode).chg_IUInt(:,1))+ste(data.exptData.(mode).chg_IUInt(:,1))*[1,1], '-', 'Color', colors.(fld1), 'LineWidth', 1);
%         plot(0.2, mean(data.exptData.(mode).chg_IUInt(:,2)), 'o', 'Color', colors.(fld2), 'LineWidth', 1);
%         plot([0.2,0.2], mean(data.exptData.(mode).chg_IUInt(:,2))+ste(data.exptData.(mode).chg_IUInt(:,2))*[-1,1], '-', 'Color', colors.(fld2), 'LineWidth', 1);
%         plot([0.2+barW/2,0.2-barW/2], mean(data.exptData.(mode).chg_IUInt(:,2))+ste(data.exptData.(mode).chg_IUInt(:,2))*[-1,-1], '-', 'Color', colors.(fld2), 'LineWidth', 1);
%         plot([0.2+barW/2,0.2-barW/2], mean(data.exptData.(mode).chg_IUInt(:,2))+ste(data.exptData.(mode).chg_IUInt(:,2))*[1,1], '-', 'Color', colors.(fld2), 'LineWidth', 1);
%         set(gca, 'XLim', [-0.6,0.6], 'XTick', [-0.2,0.2], 'XTickLabel', {'Down','Up'});
%         set(gca,'YLim',y_lim);
%         ylabel('Change in iuInt (s)');
%         
%         xs=get(gca,'XLim'); ys=get(gca,'YLim');
%         if bPos==0
%             ezlegend([xs(1)+0.05*range(xs),ys(1)+0.725*range(ys),0.775*range(xs),0.25*range(ys)],...
%                 0.3,{['Expt. data (',fld1,', \pm1 SEM)'], ['Simulation (',fld1,')'], ['Expt. data (',fld2,', \pm1 SEM)'], ['Simulation (',fld2,')']},...
%                 {colors.(fld1),colors.(fld1),colors.(fld2),colors.(fld2)},...
%                 [12,12,12,12],{'-o-','-','-o-','-'},{colors.(fld1),colors.(fld1),colors.(fld2),colors.(fld2)},[0.5,resLW,0.5,resLW],[0,0,0,0]);
%             y_s = [chg_iyInt,  mean(data.exptData.(mode).chg_IYInt(:,1))+[-1,1]*ste(data.exptData.(mode).chg_IYInt(:,1)),...
%                                mean(data.exptData.(mode).chg_IYInt(:,2))+[-1,1]*ste(data.exptData.(mode).chg_IYInt(:,2))];
%         else
%             ezlegend([xs(1)+0.05*range(xs),ys(1)+0.725*range(ys),0.775*range(xs),0.25*range(ys)],...
%                 0.3,{['Expt. data (',fld1,', \pm1 SEM)'], ['Simulation sqDIVA (',fld1,')'], ['Simulation old model (',fld1,')']...
%                      ['Expt. data (',fld2,', \pm1 SEM)'], ['Simulation sqDIVA (',fld2,')'], ['Simulation old model (',fld2,')']},...
%                 {colors.(fld1),colors.(fld1),colors.(fld1),colors.(fld2),colors.(fld2),colors.(fld2)},...
%                 [11,11,11,11,11,11],{'-o-','-','-','-o-','-','-'},...
%                 {colors.(fld1),colors.(fld1),colors.(fld1),colors.(fld2),colors.(fld2),colors.(fld2)},...
%                 [0.5,resLW,resLW/2,0.5,resLW,resLW/2],[0,0,0,0,0,0]);
%             y_s = [chg_iyInt, sChg_iyInt_p.(mode).(fld1), sChg_iyInt_p.(mode).(fld2), ...
%                    mean(data.exptData.(mode).chg_IYInt(:,1))+[-1,1]*ste(data.exptData.(mode).chg_IYInt(:,1)),...
%                    mean(data.exptData.(mode).chg_IYInt(:,2))+[-1,1]*ste(data.exptData.(mode).chg_IYInt(:,2))];
%         end
%         y_lim=[min(y_s) - 0.1 * range(y_s), max(y_s) + 0.1 * range(y_s)];
% 
%         subplot('Position',[0.675,0.16,0.3,0.37]);
%         set(gca,'FontSize',fs);
%         if bPos==0
%             bar(-.2, chg_iyInt(1), 0.4, 'EdgeColor', colors.(fld1), 'FaceColor', 'none', 'LineWidth', resLW);
%             hold on;
%             bar(0.2, chg_iyInt(2), 0.4, 'EdgeColor', colors.(fld2), 'FaceColor', 'none', 'LineWidth', resLW);
%         else
%             bar(-.3, chg_iyInt(1), 0.2, 'EdgeColor', colors.(fld1), 'FaceColor', 'none', 'LineWidth', resLW);
%             hold on;
%             bar(-.1, sChg_iyInt_p.(mode).(fld1), 0.2, 'EdgeColor', colors.(fld1), 'FaceColor', 'none', 'LineWidth', resLW/2);
%             bar(0.1, chg_iyInt(2), 0.2, 'EdgeColor', colors.(fld2), 'FaceColor', 'none', 'LineWidth', resLW);
%             bar(0.3, sChg_iyInt_p.(mode).(fld2), 0.2, 'EdgeColor', colors.(fld2), 'FaceColor', 'none', 'LineWidth', resLW/2);
%         end        
%         plot(-.2, mean(data.exptData.(mode).chg_IYInt(:,1)), 'o', 'Color', colors.(fld1), 'LineWidth', 1);
%         plot([-.2,-.2], mean(data.exptData.(mode).chg_IYInt(:,1))+ste(data.exptData.(mode).chg_IYInt(:,1))*[-1,1], '-', 'Color', colors.(fld1), 'LineWidth', 1);
%         plot([-.2+barW/2,-.2-barW/2], mean(data.exptData.(mode).chg_IYInt(:,1))+ste(data.exptData.(mode).chg_IYInt(:,1))*[-1,-1], '-', 'Color', colors.(fld1), 'LineWidth', 1);
%         plot([-.2+barW/2,-.2-barW/2], mean(data.exptData.(mode).chg_IYInt(:,1))+ste(data.exptData.(mode).chg_IYInt(:,1))*[1,1], '-', 'Color', colors.(fld1), 'LineWidth', 1);
%         plot(0.2, mean(data.exptData.(mode).chg_IYInt(:,2)), 'o', 'Color', colors.(fld2), 'LineWidth', 1);
%         plot([0.2,0.2], mean(data.exptData.(mode).chg_IYInt(:,2))+ste(data.exptData.(mode).chg_IYInt(:,2))*[-1,1], '-', 'Color', colors.(fld2), 'LineWidth', 1);
%         plot([0.2+barW/2,0.2-barW/2], mean(data.exptData.(mode).chg_IYInt(:,2))+ste(data.exptData.(mode).chg_IYInt(:,2))*[-1,-1], '-', 'Color', colors.(fld2), 'LineWidth', 1);
%         plot([0.2+barW/2,0.2-barW/2], mean(data.exptData.(mode).chg_IYInt(:,2))+ste(data.exptData.(mode).chg_IYInt(:,2))*[1,1], '-', 'Color', colors.(fld2), 'LineWidth', 1);
%         set(gca, 'XLim', [-0.6,0.6], 'XTick', [-0.2,0.2], 'XTickLabel', {'Down','Up'});
%         set(gca,'YLim',y_lim);
%         ylabel('Change in iyInt (s)');
        
    end
end

%% if bPos, plot the difference in the total curve fitting errors of posFBControl and sqDIVA
if bPos
end 
return