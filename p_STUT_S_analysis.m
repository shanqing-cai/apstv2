function p_STUT_S_analysis(varargin)
%% CONFIG
colors.PFS = [0, 0, 0];
colors.PWS = [1, 0, 0];

colors.accel=[1,0.5,0];
colors.decel=[0,0.5,0];
colors.up=[1,0,0];
colors.down=[0,0,1];
FRAME_DUR = 16/12e3;

IOA_WIN = 0.5;

fontSize = 12;

if isequal(getHostName, 'CNS-PC34')
    stutDataBookFN = 'C:\Users\scai\Dropbox\STUT\SubjectDataBase-1.xls';
    dataBookFN_AS = 'E:\DATA\APSTV\Subjects-APSTV+SLAP+APAT.xls';
else
    stutDataBookFN = 'C:\Users\systemxp\Documents\My Dropbox\STUT\SubjectDataBase-1.xls';
    dataBookFN_AS = 'E:\DATA\APSTV\Subjects-APSTV+SLAP+APAT.xls';
end


%% Loading data
% dataSet_fn = 'e:\speechres\apstv2\mcode\p_STUT_S_analysis_ds.mat';
dataSet_fn = fullfile(cdds, 'apstv2/mcode/p_STUT_S_analysis_ds.mat');

bReload = ~isempty(fsic(varargin, 'reload'));
fprintf('bReload = %d.\n', bReload);

if bReload == 1
    [ds.a_uF2_PFS, ds.a_yF2_PFS, ds.normT2_bw_uy_f2_PFS, ds.normT3_bw_yu_f2_PFS, ...
     ds.IUInt_PFS, ds.IYInt_PFS, ds.IU2Int_PFS, ds.IY2Int_PFS, ds.IU3Int_PFS, ds.IY3Int_PFS, ...
     ds.avgPertShiftF2_IOA_PFS, ds.avgChgTrajF2_IOA_PFS, ds.chgTrajF2_FTN_PFS, ds.chgTrajF2_IOA_PFS, ...
     ds.avgPertShiftF2_FTN_PFS, ds.avgChgTrajF2_FTN_PFS, ...
     ds.avgChgTrajF2_IOA_contrast_PFS, ds.avgChgTrajF2_FTN_contrast_PFS, ...
     ds.nTotTrials_PFS, ds.nDscd_all_PFS, ds.nDscd_prodErr_PFS, ...
     ds.subjIDs_PFS] = ...
                p_APSTV2DataAnalysis('apstv_expanded', 'output', ...
                                     {'a_uF2', 'a_yF2', 'normT2_bw_uy_f2', 'normT3_bw_yu_f2', ...
                                      'a_IUInt', 'a_IYInt', 'a_IU2Int', 'a_IY2Int', 'a_IU3Int', 'a_IY3Int' ...
                                      'avgPertShiftF2_IOA', 'avgChgTrajF2_IOA', 'chgTrajF2_FTN', 'chgTrajF2_IOA', ...
                                      'avgPertShiftF2_FTN', 'avgChgTrajF2_FTN', ...
                                      'avgChgTrajF2_IOA_contrast', 'avgChgTrajF2_FTN_contrast', ...
                                      'nTotTrials', 'nDscd_all', 'nDscd_prodErr', ...
                                      'subjIDs'});
    ds.F2Pert_upDown_PFS = [];    
    ds.F2Pert_up_PFS = [];
    ds.F2Pert_down_PFS = [];
    for i1 = 1 : numel(ds.subjIDs_PFS)
        if isequal(ds.subjIDs_PFS{i1}(1 : 5), 'AS_PS') || isequal(ds.subjIDs_PFS{i1}(1 : 6), 'APSTV_')
            expDir = fullfile('E:\DATA\APSTV', ds.subjIDs_PFS{i1});
        else
            expDir = fullfile('E:\STUT_DATA\', ds.subjIDs_PFS{i1}(1 : end - 2), 'APSTV2_STUT_S');
        end
        
        s_pertStats = APSTV2_stats(expDir, 'main', 'rep1', 'main', 'rep20');
        ds.F2Pert_upDown_PFS(i1) = nanmean(s_pertStats.f2minChange.up) - nanmean(s_pertStats.f2minChange.down);
        ds.F2Pert_up_PFS(i1) = nanmean(s_pertStats.f2minChange.up);
        ds.F2Pert_down_PFS(i1) = nanmean(s_pertStats.f2minChange.down);
        ds.F2Pert_dur_PFS(i1) = nanmean([s_pertStats.dur.up, s_pertStats.dur.down]);
        ds.F2Pert_uTime_PFS(i1) = nanmean([s_pertStats.uTime.up, s_pertStats.uTime.down]);
        close all hidden;
    end

    [ds.a_uF2_PWS, ds.a_yF2_PWS, ds.normT2_bw_uy_f2_PWS, ds.normT3_bw_yu_f2_PWS, ...
     ds.IUInt_PWS, ds.IYInt_PWS, ds.IU2Int_PWS, ds.IY2Int_PWS, ds.IU3Int_PWS, ds.IY3Int_PWS, ...
     ds.avgPertShiftF2_IOA_PWS, ds.avgChgTrajF2_IOA_PWS, ds.chgTrajF2_FTN_PWS, ds.chgTrajF2_IOA_PWS, ...    
     ds.avgPertShiftF2_FTN_PWS, ds.avgChgTrajF2_FTN_PWS, ...
     ds.avgChgTrajF2_IOA_contrast_PWS, ds.avgChgTrajF2_FTN_contrast_PWS, ...
     ds.nTotTrials_PWS, ds.nDscd_all_PWS, ds.nDscd_prodErr_PWS, ...
     ds.subjIDs_PWS] = ...
                p_APSTV2DataAnalysis('apstv2_stut_s.PWS', 'output', ...
                                     {'a_uF2', 'a_yF2', 'normT2_bw_uy_f2', 'normT3_bw_yu_f2', ...
                                      'a_IUInt', 'a_IYInt', 'a_IU2Int', 'a_IY2Int', 'a_IU3Int', 'a_IY3Int', ...
                                      'avgPertShiftF2_IOA', 'avgChgTrajF2_IOA', 'chgTrajF2_FTN', 'chgTrajF2_IOA', ...
                                      'avgPertShiftF2_FTN', 'avgChgTrajF2_FTN', ...
                                      'avgChgTrajF2_IOA_contrast', 'avgChgTrajF2_FTN_contrast', ...
                                      'nTotTrials', 'nDscd_all', 'nDscd_prodErr', ...
                                      'subjIDs'});
                                  
    ds.F2Pert_upDown_PWS = [];
    ds.F2Pert_up_PWS = [];
    ds.F2Pert_down_PWS = [];
    for i1 = 1 : numel(ds.subjIDs_PWS)
        if isequal(ds.subjIDs_PWS{i1}(1 : 5), 'AS_PS') || isequal(ds.subjIDs_PWS{i1}(1 : 6), 'APSTV_')
            expDir = fullfile('E:\DATA\APSTV', ds.subjIDs_PWS{i1});
        else
            expDir = fullfile('E:\STUT_DATA\', ds.subjIDs_PWS{i1}(1 : end - 2), 'APSTV2_STUT_S');
        end
        
        s_pertStats = APSTV2_stats(expDir, 'main', 'rep1', 'main', 'rep20');
        ds.F2Pert_upDown_PWS(i1) = nanmean(s_pertStats.f2minChange.up) - nanmean(s_pertStats.f2minChange.down);
        ds.F2Pert_up_PWS(i1) = nanmean(s_pertStats.f2minChange.up);
        ds.F2Pert_down_PWS(i1) = nanmean(s_pertStats.f2minChange.down);
        ds.F2Pert_dur_PWS(i1) = nanmean([s_pertStats.dur.up, s_pertStats.dur.down]);
        ds.F2Pert_uTime_PWS(i1) = nanmean([s_pertStats.uTime.up, s_pertStats.uTime.down]);
        close all hidden;
    end
    
    save(dataSet_fn, 'ds');
else
    if isfile(dataSet_fn)
        load(dataSet_fn);
        fprintf('Data set (ds) loaded from %s\n', dataSet_fn);
    else
        error(sprintf('Data set file %s does not exist. Run in "reload" mode first.', dataSet_fn));
    end
end

% F2 Perturbation magnitudes 
% for i1 = 1 : numel(ds.avgPertShiftF2_FTN_PWS.up)
%     
% end

%% Print subject demographic stats
[ages.PFS, genders.PFS] = get_STUT_subjDemogInfo(stutDataBookFN, ds.subjIDs_PFS, ...
                                                 '--dataBookFN_AS', dataBookFN_AS);
[ages.PWS, genders.PWS, SSI4.PWS] = get_STUT_subjDemogInfo(stutDataBookFN, ds.subjIDs_PWS);

fprintf(1, '--- Demographic summary ---\n');
fprintf(1, 'N(PWS) = %d; N(PFS) = %d\n', ...
        numel(ds.subjIDs_PWS), numel(ds.subjIDs_PFS));
fprintf(1, '-- Age: --\n');
fprintf(1, '\tPWS: mean = %.2f, SD = %.2f\n', ...
        mean(ages.PWS), std(ages.PWS));
fprintf(1, '\tPFS: mean = %.2f, SD = %.2f\n', ...
        mean(ages.PFS), std(ages.PFS));
[age_h, age_p] = ttest2(ages.PWS, ages.PFS);
fprintf(1, 'ttest2: p = %f\n', age_p);

fprintf(1, '-- Gender: --\n');
fprintf(1, '\tPWS: %dF%dM\n', ...
        numel(find(genders.PWS == 0)), numel(find(genders.PWS == 1)));
fprintf(1, '\tPFS: %dF%dM\n', ...
        numel(find(genders.PFS == 0)), numel(find(genders.PFS == 1)));
gender_p = chi2test([numel(find(genders.PWS == 0)), numel(find(genders.PWS == 1)); ...
              numel(find(genders.PFS == 0)), numel(find(genders.PFS == 1))]);
fprintf(1, 'ttest2: p = %f\n\n', gender_p);

fprintf(1, '-- SSI-4 scores (PWS) --\n');
fprintf(1, '\tMedian = %.2f; IQR = %.2f; mean = %.2f; SD = %.2f; range = %.2f - %.2f\n\n', ...
        median(SSI4.PWS), iqr(SSI4.PWS), mean(SSI4.PWS), std(SSI4.PWS), ...
        min(SSI4.PWS), max(SSI4.PWS));

    
%% Print stats about discard:
groups = {'PWS', 'PFS'};
for i1 = 1 : numel(groups)
    grp = groups{i1};
    fprintf('%s: %d of %d (%.4f%%) trials discarded,\n', grp, ...
            sum(ds.(['nDscd_all_', (grp)])), sum(ds.(['nTotTrials_', (grp)])), ...
            1e2 * sum(ds.(['nDscd_all_', (grp)])) / sum(ds.(['nTotTrials_', (grp)])));
    fprintf('\t%d (%.4f%%) because of dysfluencies or speech error.\n', ...
            sum(ds.(['nDscd_prodErr_', (grp)])), ...
            1e2 * sum(ds.(['nDscd_prodErr_', (grp)])) / sum(ds.(['nTotTrials_', (grp)])));
end


%% Perturbation stats: F2 perturbation
fprintf(1, '--- Pert stats: Up-Down ---\n');
F2Pert_upDown.PFS = ds.F2Pert_upDown_PFS;
F2Pert_upDown.PWS = ds.F2Pert_upDown_PWS;
metaPlot_2grp(F2Pert_upDown, 'PFS', 'PWS', 'F2 perturbation: Up - Down', 'F2 perturbation: Up - Down (Hz)', colors, 'zeroLine');

fprintf(1, '--- Pert stats: Up ---\n');
F2Pert_up.PFS = ds.F2Pert_up_PFS;
F2Pert_up.PWS = ds.F2Pert_up_PWS;
metaPlot_2grp(F2Pert_up, 'PFS', 'PWS', 'F2 perturbation: Up', 'F2 perturbation: Up (Hz)', colors, 'zeroLine');

fprintf(1, '--- Pert stats: Up ---\n');
F2Pert_down.PFS = ds.F2Pert_down_PFS;
F2Pert_down.PWS = ds.F2Pert_down_PWS;
metaPlot_2grp(F2Pert_down, 'PFS', 'PWS', 'F2 perturbation: Down', 'F2 perturbation: Down (Hz)', colors, 'zeroLine');

%% Results: F2 compensation trajectories: FTN (full time-normalized)
XLim = [0, 1500];
figure('color', 'w');
for i1 = 1 : 2    
    if i1 == 1; fld = 'up'; else; fld = 'down'; end
    metaTracePlot_(ds.avgPertShiftF2_FTN_PFS, ds.avgPertShiftF2_FTN_PWS, fld, XLim, 1, colors);
end

[p_FTN.down, p_FTN.up, p_FTN.contrast, FDR_thresh.down, FDR_thresh.up, FDR_thresh.contrast] = ...
        compute_FTN_pVals(ds.chgTrajF2_FTN_PFS, ds.chgTrajF2_FTN_PWS, 0.05, ...
        '--perm', 0.05, [250, 499], 100, ...
        '--permfile', [mfilename, '_%s_%d.mat']);

for i1 = 1 : 2    
    if i1 == 1; fld = 'up'; else; fld = 'down'; end
    figure('Color', 'white', 'Position', [100, 100, 1200, 400], ...
           'Name', ['FTN F2 change trajectorys: ', fld]);
    set(gca, 'FontSize', 20);
    
    metaTracePlot_(ds.avgChgTrajF2_FTN_PFS, ds.avgChgTrajF2_FTN_PWS, fld, XLim, 1, colors);
    taxis0 = linspace(0, 1500, length(p_FTN.(fld)));
    
    set(gca, 'YLim', [-50, 50]);
    xs = get(gca, 'XLim'); ys = get(gca, 'YLim'); 
    
    label_dx = 25;
    label_dy = 0.06;
    
    plot([0, 0], ys, '-', 'Color', [0.5, 0.5, 0.5]);
    text(0 - label_dx, ys(1) - label_dy * range(ys), '[i]', 'FontSize', 20);
    
    plot([250, 250], ys, '-', 'Color', [0.5, 0.5, 0.5]);
    text(250 - label_dx, ys(1) - label_dy * range(ys), '[u]_1', 'FontSize', 20);
    
    plot([500, 500], ys, '-', 'Color', [0.5, 0.5, 0.5]);
    text(500 - label_dx, ys(1) - label_dy * range(ys), '[j]_1', 'FontSize', 20);
    
    plot([750, 750], ys, '-', 'Color', [0.5, 0.5, 0.5]);
    text(750 - label_dx, ys(1) - label_dy * range(ys), '[u]_2', 'FontSize', 20);
    
    plot([1000, 1000], ys, '-', 'Color', [0.5, 0.5, 0.5]);
    text(1000 - label_dx, ys(1) - label_dy * range(ys), '[j]_2', 'FontSize', 20);
    
    plot([1250, 1250], ys, '-', 'Color', [0.5, 0.5, 0.5]);
    text(1250 - label_dx, ys(1) - label_dy * range(ys), '[u]_3', 'FontSize', 20);
    text(1500 - label_dx, ys(1) - label_dy * range(ys), '[j]_3', 'FontSize', 20);
    
%     set(gca, 'XLim', [0, 1500], 'XTick', [0 : 250 : 1500], 'XTickLabel', ...
%         {'[i]', '[u]_1', '[j]_1', '[u]_2', '[j]_2', '[u]_3', '[j]_3'});
    set(gca, 'XLim', [0, 1500], 'XTick', [0 : 250 : 1500], 'XTickLabel', {});
    xlabel('Piecewise normalized time'); ylabel('F2 change (Hz)');
    
    draw_xy_axes
    draw_sgn_bar(taxis0, p_FTN.(fld), 0.05, FDR_thresh.(fld), ...
        ys(1) + 0.20 * range(ys), 0.02 * range(ys), 'k', [0.5, 0.5, 0.5], 'k');
end

% ds.chgTrajF2_FTN_PFS.up = ds.chgTrajF2_FTN_PFS.contrast;
% ds.chgTrajF2_FTN_PWS.up = ds.chgTrajF2_FTN_PWS.contrast;
% ds.chgTrajF2_FTN_PFS.down = ds.chgTrajF2_FTN_PFS.contrast;
% ds.chgTrajF2_FTN_PWS.down = ds.chgTrajF2_FTN_PWS.contrast;
% [p_FTN_con.down, p_FTN_con.up, FDR_con_thresh.down, FDR_con_thresh.up] = ...
%         compute_FTN_pVals(ds.chgTrajF2_FTN_PFS, ds.chgTrajF2_FTN_PWS, 0.05);

%% Temporal parameters of the Up and Down perturbations 
F2Pert_dur.PFS = ds.F2Pert_dur_PFS;
F2Pert_dur.PWS = ds.F2Pert_dur_PWS;
metaPlot_2grp(F2Pert_dur, 'PFS', 'PWS', 'Perturbation duration', 'Perturbation duration (s)', colors);

F2Pert_uTime.PFS = ds.F2Pert_uTime_PFS;
F2Pert_uTime.PWS = ds.F2Pert_uTime_PWS;
metaPlot_2grp(F2Pert_uTime, 'PFS', 'PWS', 'Perturbation maximum timing', 'Perturbation maximum timing (s)', colors);

%% Baseline speaking rates
IUInt_noPert.PFS = ds.IUInt_PFS(:, 1)';
IUInt_noPert.PWS = ds.IUInt_PWS(:, 1)';
metaPlot_2grp(IUInt_noPert, 'PFS', 'PWS', 'IUInt: noPert baseline', 'IYInt: noPert (s)', colors, 'zeroLine');

IYInt_noPert.PFS = ds.IYInt_PFS(:, 1)';
IYInt_noPert.PWS = ds.IYInt_PWS(:, 1)';
metaPlot_2grp(IYInt_noPert, 'PFS', 'PWS', 'IYInt: noPert baseline', 'IYInt: noPert (s)', colors, 'zeroLine');

IY3Int_noPert.PFS = ds.IY3Int_PFS(:, 1)';
IY3Int_noPert.PWS = ds.IY3Int_PWS(:, 1)';
metaPlot_2grp(IYInt_noPert, 'PFS', 'PWS', 'IY3Int: noPert baseline', 'IY3Int: noPert (s)', colors, 'zeroLine');

%% Relation between speaking rate (local) and IUInt compensation
figure('color', 'w');
plot(ds.IUInt_PFS(:, 1), ds.IUInt_PFS(:, 3) - ds.IUInt_PFS(:, 2), 'o', 'Color', colors.PFS);
hold on;
plot(ds.IUInt_PWS(:, 1), ds.IUInt_PWS(:, 3) - ds.IUInt_PWS(:, 2), 'o', 'Color', colors.PWS);
xlabel('IUInt'); ylabel('IUInt compensation: Up - Down (s)');

figure('color', 'w');
plot(ds.IYInt_PFS(:, 1), ds.IYInt_PFS(:, 3) - ds.IYInt_PFS(:, 2), 'o', 'Color', colors.PFS);
hold on;
plot(ds.IYInt_PWS(:, 1), ds.IYInt_PWS(:, 3) - ds.IYInt_PWS(:, 2), 'o', 'Color', colors.PWS);
xlabel('IYInt'); ylabel('IYInt compensation: Up - Down (s)');

%% Systematic timing change comparison
chg_IUInt.PFS = ds.IUInt_PFS(:, 2:3) - repmat(ds.IUInt_PFS(:, 1), 1, 2);
chg_IYInt.PFS = ds.IYInt_PFS(:, 2:3) - repmat(ds.IYInt_PFS(:, 1), 1, 2);
chg_IU2Int.PFS = ds.IU2Int_PFS(:, 2:3) - repmat(ds.IU2Int_PFS(:, 1), 1, 2);
chg_IY2Int.PFS = ds.IY2Int_PFS(:, 2:3) - repmat(ds.IY2Int_PFS(:, 1), 1, 2);
chg_IU3Int.PFS = ds.IU3Int_PFS(:, 2:3) - repmat(ds.IU3Int_PFS(:, 1), 1, 2);
chg_IY3Int.PFS = ds.IY3Int_PFS(:, 2:3) - repmat(ds.IY3Int_PFS(:, 1), 1, 2);

chg_IUInt.PWS = ds.IUInt_PWS(:, 2:3) - repmat(ds.IUInt_PWS(:, 1), 1, 2);
chg_IYInt.PWS = ds.IYInt_PWS(:, 2:3) - repmat(ds.IYInt_PWS(:, 1), 1, 2);
chg_IU2Int.PWS = ds.IU2Int_PWS(:, 2:3) - repmat(ds.IU2Int_PWS(:, 1), 1, 2);
chg_IY2Int.PWS = ds.IY2Int_PWS(:, 2:3) - repmat(ds.IY2Int_PWS(:, 1), 1, 2);
chg_IU3Int.PWS = ds.IU3Int_PWS(:, 2:3) - repmat(ds.IU3Int_PWS(:, 1), 1, 2);
chg_IY3Int.PWS = ds.IY3Int_PWS(:, 2:3) - repmat(ds.IY3Int_PWS(:, 1), 1, 2);

[ps_1, ps_2, ps_12, FDR_p_thresh_1, FDR_p_thresh_2, FDR_p_thresh_12] ...
    = compare_tInt_chgs(chg_IUInt, chg_IYInt, chg_IU2Int, ...
                        chg_IY2Int, chg_IU3Int, chg_IY3Int, 0.05);
colors.accel = colors.down;
colors.decel = colors.up;
groups = {'PFS', 'PWS'};
figure('Position', [50, 75, 1200, 600], 'color', 'w');
for i1 = 1 : 2
    subplot('Position', [0.1, 0.15 + (2 - i1) * 0.4, 0.8, 0.4]);
    grp = groups{i1};
    if i1 == 2
        t_int_sum(chg_IUInt.(grp), chg_IYInt.(grp), chg_IU2Int.(grp), chg_IY2Int.(grp), chg_IU3Int.(grp), chg_IY3Int.(grp), colors,...
                'windowWidth', 500, 'windowHeight', 260, 'fontSize', fontSize, 'upDown', 'noFigure', 'YLim', [-10, 10]);
    else
        t_int_sum(chg_IUInt.(grp), chg_IYInt.(grp), chg_IU2Int.(grp), chg_IY2Int.(grp), chg_IU3Int.(grp), chg_IY3Int.(grp), colors,...
                'windowWidth', 500, 'windowHeight', 260, 'fontSize', fontSize, 'upDown', 'noFigure', 'noLabel', 'YLim', [-10, 10]);
    end
    
    if i1 == 2
        ys = get(gca, 'YLim');
        
        text(0.05, ys(2) - 0.045 * range(ys), 'Down comp.: ', 'color', colors.down);
        text(0.05, ys(2) - 0.090 * range(ys), 'Up comp.: ', 'color', colors.up);
        text(0.05, ys(2) - 0.135 * range(ys), 'Up-down comp.: ', 'color', 'k');
        for k1 = 1 : 6
            fw1 = ps_1(k1) < FDR_p_thresh_1;
            fw2 = ps_2(k1) < FDR_p_thresh_2;
            fw12 = ps_12(k1) < FDR_p_thresh_12;
            if fw1 == 1;    fw1 = 'Bold';   else, fw1 = 'Normal'; end
            if fw2 == 1;    fw2 = 'Bold';   else, fw2 = 'Normal'; end
            if fw12 == 1;    fw12 = 'Bold';   else, fw12 = 'Normal'; end
            text(k1 - 0.1, ys(2) - 0.045 * range(ys), sprintf('%.4f', ps_1(k1)), 'color', colors.down, 'FontWeight', fw1);
            text(k1 - 0.1, ys(2) - 0.090 * range(ys), sprintf('%.4f', ps_2(k1)), 'color', colors.up, 'FontWeight', fw2);
            text(k1 - 0.1, ys(2) - 0.135 * range(ys), sprintf('%.4f', ps_12(k1)), 'color', 'k', 'FontWeight', fw12);
        end
    end
end

%% Results: F2 compensation trajectories: un-normalized (real) time
figure('color', 'w');
XLim = [0, 500];
for i1 = 1 : 2
    if i1 == 1; fld = 'up'; else; fld = 'down'; end
    metaTracePlot_(ds.avgPertShiftF2_IOA_PFS, ds.avgPertShiftF2_IOA_PWS, fld, XLim, FRAME_DUR * 1e3, colors);
end

IOA_N = round(IOA_WIN / FRAME_DUR);
figure('Position', [50, 100, 800, 400], 'color', 'w');
for i1 = 1 : 2    
    if i1 == 1; fld = 'up'; else; fld = 'down'; end
    subplot('Position', [0.1 + 0.45 * (i1 - 1), 0.15, 0.4, 0.8]);
    set(gca, 'FontSize', fontSize);
    metaTracePlot_(ds.avgChgTrajF2_IOA_PFS, ds.avgChgTrajF2_IOA_PWS, fld, XLim, FRAME_DUR * 1e3, colors);
    set(gca, 'YLim', [-60, 60]);
    
    [p_IOA.down, p_IOA.up, FDR_thresh.down, FDR_thresh.up] = ...
        compute_IOA_pVals(ds.chgTrajF2_IOA_PFS, ds.chgTrajF2_IOA_PWS, IOA_N, 0.05);
    taxis0 = 0 : FRAME_DUR * 1e3 : FRAME_DUR * 1e3 * (IOA_N - 1);
    xs = get(gca, 'XLim'); ys = get(gca, 'YLim'); 
    draw_sgn_bar(taxis0, p_IOA.(fld), 0.05, FDR_thresh.(fld), ...
        ys(1) + 0.10 * range(ys), 0.02 * range(ys), 'k', 0.75 + 0.25 * colors.(fld), colors.(fld));
    
    draw_xy_axes;
    xlabel('Time from [i] (ms)');
    if i1 == 1
        ylabel('F1 change from noPert (Hz) (Mean\pm1 SEM)');
        
        ezlegend([xs(1) + 0.05 * range(xs), ys(2) - 0.3 * range(ys), 0.75 * range(xs), 0.25 * range(ys)], ...
             0.4, {'PWS (mean \pm 1 SEM)', 'PFS (mean\pm 1 SEM)'}, {colors.PWS, colors.PFS}, ...
             [fontSize - 1, fontSize - 1], {'-', '-'}, {colors.PWS, colors.PFS}, ...
             [1, 1], [2, 2]);   
    end
    
end

ds.chgTrajF2_IOA_PFS.up = ds.chgTrajF2_IOA_PFS.contrast;
ds.chgTrajF2_IOA_PWS.up = ds.chgTrajF2_IOA_PWS.contrast;
ds.chgTrajF2_IOA_PFS.down = ds.chgTrajF2_IOA_PFS.contrast;
ds.chgTrajF2_IOA_PWS.down = ds.chgTrajF2_IOA_PWS.contrast;
[p_IOA_con.down, p_IOA_con.up, FDR_con_thresh.down, FDR_con_thresh.up] = ...
        compute_IOA_pVals(ds.chgTrajF2_IOA_PFS, ds.chgTrajF2_IOA_PWS, IOA_N, 0.05);

figure('Position', [50, 100, 400, 400], 'color', 'w');
set(gca, 'FontSize', fontSize);
metaTracePlot_(ds.avgChgTrajF2_IOA_contrast_PFS, ds.avgChgTrajF2_IOA_contrast_PWS, [], XLim, FRAME_DUR * 1e3, colors);
draw_xy_axes;
xlabel('Time from [i] (ms)');
ylabel('Down-Up F2 response contrast (Hz)');
xs = get(gca, 'XLim'); ys = get(gca, 'YLim'); 
draw_sgn_bar(taxis0, p_IOA_con.(fld), 0.05, FDR_con_thresh.(fld), ...
        ys(1) + 0.05 * range(ys), 0.02 * range(ys), 'k', [0.5, 0.5, 0.5], 'k');



%% Between-group FTN comparison
figure('Color', 'white', 'Position', [100, 100, 1600, 700], ...
           'Name', ['FTN F2 change trajectorys: ', fld]);
subplot('Position', [0.1, 0.35, 0.8, 0.55]);
set(gca, 'FontSize', 30);
metaTracePlot_(ds.avgChgTrajF2_FTN_contrast_PFS, ds.avgChgTrajF2_FTN_contrast_PWS, [], XLim, 1, colors);    

ylabel('Down-Up F2 contrast (Hz)');
set(gca, 'YLim', [-25, 55]);
set(gca, 'XTick', []);
draw_xy_axes;

xs = get(gca, 'XLim'); ys = get(gca, 'YLim'); 
draw_sgn_bar(taxis0, p_FTN.contrast, 0.05, FDR_thresh.contrast, ...
        ys(1) + 0.05 * range(ys), 0.02 * range(ys), 'k', [0.5, 0.5, 0.5], 'k');
plot([0, 0], ys, '-', 'Color', [0.5, 0.5, 0.5]);
text(0 - label_dx, ys(1) - label_dy * range(ys), '[i]', 'FontSize', 20);
plot([250 - 2, 250 - 2], ys, '-', 'Color', [0.5, 0.5, 0.5]);
text(250 - label_dx, ys(1) - label_dy * range(ys), '[u]_1', 'FontSize', 20);
plot([500 - 2, 500 - 2], ys, '-', 'Color', [0.5, 0.5, 0.5]);
text(500 - label_dx, ys(1) - label_dy * range(ys), '[j]_1', 'FontSize', 20);
plot([750 - 2, 750 - 2], ys, '-', 'Color', [0.5, 0.5, 0.5]);
text(750 - label_dx, ys(1) - label_dy * range(ys), '[u]_2', 'FontSize', 20);
plot([1000 - 3, 1000 - 3], ys, '-', 'Color', [0.5, 0.5, 0.5]);
text(1000 - label_dx, ys(1) - label_dy * range(ys), '[j]_2', 'FontSize', 20);
plot([1250 - 2, 1250 - 2], ys, '-', 'Color', [0.5, 0.5, 0.5]);
text(1250 - label_dx, ys(1) - label_dy * range(ys), '[u]_3', 'FontSize', 20);
    text(1500 - label_dx, ys(1) - label_dy * range(ys), '[j]_3', 'FontSize', 20);

set(gca, 'XLim', [0, 1500], 'XTick', [0 : 250 : 1500], 'XTickLabel', {});
    
subplot('Position', [0.1, 0.1, 0.8, 0.2]);
set(gca, 'FontSize', 15);
plot(taxis0, p_FTN.contrast);
hold on;
plot([0, 1500], [0.05, 0.05], '--', 'Color', [0.5, 0.5, 0.5]);

plot([0, 0], ys, '-', 'Color', [0.5, 0.5, 0.5]);
plot([250, 250], ys, '-', 'Color', [0.5, 0.5, 0.5]);
plot([500, 500], ys, '-', 'Color', [0.5, 0.5, 0.5]);
plot([750, 750], ys, '-', 'Color', [0.5, 0.5, 0.5]);
plot([1000, 1000], ys, '-', 'Color', [0.5, 0.5, 0.5]);
plot([1250, 1250], ys, '-', 'Color', [0.5, 0.5, 0.5]);
set(gca, 'YLim', [0, 1]);
set(gca, 'XLim', [0, 1500], 'XTick', [0 : 250 : 1500], 'XTickLabel', ...
    {'[i]', '[u]_1', '[j]_1', '[u]_2', '[j]_2', '[u]_3', '[j]_3'});
xlabel('Piecewise normalized time'); ylabel('p-value (uncorrected)');
draw_xy_axes;

% draw_sgn_bar(taxis0, p_FTN.(fld), 0.05, FDR_thresh.(fld), ...
%         ys(1) + 0.20 * range(ys), 0.02 * range(ys), 'k', [0.5, 0.5, 0.5], 'k');

% draw_sgn_bar(taxis0, p_FTN.down, 0.05, FDR_thresh.down, ...
%                  ys(1) + 0.20 * range(ys), 0.02 * range(ys), 'k', 0.75 + 0.25 * colors.down, colors.down);
% draw_sgn_bar(taxis0, p_FTN.up, 0.05, FDR_thresh.up, ...
%                  ys(1) + 0.15 * range(ys), 0.02 * range(ys), 'k', 0.75 + 0.25 * colors.up, colors.up);

%% Results: temporal changes
IUInt_upDown.PFS = (ds.IUInt_PFS(:, 3) - ds.IUInt_PFS(:, 2))';
IUInt_upDown.PWS = (ds.IUInt_PWS(:, 3) - ds.IUInt_PWS(:, 2))';
metaPlot_2grp(IUInt_upDown, 'PFS', 'PWS', 'IUInt: up - down', ...
    'IUInt change: up - down', colors, 'zeroLine');

plot_tIntChg_2grp(ds.IUInt_PFS, ds.IUInt_PWS, colors, ...
                  {'Up', 'Down'}, 'Change in [i]-[u]_1 interval (ms)', [-4, 4]);

aoctool([IUInt_noPert.PFS, IUInt_noPert.PWS], [IUInt_upDown.PFS, IUInt_upDown.PWS], ...
        [repmat({'PFS'}, 1, numel(IUInt_noPert.PFS)), repmat({'PWS'}, 1, numel(IUInt_noPert.PWS))], 0.05, ...
        'Bsaeline [i]-[u]_1 interval', '[i]-[u]_1 interval change: Up - Down');
    
IYInt_upDown.PFS = (ds.IYInt_PFS(:, 3) - ds.IYInt_PFS(:, 2))';
IYInt_upDown.PWS = (ds.IYInt_PWS(:, 3) - ds.IYInt_PWS(:, 2))';
metaPlot_2grp(IYInt_upDown, 'PFS', 'PWS', 'IYInt: up - down', ...
    'IYInt change: up - down', colors, 'zeroLine');

aoctool([IYInt_noPert.PFS, IYInt_noPert.PWS], [IYInt_upDown.PFS, IYInt_upDown.PWS], ...
        [repmat({'PFS'}, 1, numel(IYInt_noPert.PFS)), repmat({'PWS'}, 1, numel(IYInt_noPert.PWS))], 0.05, ...
        'Bsaeline [i]-[j]_1 interval', '[i]-[j]_1 interval change: Up - Down');

plot_tIntChg_2grp(ds.IYInt_PFS, ds.IYInt_PWS, colors, ...
                  {'Up', 'Down'}, 'Change in [i]-[j]_1 interval (ms)', [-4, 4]);
    
%%
F2Pert_uTime.PFS = ds.F2Pert_uTime_PFS;
F2Pert_uTime.PWS = ds.F2Pert_uTime_PWS;
aoctool([F2Pert_uTime.PFS, F2Pert_uTime.PWS], [IUInt_upDown.PFS, IUInt_upDown.PWS], ...
        [repmat({'PFS'}, 1, numel(IUInt_noPert.PFS)), repmat({'PWS'}, 1, numel(IUInt_noPert.PWS))], 0.05, ...
        'Perturbation maximum lag (s)', '[i]-[u]_1 interval change: Up - Down');
aoctool([F2Pert_uTime.PFS, F2Pert_uTime.PWS], [IYInt_upDown.PFS, IYInt_upDown.PWS], ...
        [repmat({'PFS'}, 1, numel(IUInt_noPert.PFS)), repmat({'PWS'}, 1, numel(IUInt_noPert.PWS))], 0.05, ...
        'Perturbation maximum lag (s)', '[i]-[j]_1 interval change: Up - Down');

%% Results: spatial (F2 magnitude) changes
uF2_downUp.PFS = (ds.a_uF2_PFS(:, 2) - ds.a_uF2_PFS(:, 3))';
uF2_downUp.PWS = (ds.a_uF2_PWS(:, 2) - ds.a_uF2_PWS(:, 3))';
metaPlot_2grp(uF2_downUp, 'PFS', 'PWS', 'uF2: down - up', ...
    'uF2 change: down - up (Hz)', colors, 'zeroLine');



yF2_downUp.PFS = (ds.a_yF2_PFS(:, 2) - ds.a_yF2_PFS(:, 3))';
yF2_downUp.PWS = (ds.a_yF2_PWS(:, 2) - ds.a_yF2_PWS(:, 3))';
metaPlot_2grp(yF2_downUp, 'PFS', 'PWS', 'yF2: down - up', ...
    'yF2 change: down - up (Hz)', colors, 'zeroLine');



uyMidF2_downUp.PFS = (ds.normT2_bw_uy_f2_PFS{2}(:, 2) - ds.normT2_bw_uy_f2_PFS{2}(:, 3))';
uyMidF2_downUp.PWS = (ds.normT2_bw_uy_f2_PWS{2}(:, 2) - ds.normT2_bw_uy_f2_PWS{2}(:, 3))';
metaPlot_2grp(uyMidF2_downUp, 'PFS', 'PWS', 'uyMidF2: down - up', ...
    'uyMidF2: down - up (Hz)', colors, 'zeroLine');

uyMidF2_downUp.PFS = (ds.normT2_bw_uy_f2_PFS{2}(:, 2) - ds.normT2_bw_uy_f2_PFS{2}(:, 3))';
uyMidF2_downUp.PWS = (ds.normT2_bw_uy_f2_PWS{2}(:, 2) - ds.normT2_bw_uy_f2_PWS{2}(:, 3))';
metaPlot_2grp(uyMidF2_downUp, 'PFS', 'PWS', 'uyMidF2: down - up', ...
    'uyMidF2: down - up (Hz)', colors, 'zeroLine');

yuMidF2_downUp.PFS = (ds.normT3_bw_yu_f2_PFS{2}(:, 2) - ds.normT3_bw_yu_f2_PFS{2}(:, 3))';
yuMidF2_downUp.PWS = (ds.normT3_bw_yu_f2_PWS{2}(:, 2) - ds.normT3_bw_yu_f2_PWS{2}(:, 3))';
metaPlot_2grp(yuMidF2_downUp, 'PFS', 'PWS', 'yuMidF2: down - up', ...
    'yuMidF2: down - up (Hz)', colors, 'zeroLine');

u2F2_downUp.PFS = (ds.normT3_bw_yu_f2_PFS{end}(:, 2) - ds.normT3_bw_yu_f2_PFS{end}(:, 3))';
u2F2_downUp.PWS = (ds.normT3_bw_yu_f2_PWS{end}(:, 2) - ds.normT3_bw_yu_f2_PWS{end}(:, 3))';
metaPlot_2grp(u2F2_downUp, 'PFS', 'PWS', 'u2F2: down - up', ...
    'u2F2: down - up (Hz)', colors, 'zeroLine');

return

%%
function metaTracePlot_(traj_PFS, traj_PWS, fld, XLim, FRAME_DUR, colors)
if ~isempty(fld)
    plot(XLim, [0, 0], '-', 'Color', [0.5, 0.5, 0.5]); hold on;
    tAxis = 0 : FRAME_DUR : FRAME_DUR * (size(traj_PFS.(fld), 1) - 1);
    plot(tAxis, traj_PFS.(fld)(:, 1), 'LineWidth', 1, 'Color', colors.PFS); hold on;
    plot_sd_t(tAxis, traj_PFS.(fld)(:, 1), traj_PFS.(fld)(:, 2) ./ sqrt(traj_PFS.(fld)(:, 3)), 0.5 + 0.5 * colors.PFS, 'patch');

    % plot(tAxis, traj_PFS.(fld)(:, 1) - traj_PFS.(fld)(:, 2) ./ sqrt(traj_PFS.(fld)(:, 3)), ...
    %     'LineWidth', 0.5, 'Color', colors.(fld));
    % plot(tAxis, traj_PFS.(fld)(:, 1) + traj_PFS.(fld)(:, 2) ./ sqrt(traj_PFS.(fld)(:, 3)), ...
    %     'LineWidth', 0.5, 'Color', colors.(fld));

    tAxis = 0 : FRAME_DUR : FRAME_DUR * (size(traj_PWS.(fld), 1) - 1);
    plot(tAxis, traj_PWS.(fld)(:, 1), '-', 'LineWidth', 1, 'Color', colors.PWS); hold on;
    plot_sd_t(tAxis, traj_PWS.(fld)(:, 1), traj_PWS.(fld)(:, 2) ./ sqrt(traj_PWS.(fld)(:, 3)), 0.5 + 0.5 * colors.PWS, 'patch');

    % plot(tAxis, traj_PWS.(fld)(:, 1) - traj_PWS.(fld)(:, 2) ./ sqrt(traj_PWS.(fld)(:, 3)), ...
    %     '--', 'LineWidth', 0.5, 'Color', colors.(fld));
    % plot(tAxis, traj_PWS.(fld)(:, 1) + traj_PWS.(fld)(:, 2) ./ sqrt(traj_PWS.(fld)(:, 3)), ...
    %     '--', 'LineWidth', 0.5, 'Color', colors.(fld));
    % tAxis = 0 : FRAME_DUR : FRAME_DUR * (size(traj_PWS.down, 1) - 1);

    set(gca, 'XLim', XLim);
else
    plot(XLim, [0, 0], '-', 'Color', [0.5, 0.5, 0.5]); hold on;
    tAxis = 0 : FRAME_DUR : FRAME_DUR * (size(traj_PFS, 1) - 1);
    plot(tAxis, traj_PFS(:, 1), 'LineWidth', 1, 'Color', colors.PFS); hold on;
    plot_sd_t(tAxis, traj_PFS(:, 1), traj_PFS(:, 2) ./ sqrt(traj_PFS(:, 3)), 0.5 + 0.5 * colors.PFS, 'patch');

    tAxis = 0 : FRAME_DUR : FRAME_DUR * (size(traj_PWS, 1) - 1);
    plot(tAxis, traj_PWS(:, 1), '-', 'LineWidth', 1, 'Color', colors.PWS); hold on;
    plot_sd_t(tAxis, traj_PWS(:, 1), traj_PWS(:, 2) ./ sqrt(traj_PWS(:, 3)), 0.5 + 0.5 * colors.PWS, 'patch');

    set(gca, 'XLim', XLim);
end
 return