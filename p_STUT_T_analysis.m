function p_STUT_T_analysis(varargin)
%%
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
dataSet_fn = 'e:\speechres\apstv2\mcode\p_STUT_T_analysis_ds.mat';

bReload = ~isempty(fsic(varargin, 'reload'));
fprintf('bReload = %d.\n', bReload);

if bReload == 1
    [ds.a_uF2_PFS, ds.a_yF2_PFS, ds.normT2_bw_uy_f2_PFS, ds.normT3_bw_yu_f2_PFS, ...
     ds.IUInt_PFS, ds.IYInt_PFS, ds.IU2Int_PFS, ds.IY2Int_PFS, ds.IU3Int_PFS, ds.IY3Int_PFS, ...
     ds.avgPertShiftF2_IOA_PFS, ds.avgChgTrajF2_IOA_PFS, ds.chgTrajF2_FTN_PFS, ds.chgTrajF2_IOA_PFS, ...
     ds.avgPertShiftF2_FTN_PFS, ds.avgChgTrajF2_FTN_PFS, ...
     ds.nTotTrials_PFS, ds.nDscd_all_PFS, ds.nDscd_prodErr_PFS, ...
     ds.subjIDs_PFS] = ...
                p_APSTV2DataAnalysis('apstv2t_expanded', 'output', ...
                                     {'a_uF2', 'a_yF2', 'normT2_bw_uy_f2', 'normT3_bw_yu_f2', ...
                                      'a_IUInt', 'a_IYInt', 'a_IU2Int', 'a_IY2Int', 'a_IU3Int', 'a_IY3Int', ...
                                      'avgPertShiftF2_IOA', 'avgChgTrajF2_IOA', 'chgTrajF2_FTN', 'chgTrajF2_IOA', ...
                                      'avgPertShiftF2_FTN', 'avgChgTrajF2_FTN', ...
                                      'nTotTrials', 'nDscd_all', 'nDscd_prodErr', ...
                                      'subjIDs'});
    ds.tPert_decelAccel_PFS = [];
    for i1 = 1 : numel(ds.subjIDs_PFS)
        if isequal(ds.subjIDs_PFS{i1}(1 : 8), 'APSTV2T_')
            expDir = fullfile('E:\DATA\APSTV2', ds.subjIDs_PFS{i1});
        else
            expDir = fullfile('E:\STUT_DATA\', ds.subjIDs_PFS{i1}(1 : end - 2), 'APSTV2_STUT_T');
        end
        s_pertStats = APSTV2_stats(expDir, 'main', 'rep1', 'main', 'rep20');
        ds.tPert_decelAccel_PFS(i1) = nanmean(s_pertStats.uTimeShift.decel) - nanmean(s_pertStats.uTimeShift.accel);
        ds.F2Pert_dur_PFS(i1) = nanmean([s_pertStats.dur.accel, s_pertStats.dur.decel]);
        ds.F2Pert_uTime_PFS(i1) = nanmean([s_pertStats.uTime.accel, s_pertStats.uTime.decel]);
        close all hidden;
    end
    
%     save(dataSet_fn, 'ds');
%     clear('ds');

    [ds.a_uF2_PWS, ds.a_yF2_PWS, ds.normT2_bw_uy_f2_PWS, ds.normT3_bw_yu_f2_PWS, ...
     ds.IUInt_PWS, ds.IYInt_PWS, ds.IU2Int_PWS, ds.IY2Int_PWS, ds.IU3Int_PWS, ds.IY3Int_PWS, ...
     ds.avgPertShiftF2_IOA_PWS, ds.avgChgTrajF2_IOA_PWS, ds.chgTrajF2_FTN_PWS, ds.chgTrajF2_IOA_PWS, ...
     ds.avgPertShiftF2_FTN_PWS, ds.avgChgTrajF2_FTN_PWS, ...
     ds.nTotTrials_PWS, ds.nDscd_all_PWS, ds.nDscd_prodErr_PWS, ...
     ds.subjIDs_PWS] = ...
                p_APSTV2DataAnalysis('apstv2_stut_t.PWS', 'output', ...
                                     {'a_uF2', 'a_yF2', 'normT2_bw_uy_f2', 'normT3_bw_yu_f2', ...
                                      'a_IUInt', 'a_IYInt', 'a_IU2Int', 'a_IY2Int', 'a_IU3Int', 'a_IY3Int', ...
                                      'avgPertShiftF2_IOA', 'avgChgTrajF2_IOA', 'chgTrajF2_FTN', 'chgTrajF2_IOA', ...
                                      'avgPertShiftF2_FTN', 'avgChgTrajF2_FTN', ...
                                      'nTotTrials', 'nDscd_all', 'nDscd_prodErr', ...
                                      'subjIDs'});
    ds.tPert_decelAccel_PWS = [];
    for i1 = 1 : numel(ds.subjIDs_PWS)
        if isequal(ds.subjIDs_PWS{i1}(1 : 8), 'APSTV2T_')
            expDir = fullfile('E:\DATA\APSTV2', ds.subjIDs_PWS{i1});
        else
            expDir = fullfile('E:\STUT_DATA\', ds.subjIDs_PWS{i1}(1 : end - 2), 'APSTV2_STUT_T');
        end
        s_pertStats = APSTV2_stats(expDir, 'main', 'rep1', 'main', 'rep20');
        ds.tPert_decelAccel_PWS(i1) = nanmean(s_pertStats.uTimeShift.decel) - nanmean(s_pertStats.uTimeShift.accel);
        ds.F2Pert_dur_PWS(i1) = nanmean([s_pertStats.dur.accel, s_pertStats.dur.decel]);
        ds.F2Pert_uTime_PWS(i1) = nanmean([s_pertStats.uTime.accel, s_pertStats.uTime.decel]);
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
[~, age_p] = ttest2(ages.PWS, ages.PFS);
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

%% Results: F2 compensation trajectories: un-normalized (real) time
XLim = [0, 500];
for i1 = 1 : 2
    if i1 == 1; fld = 'accel'; else, fld = 'decel'; end
    metaTracePlot_(ds.avgPertShiftF2_IOA_PFS, ds.avgPertShiftF2_IOA_PWS, fld, XLim, FRAME_DUR * 1e3, colors);
end

IOA_N = round(IOA_WIN / FRAME_DUR);
figure('Position', [50, 100, 800, 400], 'Color', 'w');
for i1 = 1 : 2    
    if i1 == 1; fld = 'accel'; else, fld = 'decel'; end
    subplot('Position', [0.1 + 0.45 * (i1 - 1), 0.15, 0.4, 0.8]);
    set(gca, 'FontSize', fontSize);
    metaTracePlot_(ds.avgChgTrajF2_IOA_PFS, ds.avgChgTrajF2_IOA_PWS, fld, XLim, FRAME_DUR * 1e3, colors);
    set(gca, 'YLim', [-90, 90]);
    
    [p_IOA.accel, p_IOA.decel, FDR_thresh.accel, FDR_thresh.decel] = ...
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

%% Results: F2 compensation trajectories: FTN (full time-normalized)
XLim = [0, 1500];
figure;
for i1 = 1 : 2    
    if i1 == 1; fld = 'accel'; else, fld = 'decel'; end
    metaTracePlot_(ds.avgPertShiftF2_FTN_PFS, ds.avgPertShiftF2_FTN_PWS, fld, XLim, 1, colors);
end

[p_FTN.accel, p_FTN.decel, FDR_thresh.accel, FDR_thresh.decel] = ...
        compute_FTN_pVals(ds.chgTrajF2_FTN_PFS, ds.chgTrajF2_FTN_PWS, 0.05);


for i1 = 1 : 2    
    if i1 == 1; fld = 'accel'; else, fld = 'decel'; end
    figure('Color', 'white', 'Position', [100, 100, 1200, 700], ...
           'Name', ['FTN F2 change trajectorys: ', fld]);
    set(gca, 'FontSize', 20);
    
    metaTracePlot_(ds.avgChgTrajF2_FTN_PFS, ds.avgChgTrajF2_FTN_PWS, fld, XLim, 1, colors);
    taxis0 = linspace(0, 1500, length(p_FTN.(fld)));
    
    set(gca, 'YLim', [-50, 50]);
    xs = get(gca, 'XLim'); ys = get(gca, 'YLim'); 
    plot([0, 0], ys, '-', 'Color', [0.5, 0.5, 0.5]);
    plot([250, 250], ys, '-', 'Color', [0.5, 0.5, 0.5]);
    plot([500, 500], ys, '-', 'Color', [0.5, 0.5, 0.5]);
    plot([750, 750], ys, '-', 'Color', [0.5, 0.5, 0.5]);
    plot([1000, 1000], ys, '-', 'Color', [0.5, 0.5, 0.5]);
    plot([1250, 1250], ys, '-', 'Color', [0.5, 0.5, 0.5]);
    set(gca, 'XLim', [0, 1500], 'XTick', [0 : 250 : 1500], 'XTickLabel', ...
        {'[i]', '[u]_1', '[j]_1', '[u]_2', '[j]_2', '[u]_3', '[j]_3'});
    xlabel('Piecewise normalized time'); ylabel('F2 change (Hz)');
    
    draw_xy_axes;
%     draw_sgn_bar(taxis0, p_FTN.(fld), 0.05, FDR_thresh.(fld), ...
%         ys(1) + 0.20 * range(ys), 0.02 * range(ys), 'k', 0.75 + 0.25 * colors.(fld), colors.(fld));
        
    text(xs(1) + 0.05 * range(xs), ys(2) - 0.05 * range(ys), ...
         [upper(fld(1)), fld(2 : end)], 'FontSize', 20, 'FontWeight', 'Bold');
end

%% Perturbation stats: uTime changes
F2Pert_dur.PFS = ds.F2Pert_dur_PFS * 1e3;
F2Pert_dur.PWS = ds.F2Pert_dur_PWS * 1e3;
metaPlot_2grp(F2Pert_dur, 'PFS', 'PWS', 'Perturbation duration', 'Perturbation duration (ms)', colors, 'zeroLine');

tPert_decelAccel.PFS = ds.tPert_decelAccel_PFS * 1e3;
tPert_decelAccel.PWS = ds.tPert_decelAccel_PWS * 1e3;
metaPlot_2grp(tPert_decelAccel, 'PFS', 'PWS', '[u]_1 timing shifts _in perturbation_', '[u]_1 timing shifts: Decel - Accel (ms)', colors, 'zeroLine');

%% Baseline speaking rates
IUInt_noPert.PFS = ds.IUInt_PFS(:, 1)';
IUInt_noPert.PWS = ds.IUInt_PWS(:, 1)';
metaPlot_2grp(IUInt_noPert, 'PFS', 'PWS', 'IUInt: noPert baseline', 'IYInt: noPert (s)', colors, 'zeroLine');

IYInt_noPert.PFS = ds.IYInt_PFS(:, 1)';
IYInt_noPert.PWS = ds.IYInt_PWS(:, 1)';
metaPlot_2grp(IYInt_noPert, 'PFS', 'PWS', 'IYInt: noPert baseline', 'IYInt: noPert (s)', colors, 'zeroLine');

IY3Int_noPert.PFS = ds.IY3Int_PFS(:, 1)';
IY3Int_noPert.PWS = ds.IY3Int_PWS(:, 1)';
metaPlot_2grp(IY3Int_noPert, 'PFS', 'PWS', 'IY3Int: noPert baseline', 'IY3Int: noPert (s)', colors, 'zeroLine');

%% Relation between speaking rate (local) and IUInt compensation
figure;
plot(ds.IUInt_PFS(:, 1), ds.IUInt_PFS(:, 3) - ds.IUInt_PFS(:, 2), 'o', 'Color', colors.PFS);
hold on;
plot(ds.IUInt_PWS(:, 1), ds.IUInt_PWS(:, 3) - ds.IUInt_PWS(:, 2), 'o', 'Color', colors.PWS);
xlabel('noPert (baseline) IUInt (s)'); ylabel('IUInt contrast: Decel - Accel');

% figure;
% plot(ds.IUInt_PFS(:, 1), ds.I_UYMidF_Int_PFS(:, 3) - ds.I_UYMidF_Int_PFS(:, 2), 'o', 'Color', colors.PFS);
% hold on;
% plot(ds.IUInt_PWS(:, 1), ds.I_UYMidF_Int_PWS(:, 3) - ds.I_UYMidF_Int_PWS(:, 2), 'o', 'Color', colors.PWS);
% xlabel('noPert (baseline) IUInt (s)'); ylabel(strrep('I_UYMidF_Int contrast: Decel - Accel', '_', '\_'));
    
figure;
plot(ds.IYInt_PFS(:, 1), ds.IYInt_PFS(:, 3) - ds.IYInt_PFS(:, 2), 'o', 'Color', colors.PFS);
hold on;
plot(ds.IYInt_PWS(:, 1), ds.IYInt_PWS(:, 3) - ds.IYInt_PWS(:, 2), 'o', 'Color', colors.PWS);
xlabel('noPert (baseline) IYInt (s)'); ylabel(strrep('IYInt contrast: Decel - Accel', '_', '\_'));

%% 
IUInt_decelAccel.PFS = (ds.IUInt_PFS(:, 3) - ds.IUInt_PFS(:, 2))' * 1e3;
IUInt_decelAccel.PWS = (ds.IUInt_PWS(:, 3) - ds.IUInt_PWS(:, 2))' * 1e3;
metaPlot_2grp(IUInt_decelAccel, 'PFS', 'PWS', 'IUInt: decel - accel', ...
    '[u]_1 timing change (Decel - Accel) (ms)', colors, 'zeroLine');

aoctool([IUInt_noPert.PFS, IUInt_noPert.PWS], [IUInt_decelAccel.PFS, IUInt_decelAccel.PWS], ...
        [repmat({'PFS'}, 1, numel(IUInt_decelAccel.PFS)), repmat({'PWS'}, 1, numel(IUInt_decelAccel.PWS))], 0.05, ...
        'Bsaeline [i]-[u]_1 interval', '[i]-[u]_1 interval change: Decel - Accel');
aoctool([tPert_decelAccel.PFS, tPert_decelAccel.PWS], [IUInt_decelAccel.PFS, IUInt_decelAccel.PWS], ...
        [repmat({'PFS'}, 1, numel(IUInt_decelAccel.PFS)), repmat({'PWS'}, 1, numel(IUInt_decelAccel.PWS))], 0.05, ...
        'Bsaeline [i]-[u]_1 interval', '[i]-[u]_1 interval change: Decel - Accel'); 

IUInt_decelNone.PFS = (ds.IUInt_PFS(:, 3) - ds.IUInt_PFS(:, 1))' * 1e3;
IUInt_decelNone.PWS = (ds.IUInt_PWS(:, 3) - ds.IUInt_PWS(:, 1))' * 1e3;
metaPlot_2grp(IUInt_decelNone, 'PFS', 'PWS', 'IUInt: decel - noPert', ...
    'IUInt change: decel - noPert', colors, 'zeroLine');

% aoctool([IUInt_noPert.PFS, IUInt_noPert.PWS], [IUInt_decelNone.PFS, IUInt_decelNone.PWS], ...
%         [repmat({'PFS'}, 1, numel(IUInt_decelNone.PFS)), repmat({'PWS'}, 1, numel(IUInt_decelNone.PWS))], 0.05, ...
%         'Bsaeline [i]-[u]_1 interval', '[i]-[u]_1 interval change: Decel - None');

IYInt_decelNone.PFS = (ds.IYInt_PFS(:, 3) - ds.IYInt_PFS(:, 1))' * 1e3;
IYInt_decelNone.PWS = (ds.IYInt_PWS(:, 3) - ds.IYInt_PWS(:, 1))' * 1e3;
IYInt_decelAccel.PFS = (ds.IYInt_PFS(:, 3) - ds.IYInt_PFS(:, 2))' * 1e3;
IYInt_decelAccel.PWS = (ds.IYInt_PWS(:, 3) - ds.IYInt_PWS(:, 2))' * 1e3;
metaPlot_2grp(IYInt_decelAccel, 'PFS', 'PWS', 'IYInt: decel - accel', ...
    '[j]_1 timing change (Decel - Accel) (ms)', colors, 'zeroLine');

% --- Average change in [u]_1 timing and [j]_1 timing ---
IUIYInt_decelAccel.PFS = mean([IUInt_decelAccel.PFS; IYInt_decelAccel.PFS]);
IUIYInt_decelAccel.PWS = mean([IUInt_decelAccel.PWS; IYInt_decelAccel.PWS]);
metaPlot_2grp(IUIYInt_decelAccel, 'PFS', 'PWS', 'Avg IU-IY Ints: decel - accel', ...
    'Avg. [u]_1 and [j]_1 timing change (ms)', colors, 'zeroLine');
% --- ~ Average change in [u]_1 timing and [j]_1 timing ---

IU2Int_decelNone.PFS = (ds.IU2Int_PFS(:, 3) - ds.IU2Int_PFS(:, 1))' * 1e3;
IU2Int_decelNone.PWS = (ds.IU2Int_PWS(:, 3) - ds.IU2Int_PWS(:, 1))' * 1e3;
IU2Int_decelAccel.PFS = (ds.IU2Int_PFS(:, 3) - ds.IU2Int_PFS(:, 2))' * 1e3;
IU2Int_decelAccel.PWS = (ds.IU2Int_PWS(:, 3) - ds.IU2Int_PWS(:, 2))' * 1e3;
metaPlot_2grp(IU2Int_decelAccel, 'PFS', 'PWS', 'IU2Int: decel - accel', ...
    '[u]_2 timing change (ms)', colors, 'zeroLine');

IU2Int_decelNone.PFS = (ds.IU2Int_PFS(:, 3) - ds.IU2Int_PFS(:, 1))' * 1e3;
IU2Int_decelNone.PWS = (ds.IU2Int_PWS(:, 3) - ds.IU2Int_PWS(:, 1))' * 1e3;
metaPlot_2grp(IU2Int_decelNone, 'PFS', 'PWS', 'IU2Int: decel - noPert', ...
    '[u]_2 timing change (ms)', colors, 'zeroLine');

IY2Int_decelAccel.PFS = (ds.IY2Int_PFS(:, 3) - ds.IY2Int_PFS(:, 2))' * 1e3;
IY2Int_decelAccel.PWS = (ds.IY2Int_PWS(:, 3) - ds.IY2Int_PWS(:, 2))' * 1e3;
metaPlot_2grp(IY2Int_decelAccel, 'PFS', 'PWS', 'IY2Int: decel - accel', ...
    '[j]_2 timing change (ms)', colors, 'zeroLine');

IY2Int_decelNone.PFS = (ds.IY2Int_PFS(:, 3) - ds.IY2Int_PFS(:, 1))' * 1e3;
IY2Int_decelNone.PWS = (ds.IY2Int_PWS(:, 3) - ds.IY2Int_PWS(:, 1))' * 1e3;
metaPlot_2grp(IY2Int_decelNone, 'PFS', 'PWS', 'IY2Int: decel - noPert', ...
    '[j]_2 timing change (ms)', colors, 'zeroLine');

IU3Int_decelNone.PFS = (ds.IU3Int_PFS(:, 3) - ds.IU3Int_PFS(:, 1))' * 1e3;
IU3Int_decelNone.PWS = (ds.IU3Int_PWS(:, 3) - ds.IU3Int_PWS(:, 1))' * 1e3;
IU3Int_decelAccel.PFS = (ds.IU3Int_PFS(:, 3) - ds.IU3Int_PFS(:, 2))' * 1e3;
IU3Int_decelAccel.PWS = (ds.IU3Int_PWS(:, 3) - ds.IU3Int_PWS(:, 2))' * 1e3;
metaPlot_2grp(IU3Int_decelAccel, 'PFS', 'PWS', 'IU2Int: decel - accel', ...
    '[u]_3 timing change (ms)', colors, 'zeroLine');

IY3Int_decelNone.PFS = (ds.IY3Int_PFS(:, 3) - ds.IY3Int_PFS(:, 1))' * 1e3;
IY3Int_decelNone.PWS = (ds.IY3Int_PWS(:, 3) - ds.IY3Int_PWS(:, 1))' * 1e3;
IY3Int_decelAccel.PFS = (ds.IY3Int_PFS(:, 3) - ds.IY3Int_PFS(:, 2))' * 1e3;
IY3Int_decelAccel.PWS = (ds.IY3Int_PWS(:, 3) - ds.IY3Int_PWS(:, 2))' * 1e3;
metaPlot_2grp(IU3Int_decelAccel, 'PFS', 'PWS', 'IY3Int: decel - accel', ...
    '[j]_3 timing change (ms)', colors, 'zeroLine');

aoctool([IYInt_noPert.PFS, IYInt_noPert.PWS], [IYInt_decelAccel.PFS, IYInt_decelAccel.PWS], ...
        [repmat({'PFS'}, 1, numel(IUInt_decelAccel.PFS)), repmat({'PWS'}, 1, numel(IYInt_decelAccel.PWS))], 0.05, ...
        'Bsaeline [i]-[j]_1 interval', '[i]-[j]_1 interval change: Decel - Accel');

IYInt_decelNoPert.PFS = (ds.IYInt_PFS(:, 3) - ds.IYInt_PFS(:, 1))';
IYInt_decelNoPert.PWS = (ds.IYInt_PWS(:, 3) - ds.IYInt_PWS(:, 1))';
metaPlot_2grp(IYInt_decelAccel, 'PFS', 'PWS', 'IYInt: decel - noPert', ...
    'IYInt change: decel - noPert', colors, 'zeroLine');

% I_UYMidF_Int_decelAccel.PFS = (ds.I_UYMidF_Int_PFS(:, 3) - ds.I_UYMidF_Int_PFS(:, 2))';
% I_UYMidF_Int_decelAccel.PWS = (ds.I_UYMidF_Int_PWS(:, 3) - ds.I_UYMidF_Int_PWS(:, 2))';
% metaPlot_2grp(I_UYMidF_Int_decelAccel, 'PFS', 'PWS', 'I_UYMidF_Int: decel - accel', ...
%     'I_UYMidF_Int change: decel - accel', colors, 'zeroLine');

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
groups = {'PFS', 'PWS'};
figure('Position', [50, 75, 1200, 600]);
for i1 = 1 : 2
    subplot('Position', [0.05, 0.15 + (2 - i1) * 0.4, 0.9, 0.4]);
    grp = groups{i1};
    if i1 == 2
        t_int_sum(chg_IUInt.(grp), chg_IYInt.(grp), chg_IU2Int.(grp), chg_IY2Int.(grp), chg_IU3Int.(grp), chg_IY3Int.(grp), colors,...
                'windowWidth', 500, 'windowHeight', 260, 'fontSize', fontSize, 'noFigure', 'YLim', [-5, 15]);
    else
        t_int_sum(chg_IUInt.(grp), chg_IYInt.(grp), chg_IU2Int.(grp), chg_IY2Int.(grp), chg_IU3Int.(grp), chg_IY3Int.(grp), colors,...
                'windowWidth', 500, 'windowHeight', 260, 'fontSize', fontSize, 'noFigure', 'noLabel', 'YLim', [-5, 15]);
    end
    
    if i1 == 2
        ys = get(gca, 'YLim');
        
        text(0.05, ys(2) - 0.045 * range(ys), 'Accel comp.: ', 'color', colors.accel);
        text(0.05, ys(2) - 0.090 * range(ys), 'Decel comp.: ', 'color', colors.decel);
        text(0.05, ys(2) - 0.135 * range(ys), 'Decel-Accel comp.: ', 'color', 'k');
        for k1 = 1 : 6
%             fw1 = ps_1(k1) < FDR_p_thresh_1;
%             fw2 = ps_2(k1) < FDR_p_thresh_2;
%             fw12 = ps_12(k1) < FDR_p_thresh_12;
            fw1 = ps_1(k1) < 0.05;
            fw2 = ps_2(k1) < 0.05;
            fw12 = ps_12(k1) < 0.05;
            if fw1 == 1;    fw1 = 'Bold';   else, fw1 = 'Normal'; end
            if fw2 == 1;    fw2 = 'Bold';   else, fw2 = 'Normal'; end
            if fw12 == 1;    fw12 = 'Bold';   else, fw12 = 'Normal'; end
            text(k1 - 0.1, ys(2) - 0.045 * range(ys), sprintf('%.4f', ps_1(k1)), 'color', colors.accel, 'FontWeight', fw1);
            text(k1 - 0.1, ys(2) - 0.090 * range(ys), sprintf('%.4f', ps_2(k1)), 'color', colors.decel, 'FontWeight', fw2);
            text(k1 - 0.1, ys(2) - 0.135 * range(ys), sprintf('%.4f', ps_12(k1)), 'color', 'k', 'FontWeight', fw12);
        end
    end
end

%% Generate long-format R table for analyses in R
chg_header = ['subject grp sdir tint value'];
chg_header_wide = ['subject grp down_iuint ', ...
                   'down_iyint down_iu2int down_iy2int down_iu3int down_iy3int', ...
                   ' up_iuint up_iyint up_iu2int up_iy2int up_iu3int up_iy3int'];  

chg_long_fn = 'temporalPert_chg.dat';
chg_wide_fn = 'temporalPert_chg_wide.dat';
scp_remote_user = 'cais';
scp_remote_host = 'ba3.mit.edu';
scp_remote_fn = ['/users/cais/STUT/AP_Paper_2/', chg_long_fn];

chg_f = fopen(chg_long_fn, 'w');
chg_f_wide = fopen(chg_wide_fn, 'w');
fprintf(chg_f, '%s\n', chg_header);
fprintf(chg_f_wide, '%s\n', chg_header_wide);
ncols = length(chg_header);
chg_mat = nan(0, ncols);
sCnt = 1;

sDirs = {'accel', 'decel'};
tInts = {'IUInt', 'IYInt', 'IU2Int', 'IY2Int', 'IU3Int', 'IY3Int'};
for i1 = 1 : numel(groups)
    t_grp = groups{i1};
    
    for i2 = 1 : numel(ds.(['subjIDs_', t_grp]))
        dat_line_wide = [];
        
        for i3 = 1 : 2 % sdir
            sDir = sDirs{i3};
            
            for i4 = 1 : 6 % tint
                tInt = tInts{i4};
                if i4 == 1      t_val = chg_IUInt.(t_grp)(i2, i3);
                elseif i4 == 2  t_val = chg_IYInt.(t_grp)(i2, i3);
                elseif i4 == 3  t_val = chg_IU2Int.(t_grp)(i2, i3);
                elseif i4 == 4  t_val = chg_IY2Int.(t_grp)(i2, i3);
                elseif i4 == 5  t_val = chg_IU3Int.(t_grp)(i2, i3);
                elseif i4 == 6  t_val = chg_IY3Int.(t_grp)(i2, i3);
                end
                
                t_line = sprintf('%d %s %s %s %f', ...
                                 sCnt, t_grp, sDir, tInt, t_val);
                fprintf(chg_f, '%s\n', t_line);
                
                dat_line_wide = [dat_line_wide, sprintf('%f ', t_val)];
            end
        end
        
        fprintf(chg_f_wide, '%d %s %s\n', sCnt, t_grp, dat_line_wide);
        
        sCnt = sCnt + 1;
    end
end

fclose(chg_f);
fclose(chg_f_wide);
fprintf('Written long-format table to file %s\n', chg_long_fn);
fprintf('Written wide-format table to file %s\n', chg_wide_fn);

% sftpfrommatlab(scp_remote_user, scp_remote_host, password, ...
%                chg_long_fn, scp_remote_fn);
               

%% Visualization: iuInt and iyInt changes
for i1 = 1 : 2 
    if i1 == 1
        t_meas = chg_IUInt;
    elseif i1 == 2
        t_meas = chg_IYInt;
    end 
    
    figure('Position', [100, 100, 200, 300]);
    barWidth = 1;
    barh(10, 1e3 * mean(t_meas.PFS(:, 1)), barWidth, 'FaceColor', colors.accel, ...
         'LineWidth', 2); 
    hold on;
    plot(1e3 * (mean(t_meas.PFS(:, 1)) + [-1, 1] * ste(t_meas.PFS(:, 1))), ...
         [10, 10], 'k-', 'LineWidth', 2);
    barh(9, 1e3 * mean(t_meas.PFS(:, 2)), barWidth, 'FaceColor', colors.decel, ...
         'LineWidth', 2);
    plot(1e3 * (mean(t_meas.PFS(:, 2)) + [-1, 1] * ste(t_meas.PFS(:, 2))), ...
         [9, 9], 'k-', 'LineWidth', 2);

    barh(7.5, 1e3 * mean(t_meas.PWS(:, 1)), barWidth, 'FaceColor', colors.accel, ...
         'LineWidth', 2); 
    plot(1e3 * (mean(t_meas.PWS(:, 1)) + [-1, 1] * ste(t_meas.PWS(:, 1))), ...
         [7.5, 7.5], 'k-', 'LineWidth', 2);
    barh(6.5, 1e3 * mean(t_meas.PWS(:, 2)), barWidth, 'FaceColor', colors.decel, ...
         'LineWidth', 2);
    plot(1e3 * (mean(t_meas.PWS(:, 2)) + [-1, 1] * ste(t_meas.PWS(:, 2))), ...
         [6.5, 6.5], 'k-', 'LineWidth', 2);

    set(gca, 'XTick', [-4 : 2 : 4]);
    set(gca, 'YLim', [5.5, 11]);
    set(gca, 'YTick', [6.5, 7.5, 9, 10], 'YTickLabel', { 'Decel', 'Accel', 'Decel', 'Accel'});
    if i1 == 1
        xlabel('Change in [i]-[u]_1 interval (ms)');
    else
        xlabel('Change in [i]-[j]_1 interval (ms)');
    end
       
end

%% RM-ANOVA: 1G1W

for i1 = 1 : 2
    if i1 == 1
        grp = 'PFS';
    else
        grp = 'PWS';
    end
    
    decelAccel_contrasts.(grp) = nan(numel(IUInt_decelAccel.(grp)), 0);
    decelNone_contrasts.(grp) = nan(numel(IUInt_decelNone.(grp)), 0);
    for i2 = 1 : 6
        if i2 == 1
            decelAccel_contrasts.(grp) = [decelAccel_contrasts.(grp), IUInt_decelAccel.(grp)'];
            decelNone_contrasts.(grp) = [decelNone_contrasts.(grp), IUInt_decelNone.(grp)'];
        elseif i2 == 2
            decelAccel_contrasts.(grp) = [decelAccel_contrasts.(grp), IYInt_decelAccel.(grp)'];
            decelNone_contrasts.(grp) = [decelNone_contrasts.(grp), IYInt_decelNone.(grp)'];
        elseif i2 == 3
            decelAccel_contrasts.(grp) = [decelAccel_contrasts.(grp), IU2Int_decelAccel.(grp)'];
            decelNone_contrasts.(grp) = [decelNone_contrasts.(grp), IU2Int_decelNone.(grp)'];
        elseif i2 == 4
            decelAccel_contrasts.(grp) = [decelAccel_contrasts.(grp), IY2Int_decelAccel.(grp)'];
            decelNone_contrasts.(grp) = [decelNone_contrasts.(grp), IY2Int_decelNone.(grp)'];
        elseif i2 == 5
            decelAccel_contrasts.(grp) = [decelAccel_contrasts.(grp), IU3Int_decelAccel.(grp)'];
            decelNone_contrasts.(grp) = [decelNone_contrasts.(grp), IU3Int_decelNone.(grp)'];
        elseif i2 == 6
            decelAccel_contrasts.(grp) = [decelAccel_contrasts.(grp), IY3Int_decelAccel.(grp)'];
            decelNone_contrasts.(grp) = [decelNone_contrasts.(grp), IY3Int_decelNone.(grp)'];
        end
    end
    
    decelAccel_contrasts_2dv.(grp) = decelAccel_contrasts.(grp)(:, 1 : 2);
    decelNone_contrasts_2dv.(grp) = decelNone_contrasts.(grp)(:, 1 : 2);
end

do_RMAOV_1G1W(decelAccel_contrasts, 0.05);
do_RMAOV_1G1W(decelNone_contrasts, 0.05);

do_RMAOV_1G1W(decelAccel_contrasts_2dv, 0.05);
do_RMAOV_1G1W(decelNone_contrasts_2dv, 0.05);

%% Correlation between compensation and SSI-4 (severity rating)
SSI4.PWS = [];
for i1 = 1 : numel(ds.subjIDs_PWS)
    sID = strrep(ds.subjIDs_PWS{i1}, '_T', '');
    
    SSI4.PWS(i1) = get_PWS_SSI4(sID);
end

figure;
plot(SSI4.PWS, decelNone_contrasts.PWS(:, 1), 'o', 'Color', colors.PWS);
hold on;
[k, r2, p] = lincorr(SSI4.PWS, decelNone_contrasts.PWS(:, 1));
xs = get(gca, 'XLim'); ys = get(gca, 'YLim');
plot(xs, xs * k(2) + k(1), '--', 'Color', colors.PWS);
set(gca, 'XLim', xs); set(gca, 'YLim', ys);
text(xs(1) + 0.05 * range(xs), ys(2) - 0.075 * range(ys), ...
     sprintf('N = %d; R^2 = %.4f; p = %.4f', numel(find(~isnan(SSI4.PWS))), r2, p));
xlabel('SSI-4');
ylabel('[i]-[u]1 interval compensation (Decel vs. noPert)'); 

% if isequal(grp, 'PWS')
%     SSI4.PWS(i1) = get_PWS_SSI4(strrep(subjIDs.(grp){i1}, '_1', ''));
% end

return


%%
function metaTracePlot_(traj_PFS, traj_PWS, fld, XLim, FRAME_DUR, colors)
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
return