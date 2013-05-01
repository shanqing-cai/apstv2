function [ps_1, ps_2, ps_12, FDR_p_thresh_1, FDR_p_thresh_2, FDR_p_thresh_12] = ...
    compare_tInt_chgs(chg_IUInt, chg_IYInt, chg_IU2Int, ...
                      chg_IY2Int, chg_IU3Int, chg_IY3Int, FDR)
ps_1 = nan(1, 6);
ps_2 = nan(1, 6);

for i1 = 1 : 2    
    t_ps = nan(1, 6);
    
    [h, t_ps(1)] = ttest2(chg_IUInt.PFS(:, i1), chg_IUInt.PWS(:, i1));
    [h, t_ps(2)] = ttest2(chg_IYInt.PFS(:, i1), chg_IYInt.PWS(:, i1));
    [h, t_ps(3)] = ttest2(chg_IU2Int.PFS(:, i1), chg_IU2Int.PWS(:, i1));
    [h, t_ps(4)] = ttest2(chg_IY2Int.PFS(:, i1), chg_IY2Int.PWS(:, i1));
    [h, t_ps(5)] = ttest2(chg_IU3Int.PFS(:, i1), chg_IU3Int.PWS(:, i1));
    [h, t_ps(6)] = ttest2(chg_IY3Int.PFS(:, i1), chg_IY3Int.PWS(:, i1));

    if i1 == 1
        ps_1 = t_ps;
    else
        ps_2 = t_ps;
    end
end

ps_12 = nan(1, 6);
[h, ps_12(1)] = ttest2(chg_IUInt.PFS(:, 2) - chg_IUInt.PFS(:, 1), chg_IUInt.PWS(:, 2) - chg_IUInt.PWS(:, 1));
[h, ps_12(2)] = ttest2(chg_IYInt.PFS(:, 2) - chg_IYInt.PFS(:, 1), chg_IYInt.PWS(:, 2) - chg_IYInt.PWS(:, 1));
[h, ps_12(3)] = ttest2(chg_IU2Int.PFS(:, 2) - chg_IU2Int.PFS(:, 1), chg_IU2Int.PWS(:, 2) - chg_IU2Int.PWS(:, 1));
[h, ps_12(4)] = ttest2(chg_IY2Int.PFS(:, 2) - chg_IY2Int.PFS(:, 1), chg_IY2Int.PWS(:, 2) - chg_IY2Int.PWS(:, 1));
[h, ps_12(5)] = ttest2(chg_IU3Int.PFS(:, 2) - chg_IU3Int.PFS(:, 1), chg_IU3Int.PWS(:, 2) - chg_IU3Int.PWS(:, 1));
[h, ps_12(6)] = ttest2(chg_IY3Int.PFS(:, 2) - chg_IY3Int.PFS(:, 1), chg_IY3Int.PWS(:, 2) - chg_IY3Int.PWS(:, 1));

FDR_p_thresh_1 = fdr(ps_1, FDR);
FDR_p_thresh_2 = fdr(ps_2, FDR);
FDR_p_thresh_12 = fdr(ps_12, FDR);

return