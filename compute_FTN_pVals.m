function [p_FTN_down, p_FTN_up, p_FTN_contrast, FDR_thresh_down, FDR_thresh_up, FDR_thresh_contrast] = ...
            compute_FTN_pVals(chgTrajF2_FTN_PFS, chgTrajF2_FTN_PWS, FDR)
flds = fields(chgTrajF2_FTN_PFS);
if ~isempty(fsic(flds, 'down'))
    Fld1 = 'down';
    Fld2 = 'up';
else
    Fld1 = 'accel';
    Fld2 = 'decel';
end
        
nSubjs.PFS = numel(chgTrajF2_FTN_PFS.(Fld1));
nSubjs.PWS = numel(chgTrajF2_FTN_PWS.(Fld2));
N = numel(chgTrajF2_FTN_PFS.(Fld1){1});

p_FTN.(Fld1) = nan(N, 1);
p_FTN.(Fld2) = nan(N, 1);

perts = fields(p_FTN);
perts{end + 1} = 'contrast';
groups = {'PFS', 'PWS'};

chgTrajF2_FTN.PFS = chgTrajF2_FTN_PFS;
chgTrajF2_FTN.PWS = chgTrajF2_FTN_PWS;

for i1 = 1 : numel(perts)
    pert = perts{i1};
    
    t_mat = struct;
    for i2 = 1 : numel(groups)
        grp = groups{i2};
        t_mat.(grp) = nan(N, nSubjs.(grp));
    
        for i3 = 1 : nSubjs.(grp)
            t_mat.(grp)(:, i3) = chgTrajF2_FTN.(grp).(pert){i3};
        end
    end
        
    for i2 = 1 : N
        [h, p_FTN.(pert)(i2)] = ttest2(t_mat.PWS(i2, :), t_mat.PFS(i2, :));
    end
end

p_FTN_down = p_FTN.(Fld1);
p_FTN_up = p_FTN.(Fld2);
p_FTN_contrast = p_FTN.contrast;

FDR_thresh_down = fdr(p_FTN.(Fld1), FDR);
FDR_thresh_up = fdr(p_FTN.(Fld2), FDR);
FDR_thresh_contrast = fdr(p_FTN.contrast, FDR);


return