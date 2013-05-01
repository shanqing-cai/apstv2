function [p_FTN_down, p_FTN_up, FDR_thresh_down, FDR_thresh_up] = ...
            compute_IOA_pVals(chgTrajF2_IOA_PFS, chgTrajF2_IOA_PWS, IOA_N, FDR)
flds = fields(chgTrajF2_IOA_PFS);
if ~isempty(fsic(flds, 'down'))
    Fld1 = 'down';
    Fld2 = 'up';
else
    Fld1 = 'accel';
    Fld2 = 'decel';
end

nSubjs.PFS = numel(chgTrajF2_IOA_PFS.(Fld1));
nSubjs.PWS = numel(chgTrajF2_IOA_PWS.(Fld2));
% IOA_N = numel(chgTrajF2_IOA_PFS.(Fld1){1});

p_FTN.(Fld1) = nan(IOA_N, 1);
p_FTN.(Fld2) = nan(IOA_N, 1);

perts = fields(p_FTN);
groups = {'PFS', 'PWS'};

chgTrajF2_IOA.PFS = chgTrajF2_IOA_PFS;
chgTrajF2_IOA.PWS = chgTrajF2_IOA_PWS;

for i1 = 1 : numel(perts)
    pert = perts{i1};
    
    t_mat = struct;
    for i2 = 1 : numel(groups)
        grp = groups{i2};
        t_mat.(grp) = nan(IOA_N, nSubjs.(grp));
    
        for i3 = 1 : nSubjs.(grp)
            t_mat.(grp)(:, i3) = chgTrajF2_IOA.(grp).(pert){i3}(1 : IOA_N);
        end
    end
        
    for i2 = 1 : IOA_N
        [h, p_FTN.(pert)(i2)] = ttest2(t_mat.PWS(i2, :), t_mat.PFS(i2, :));
    end
end

p_FTN_down = p_FTN.(Fld1);
p_FTN_up = p_FTN.(Fld2);

FDR_thresh_down = fdr(p_FTN.(Fld1), 0.05);
FDR_thresh_up = fdr(p_FTN.(Fld2), 0.05);


return