function [ps_1, ps_2, ps_12, FDR_p_thresh_1, FDR_p_thresh_2, FDR_p_thresh_12, ...
          varargout] = ...
    compare_tInt_chgs(chg_IUInt, chg_IYInt, chg_IU2Int, ...
                      chg_IY2Int, chg_IU3Int, chg_IY3Int, FDR, varargin)
%% Optional input arguments
% --perm nPerm: Random permutation for nPerm times
% --permfile permMatFN: .mat file name for saving permutation data

%% Process additional input arguments 
if ~isempty(fsic(varargin, '--perm'))
    nPerm = varargin{fsic(varargin, '--perm') + 1};
    assert(nPerm > 0 && round(nPerm) == nPerm);
    
    if ~isempty(fsic(varargin, '--permfile'))
        permMatWC = varargin{fsic(varargin, '--permfile') + 1};
        assert(length(strfind(permMatWC, '%d')) == 1)
        permMatFN = sprintf(permMatWC, nPerm);
    else
        permMatFN = '';
    end
    
else
    nPerm = 0;
    
end    

%%
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

%% Perform random permutation test (optional)
groups = fields(chg_IUInt);
if nPerm > 0
    for i1 = 1 : numel(groups)
        grp = groups{i1};
        mat.contr.(grp) = [chg_IUInt.(grp)(:, 2) - chg_IUInt.(grp)(:, 1), ...
                           chg_IYInt.(grp)(:, 2) - chg_IYInt.(grp)(:, 1), ...
                           chg_IU2Int.(grp)(:, 2) - chg_IU2Int.(grp)(:, 1), ...
                           chg_IY2Int.(grp)(:, 2) - chg_IY2Int.(grp)(:, 1), ...
                           chg_IU3Int.(grp)(:, 2) - chg_IU3Int.(grp)(:, 1), ...
                           chg_IY3Int.(grp)(:, 2) - chg_IY3Int.(grp)(:, 1)];
    end
    
    n_PWS = size(mat.contr.PWS, 1);
    n_PFS = size(mat.contr.PFS, 1);
    amat.contr = cat(1, mat.contr.PWS, mat.contr.PFS);
    
    pmin_perm = nan(1, nPerm);
    
    if ~isfile(permMatFN)
        permMsg = 'Performing permutation test on [contr]';
        nc = print_progress_bar(0, nPerm, permMsg);
        for i1 = 1 : nPerm
            perm_amat = amat.contr(randperm(n_PWS + n_PFS), :);

            perm_mat_PWS = perm_amat(1 : n_PWS, :);
            perm_mat_PFS = perm_amat(n_PWS + 1 : n_PWS + n_PFS, :);

            t_ps = nan(1, size(perm_amat, 2));       
            for i2 = 1 : size(perm_amat, 2)
                [~, t_p, ~, t_stats] = ttest2(perm_mat_PWS(:, i2), perm_mat_PFS(:, i2));
    %             if t_stats.tstat < 0
    %                 t_ps(i2) = t_p;
    %             end
                t_ps(i2) = t_p;
            end

            pmin_perm(i1) = nanmin(t_ps);

            if mod(i1, round(nPerm / 10)) == 0
                for k1 = 1 : nc; fprintf(1, '\b'); end 
                print_progress_bar(i1, nPerm, permMsg);
            end
        end
        fprintf(1, '\n');
        
        if ~isempty(permMatFN)
            save(permMatFN, 'pmin_perm');
            check_file(permMatFN);
            fprintf(1, 'INFO: %s: [contr]: Saved permutation data to file %s\n', ...
                    mfilename, permMatFN);
        end
    else
        % -- File already exists, load it -- %
        fprintf(1, 'INFO: %s: [contr]: loading pre-existing permutation data from file %s\n', ...
                mfilename, permMatFN);
        load(permMatFN);
    end
    
    fprintf(1, 'Permutation test on [contr]:\n');
    corrps.contr = nan(1, numel(ps_12));
    for i1 = 1 : numel(ps_12)
        corrps.contr(i1) = numel(find(pmin_perm <= ps_12(i1))) / nPerm;
        fprintf(1, '\tps_12(%d) = %f; corrp = %f\n', i1, ps_12(i1), corrps.contr(i1))
    end
    fprintf(1, '\n');
    
    varargout{1} = corrps;        
end

return