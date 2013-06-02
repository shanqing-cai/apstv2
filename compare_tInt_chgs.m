function [ps_1, ps_2, ps_12, FDR_p_thresh_1, FDR_p_thresh_2, FDR_p_thresh_12, ...
          varargout] = ...
    compare_tInt_chgs(chg_IUInt, chg_IYInt, chg_IU2Int, ...
                      chg_IY2Int, chg_IU3Int, chg_IY3Int, FDR, varargin)
%% Optional input arguments
% --perm nPerm: Random permutation for nPerm times
% --permfile permMatFN: .mat file name for saving permutation data

% -- DEBUG -- %
% swap = chg_IUInt;
% chg_IUInt = chg_IY3Int;
% chg_IY3Int = swap;

%% Process additional input arguments 
if ~isempty(fsic(varargin, '--perm'))
    nPerm = varargin{fsic(varargin, '--perm') + 1};
    assert(nPerm > 0 && round(nPerm) == nPerm);
    
    if ~isempty(fsic(varargin, '--permfile'))
        permMatWC = varargin{fsic(varargin, '--permfile') + 1};
        assert(length(strfind(permMatWC, '%d')) == 1)
        permMatFN = sprintf(permMatWC, nPerm);
        
        b_permMat_exists = isfile(permMatFN);
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
        mat.accel.(grp) = [chg_IUInt.(grp)(:, 1), ...
                           chg_IYInt.(grp)(:, 1), ...
                           chg_IU2Int.(grp)(:, 1), ...
                           chg_IY2Int.(grp)(:, 1), ...
                           chg_IU3Int.(grp)(:, 1), ...
                           chg_IY3Int.(grp)(:, 1)];
        mat.decel.(grp) = [chg_IUInt.(grp)(:, 2), ...
                           chg_IYInt.(grp)(:, 2), ...
                           chg_IU2Int.(grp)(:, 2), ...
                           chg_IY2Int.(grp)(:, 2), ...
                           chg_IU3Int.(grp)(:, 2), ...
                           chg_IY3Int.(grp)(:, 2)];
        mat.contr.(grp) = [chg_IUInt.(grp)(:, 2) - chg_IUInt.(grp)(:, 1), ...
                           chg_IYInt.(grp)(:, 2) - chg_IYInt.(grp)(:, 1), ...
                           chg_IU2Int.(grp)(:, 2) - chg_IU2Int.(grp)(:, 1), ...
                           chg_IY2Int.(grp)(:, 2) - chg_IY2Int.(grp)(:, 1), ...
                           chg_IU3Int.(grp)(:, 2) - chg_IU3Int.(grp)(:, 1), ...
                           chg_IY3Int.(grp)(:, 2) - chg_IY3Int.(grp)(:, 1)];
    end
    
    pmin_perm = struct;
    corrps = struct;
       
    % --- Between group comparison of contrast --- %
    flds = fields(mat);
    for i1 = 1 : numel(flds)
        fld = flds{i1};
        
        n_PWS = size(mat.(fld).PWS, 1);
        n_PFS = size(mat.(fld).PFS, 1);
        amat.(fld) = cat(1, mat.(fld).PWS, mat.(fld).PFS);

        pmin_perm.(fld) = nan(1, nPerm);

%         if ~isfile(permMatFN)
        if ~b_permMat_exists 
            permMsg = sprintf('Performing permutation test on [%s]', fld);
            nc = print_progress_bar(0, nPerm, permMsg);
            for i1 = 1 : nPerm
                perm_amat = amat.(fld)(randperm(n_PWS + n_PFS), :);

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

                pmin_perm.(fld)(i1) = nanmin(t_ps);

                if mod(i1, round(nPerm / 10)) == 0
                    for k1 = 1 : nc; fprintf(1, '\b'); end 
                    print_progress_bar(i1, nPerm, permMsg);
                end
            end
            fprintf(1, '\n');

            if ~isempty(permMatFN)
                if ~isfile(permMatFN)
                    save(permMatFN, 'pmin_perm');
                else
                    save(permMatFN, 'pmin_perm', '-append');
                end
                check_file(permMatFN);
                fprintf(1, 'INFO: %s: [%s]: Saved permutation data to file %s\n', ...
                        mfilename, fld, permMatFN);
            end
        else
            % -- File already exists, load it -- %
            fprintf(1, 'INFO: %s: [%s]: loading pre-existing permutation data from file %s\n', ...
                    mfilename, fld, permMatFN);
            load(permMatFN);
        end

        fprintf(1, 'Permutation test on [%s]:\n', fld);
        corrps.(fld) = nan(1, numel(ps_12));
        for i1 = 1 : numel(ps_12)
            corrps.(fld)(i1) = numel(find(pmin_perm.(fld) <= ps_12(i1))) / nPerm;
            
            if isequal(fld, 'accel')
                fprintf(1, '\tps_1(%d) = %f; corrp = %f\n', i1, ps_1(i1), corrps.(fld)(i1));
            elseif isequal(fld, 'decel')
                fprintf(1, '\tps_2(%d) = %f; corrp = %f\n', i1, ps_2(i1), corrps.(fld)(i1));
            else
                fprintf(1, '\tps_12(%d) = %f; corrp = %f\n', i1, ps_12(i1), corrps.(fld)(i1));
            end
        end
        fprintf(1, '\n');
    
    
    end
    
    varargout{1} = corrps;
end

return