function [p_FTN_down, p_FTN_up, p_FTN_contrast, ...
          wg_sigSegs, corrp_wg, ...
          FDR_thresh_down, FDR_thresh_up, FDR_thresh_contrast] = ...
            compute_FTN_pVals(chgTrajF2_FTN_PFS, chgTrajF2_FTN_PWS, FDR, varargin)
%% Optional input arguments
% "--perm" (permutation test): [tpwPThresh, [n0, n1], nPerm]
%       tpwPThresh: time-point-wise p-value threshold (e.g., 0.05);
%       [n0, n1]: the range in which the permutation is to be performed,
%               inclusive interval (e.g., )
%       nPerm: number of permutation iterations (e.g., 10000)
%
% "--permfile": permutation results file (.mat format)
%       Used only with --perm

%% Process optional input arguments
if ~isempty(fsic(varargin, '--perm'))        
    tpwPThresh = varargin{fsic(varargin, '--perm') + 1};
    assert(tpwPThresh > 0 & tpwPThresh <= 1);
    
    perm_ns = varargin{fsic(varargin, '--perm') + 2};
    n0 = perm_ns(1);
    n1 = perm_ns(2);
    assert(n1 > n0);
   
    nPerm = uint32(varargin{fsic(varargin, '--perm') + 3});
    assert(nPerm > 0);
    
    if ~isempty(fsic(varargin, '--permfile'))
        permMatWC = varargin{fsic(varargin, '--permfile') + 1};
        assert(length(strfind(permMatWC, '%s')) == 2);
        assert(length(strfind(permMatWC, '%d')) == 1);
    else
        permMatWC = '';
    end    
else
    nPerm = 0;
end

%%
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
t_FTN.(Fld1) = nan(N, 1);
t_FTN.(Fld2) = nan(N, 1);

perts = fields(p_FTN);
perts{end + 1} = 'contrast';
groups = {'PFS', 'PWS'};

chgTrajF2_FTN.PFS = chgTrajF2_FTN_PFS;
chgTrajF2_FTN.PWS = chgTrajF2_FTN_PWS;


p_wg = struct; % Within-group p-values
corrp_wg = struct; % Within-group corrected p-values
t_wg = struct; % Within-group t-values

wg_sigSegs = struct;

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
        [~, p_FTN.(pert)(i2), ~, t_stats] = ttest2(t_mat.PWS(i2, :), t_mat.PFS(i2, :));
        t_FTN.(pert)(i2) = t_stats.tstat;
    end
    
    % ---- Comparisons (with zero baseline) within groups ---- %
    M = size(t_mat.PFS, 1);
    p_wg.(pert) = struct;
    corrp_wg.(pert) = struct;
    wg_sigSegs.(pert) = struct;
    
    for i2 = 1 : numel(groups)
        grp = groups{i2};
        
        p_wg.(pert).(grp) = nan(1, M);
        t_wg.(pert).(grp) = nan(1, M);
        
        
        for i3 = 1 : M
            [~, t_p, ~, t_stat] = ttest(t_mat.(grp)(i3, :));
            p_wg.(pert).(grp)(i3) = t_p;
            t_wg.(pert).(grp)(i3) = t_stat.tstat;
        end
        
        % --- Optional: random permutation --- %
        % --- NOTE: Permutation for within-group tests are always on the entire
        % duration, i.e., not restricted by n0 and n1 --- %
        
        if nPerm > 0
            [idx_on, idx_off] = get_cont_stretches(p_wg.(pert).(grp) < tpwPThresh);
            
            wgSigSegLens = idx_off - idx_on;
            wg_sigSegs.(pert).(grp) = [idx_on; idx_off];
            
            if ~isempty(wgSigSegLens)
                wgSigSigns = sign(t_wg.(pert).(grp)(idx_on));
                corrp_wg.(pert).(grp) = nan(1, length(idx_on));
            end
            
            permMatFN = sprintf(permMatWC, 'WG', sprintf('%s-%s', pert, grp), nPerm);
            
            if ~isfile(permMatFN)
                wg_maxSigSegLens_perm = zeros(1, nPerm);
                nc = print_progress_bar(0, nPerm, sprintf('Performing permutation test on [WG - %s - %s]', pert, grp));
                for i3 = 1 : nPerm
                    rnd_sign = sign(rand(1, size(t_mat.(grp), 2)) - 0.5);
                    wg_perm_mat = t_mat.(grp) .* repmat(rnd_sign, size(t_mat.(grp), 1), 1);

                    wg_p_perm = nan(1, M);
                    wg_t_perm = nan(1, M);
                    for i4 = 1 : M
                        [~, t_wg_p, ~, t_wg_stats] = ttest(wg_perm_mat(i4, :));
                        wg_p_perm(i4) = t_wg_p;
                        wg_t_perm(i4) = t_wg_stats.tstat;
                    end

                    [idx_on, idx_off] = get_cont_stretches(wg_p_perm < tpwPThresh);

                    if ~isempty(idx_on)
                        wg_sigSegLens_perm = idx_off - idx_on;

                        wg_sigSigns_perm = sign(wg_t_perm(idx_on));

                        wg_sigSegLens_perm = wg_sigSegLens_perm(wg_sigSigns_perm == wgSigSigns(1));

                        if ~isempty(wg_sigSegLens_perm)
                            wg_maxSigSegLens_perm(i3) = max(wg_sigSegLens_perm);
                        end
                    end

                    if mod(i3, round(nPerm / 10)) == 0
                        for k1 = 1 : nc; fprintf(1, '\b'); end 
                        print_progress_bar(i3, nPerm, sprintf('Performing permutation test on [WG - %s - %s]', pert, grp));
                    end
                end
                
                fprintf(1, '\n');
                
                save(permMatFN, 'wg_maxSigSegLens_perm');
                check_file(permMatFN);
                fprintf(1, 'INFO: %s: [WG - %s - %s]: Saved permutation data to file %s\n', ...
                        mfilename, pert, grp, permMatFN);
            else
                % -- File already exists, load it -- %
                fprintf(1, 'INFO: %s: [WG - %s - %s]: loading pre-existing permutation data from file %s\n', ...
                        mfilename, pert, grp, permMatFN);
                load(permMatFN);
            end
        end

        fprintf('Within-group permutation test on [WG - %s - %s]:\n', pert, grp);        
        for i3 = 1 : numel(wgSigSigns)
            corrp_wg.(pert).(grp)(i3) = numel(find(wg_maxSigSegLens_perm > wgSigSegLens(i3))) / double(nPerm);
            fprintf('\tLength = %d; p = %f;\n', wgSigSegLens(i3), corrp_wg.(pert).(grp)(i3));
        end
        fprintf(1, '\n');
    end
    
    % ---- Between-group comparison ---- %
    % --- Random permutation --- %
    if nPerm > 0
        p_FTN_seg = p_FTN.(pert)(n0 : n1);
        t_FTN_seg = t_FTN.(pert)(n0 : n1);
        
        [idx_on, idx_off] = get_cont_stretches(p_FTN_seg < tpwPThresh);
        
        permMatFN = sprintf(permMatWC, 'BCG', pert, nPerm);
        if ~isempty(idx_on)
            segLens = idx_off - idx_on + 1;
            [sigSegLens, idxMaxSeg] = max(segLens);
            sigSign = sign(t_FTN_seg(idx_on(idxMaxSeg)));
            sigSegLens_perm = zeros(1, nPerm);
            
            n_PWS = size(t_mat.PWS, 2);
            n_PFS = size(t_mat.PFS, 2);
            t_mat_tot = [t_mat.PWS, t_mat.PFS];
            t_mat_tot = t_mat_tot(n0 : n1, :);
            
            if ~isfile(permMatFN) % -- File not found, generate from scratch -- %
                nc = print_progress_bar(0, nPerm, sprintf('Performing permutation test on [BGC - %s]', pert));
                for j1 = 1 : nPerm
                    t_mat_perm = t_mat_tot(:, randperm(n_PWS + n_PFS));

                    perm_mat_PWS = t_mat_perm(:, 1 : n_PWS);
                    perm_mat_PFS = t_mat_perm(:, n_PWS + 1 : n_PWS + n_PFS);

                    p_perm = nan(1, n1 - n0 + 1);
                    t_perm = nan(1, n1 - n0 + 1);
                    for j2 = 1 : n1 - n0 + 1;
                        [~, p_perm(j2), ~, t_stat] = ttest2(perm_mat_PWS(j2, :), perm_mat_PFS(j2, :));
                        t_perm(j2) = t_stat.tstat;
                    end

                    [idx_on, idx_off] = get_cont_stretches(p_perm < tpwPThresh);

                    if ~isempty(idx_on)
                        t_sigSegLens = idx_off - idx_on + 1;
                        t_signs = sign(t_perm(idx_on));

                        t_sigSegLens = t_sigSegLens(t_signs == sigSign);
                        if ~isempty(t_sigSegLens)
                            sigSegLens_perm(j1) = max(t_sigSegLens);
                        end
                    end

                    if mod(j1, round(nPerm / 10)) == 0
                        for k1 = 1 : nc; fprintf(1, '\b'); end 
                        print_progress_bar(j1, nPerm, sprintf('Performing permutation test on [BGC - %s]', pert));
                    end

                end
                
                fprintf(1, '\n');
            
                if ~isempty(permMatFN)
                    save(permMatFN, 'sigSegLens_perm');
                    check_file(permMatFN);
                    fprintf(1, 'INFO: %s: [BGC - %s]: Saved permutation data to file %s\n', ...
                            mfilename, pert, permMatFN);
                end

            else
                % -- File already exists, load it -- %
                fprintf(1, 'INFO: %s: [BGC - %s]: loading pre-existing permutation data from file %s\n', ...
                        mfilename, pert, permMatFN);
                load(permMatFN);
            end
            
            ps_perm.(pert) = numel(find(sigSegLens_perm > sigSegLens)) / double(nPerm);
            fprintf('Permutation test on [%s], p = %f\n', pert, ps_perm.(pert));
            
        else
            fprintf(1, 'Pert [BGC - %s]: No significant intervals under tpwPThresh=%f. Moving on.\n', ...
                    pert, tpwPThresh);
        end
        
        
    end
end

p_FTN_down = p_FTN.(Fld1);
p_FTN_up = p_FTN.(Fld2);
p_FTN_contrast = p_FTN.contrast;

FDR_thresh_down = fdr(p_FTN.(Fld1), FDR);
FDR_thresh_up = fdr(p_FTN.(Fld2), FDR);
FDR_thresh_contrast = fdr(p_FTN.contrast, FDR);


return