function [corrps_wg, uncorrps_wg] = perm_test_tInt_chgs(nPerm, chg_IUInt, chg_IYInt, chg_IU2Int, ...
                                         chg_IY2Int, chg_IU3Int, chg_IY3Int, varargin)
%% Optional input arguments
% --permfile: .mat file name for saving perturbation results
% 
%% Constants
perts = {'accel', 'decel', 'contr'};
nInts = 6;

%% Process optional input arguments
permMatWC = '';
if ~isempty(fsic(varargin, '--permfile'))
    permMatWC = varargin{fsic(varargin, '--permfile') + 1};
    
    assert(length(strfind(permMatWC, '%d')) == 1);
    permMatFN = sprintf(permMatWC, nPerm);    
end


%%
grps = fields(chg_IUInt);
nPerts = numel(perts);

corrps_wg = struct;
uncorrps_wg = struct;
amat = struct;

for i1 = 1 : numel(grps)
    grp = grps{i1};
    
    corrps_wg.(grp) = nan(nInts, nPerts); % Columns: {accel, decel}; Rows: the six time intervals
    uncorrps_wg.(grp) = nan(nInts, nPerts);
    
    for i2 = 1 : numel(perts)
        pert = perts{i2};
        
        if isequal(pert, 'accel') || isequal(pert, 'decel')
            amat.(pert) = [chg_IUInt.(grp)(:, i2), chg_IYInt.(grp)(:, i2), chg_IU2Int.(grp)(:, i2), ...
                           chg_IY2Int.(grp)(:, i2), chg_IU3Int.(grp)(:, i2), chg_IY3Int.(grp)(:, i2)];
        else
            amat.(pert) = [chg_IUInt.(grp)(:, 2) - chg_IUInt.(grp)(:, 1), ...
                           chg_IYInt.(grp)(:, 2) - chg_IYInt.(grp)(:, 1), ...
                           chg_IU2Int.(grp)(:, 2) - chg_IU2Int.(grp)(:, 1), ...
                           chg_IY2Int.(grp)(:, 2) - chg_IY2Int.(grp)(:, 1), ...
                           chg_IU3Int.(grp)(:, 2) - chg_IU3Int.(grp)(:, 1), ...
                           chg_IY3Int.(grp)(:, 2) - chg_IY3Int.(grp)(:, 1)];
        end
                   
        % -- Calculate uncorrected p-values -- %
        for i3 = 1 : nInts
            [~, uncorrps_wg.(grp)(i3, i2)] = ttest(amat.(pert)(:, i3));
        end
        
        if ~isfile(permMatFN)
            minps.(grp).(pert) = nan(1, nPerm);
            
            % -- Permutation -- %
            permMsg = sprintf('Performing within-group permutation test on tInt_chg [%s - %s]', grp, pert);
            nc = print_progress_bar(0, nPerm, permMsg);
            for k1 = 1 : nPerm
                rnd_sign = sign(rand(1, size(amat.(pert), 1)) - 0.5)';
                pmat = amat.(pert) .* repmat(rnd_sign, 1, nInts);

                p_pvals = nan(1, nInts);

                for i3 = 1 : nInts
                    [~, p_pvals(i3)] = ttest(pmat(:, i3));
                end

                minps.(grp).(pert)(k1) = min(p_pvals);
                
                if mod(k1, round(nPerm / 10)) == 0
                    for k2 = 1 : nc; fprintf(1, '\b'); end
                    print_progress_bar(k1, nPerm, permMsg);
                end
            end            
            fprintf(1, '\n');
            
        else
            if ~exist('minps', 'var')
                load(permMatFN);
                assert(exist('minps', 'var') == 1);
                fprintf(1, 'INFO: Loaded existing permutation data from file: %s\n', permMatFN);
            end
        end

        
        fprintf(1, '--- Within-group permutation test: [%s - %s] ---\n', grp, pert);
        for i3 = 1 : nInts
            corrps_wg.(grp)(i3, i2) = numel(find(uncorrps_wg.(grp)(i3, i2) >= minps.(grp).(pert))) / nPerm;
            fprintf(1, '\tInterval #%d: uncorr_p = %f; corr_p = %f\n', ...
                    i3, uncorrps_wg.(grp)(i3, i2), corrps_wg.(grp)(i3, i2));
        end
        fprintf(1, '\n');
    end
end

%% Save permutation data to .mat file (optional)
if ~isempty(permMatFN)
    save(permMatFN, 'minps');
    check_file(permMatFN);
    fprintf(1, 'INFO: %s: nPerm = %d: Saved permutation data to file %s\n', ...
            mfilename, nPerm, permMatFN);
end

return