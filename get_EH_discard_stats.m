function [nDscd_prodErr, nTotTrials, nDscd_all] = get_EH_discard_stats(pdata)
flds = fields(pdata.utters);

nTotTrials = 0;
nDscd_prodErr = 0;
nDscd_all = 0;

for i1 = 1 : numel(flds)
    fld = flds{i1};
    nTrials = numel(pdata.utters.(fld));
    for i2 = 1 : nTrials
        t_utter = pdata.utters.(fld){i2};
        if ~isfield(t_utter, 'bDiscard')
            continue;
        end
        nTotTrials = nTotTrials + 1;
        if t_utter.bDiscard == 1
            nDscd_all = nDscd_all + 1;
            if isfield(t_utter, 'ratingComments')
                if ~isempty(strfind(lower(t_utter.ratingComments), 'sfluent')) ...
                || ~isempty(strfind(lower(t_utter.ratingComments), 'eech err')) ...
                || ~isempty(strfind(lower(t_utter.ratingComments), 'production error'))
                    nDscd_prodErr = nDscd_prodErr + 1;
                end
            end
        end
    end
end
return