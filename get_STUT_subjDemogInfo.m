function [ages, genders, SSI4] = get_STUT_subjDemogInfo(dataBookFN, subjIDs)
[N, T] = xlsread(dataBookFN);

col1 = T(:, 1);
ages = nan(1, numel(subjIDs));
genders = nan(1, numel(subjIDs));
SSI4 = nan(1, numel(subjIDs));

for i1 = 1 : numel(subjIDs)
    t_subjID = subjIDs{i1};
    if isequal(t_subjID(end - 1 : end), '_1')
       t_subjID = t_subjID(1 : end - 2);
    end
    
    idx = fsic(col1, t_subjID); 

    if ~isempty(idx)
        str_DOB = T{idx, 5};
        str_audPertDate = T{idx, 9};
        genders(i1) = isequal(upper(T{idx, 6}(1)), 'M');
    
        if ~isempty(str_audPertDate) && ~isempty(str_DOB)
            ages(i1) = (datenum(str_audPertDate) - datenum(str_DOB)) / 365.2442;    
        end
        
        if ~isempty(strfind(t_subjID, 'PWS_'))
            if idx < size(N, 1)
                SSI4(i1) = N(idx, 15);
            end
        end
        
    end
end
return