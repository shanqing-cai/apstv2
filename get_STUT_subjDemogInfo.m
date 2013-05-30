function [ages, genders, SSI4] = get_STUT_subjDemogInfo(dataBookFN, subjIDs, varargin)
[N, T] = xlsread(dataBookFN);

col1 = T(:, 1);
ages = nan(1, numel(subjIDs));
genders = nan(1, numel(subjIDs));
SSI4 = nan(1, numel(subjIDs));

dataBookFN_AS = '';
if ~isempty(fsic(varargin, '--dataBookFN_AS'))
    dataBookFN_AS = varargin{fsic(varargin, '--dataBookFN_AS') + 1};
    check_file(dataBookFN_AS);
    
    [N_AS, T_AS] = xlsread(dataBookFN_AS, 3);
    
    col1_AS = T_AS(:, 1);    
end

for i1 = 1 : numel(subjIDs)
    t_subjID = subjIDs{i1};
    if isequal(t_subjID(end - 1 : end), '_1')
       t_subjID = t_subjID(1 : end - 2);
    end
    
    if isequal(t_subjID(end - 1 : end), '_S')
        t_subjID = t_subjID(1 : end - 2);
    end
    
    if isequal(t_subjID(end - 1 : end), '_T')
        t_subjID = t_subjID(1 : end - 2);
    end
    
    if isequal(t_subjID(1 : 5), 'AS_PS') || ...
            (length(t_subjID) >= 10 && isequal(t_subjID(1 : 9), 'APSTV_PFS'))
        if isempty(dataBookFN_AS)
            error('AS_PS* subjects are found, but --dataBookFN_AS is not supplied');
        end
        
        idx = fsic(col1_AS, t_subjID);
        
        if ~isempty(idx)
            str_DOB = T_AS{idx, 6};
            str_audPertDate = T_AS{idx, 2};
            genders(i1) = isequal(upper(T_AS{idx, 4}(1)), 'M');
        end
    else
        idx = fsic(col1, t_subjID);

        if ~isempty(idx)
            str_DOB = T{idx, 5};
            str_audPertDate = T{idx, 9};
            genders(i1) = isequal(upper(T{idx, 6}(1)), 'M');            
        end
    end
    
    if ~isempty(idx) && ~isempty(str_audPertDate) && ~isempty(str_DOB)
        ages(i1) = (datenum(str_audPertDate) - datenum(str_DOB)) / 365.2442;    
    end

    if ~isempty(strfind(t_subjID, 'PWS_'))
%         if idx < size(N, 1)
            SSI4(i1) = N(idx - 2, 17); % ## WARNING: Ad hoc ##
%         end
    end
    
    
end
return