function res = do_RMAOV_1G1W(meas, alpha)
% meas: two groups: PWS and PFS
%   meas.(grp): rows: different subjects
%               columns: different measures

tab = nan(0, 4); 

groups = {'PFS','PWS'};
subjCnt = 1;
for i1 = 1 : numel(groups)
    grp = groups{i1};
    for i2 = 1 : size(meas.(grp), 1)
        for i3 = 1 : size(meas.(grp), 2)
            tab_row = [meas.(grp)(i2, i3), i1, i3, subjCnt];
            tab = [tab; tab_row];
        end
        subjCnt = subjCnt + 1;
    end
end

BWAOV2(tab, alpha);
return