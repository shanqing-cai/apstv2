function sIDs = get_STUT_EH_subjIDs(grp)
if isequal(grp, 'PWS')
%     sIDs = {'PWS_F01_1', 'PWS_F04_1', 'PWS_F05_1', 'PWS_M01', 'PWS_M02', 'PWS_M04', ...
%             'PWS_M05', 'PWS_M03_1', 'PWS_F06', 'PWS_M08', 'PWS_M07', 'PWS_M06', ...
%             'PWS_M10', 'PWS_M09', 'PWS_M11', 'PWS_F07', 'PWS_M12', 'PWS_M14', 'PWS_M15', ...
%             'PWS_M13'};
    sIDs = {'PWS_F01_1', 'PWS_F04_1', 'PWS_F05_1', 'PWS_M01', 'PWS_M02', 'PWS_M04', ...
            'PWS_M05', 'PWS_M03_1', 'PWS_F06', 'PWS_M08', 'PWS_M07', 'PWS_M06', ...
            'PWS_M10', 'PWS_M09', 'PWS_M11', 'PWS_F07', 'PWS_M12', 'PWS_M14', ...
            'PWS_M15', 'PWS_M13', 'PWS_M16'};
elseif isequal(grp, 'PFS')
    sIDs = {'PFS_F01_1', 'PFS_F05', 'PFS_M01_1', 'PFS_M04', 'PFS_M05', 'PFS_M02', 'PFS_M06', 'PFS_M07', ...
            'PFS_M08', 'PFS_M09', 'PFS_M10', 'PFS_M03', 'PFS_F04', 'PFS_F07', ...
            'PFS_M12', 'PFS_M11', 'PFS_M14', 'PFS_M13'};
end


return