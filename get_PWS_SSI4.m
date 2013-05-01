function SSI4 = get_PWS_SSI4(sID, varargin)
% tab = {'PWS_M01', 19; 
%        'PWS_F01', 32;
%        'PWS_F02', NaN;  % Not done.
%        'PWS_F04', 19;
%        'PWS_F05', 16;      
%        'PWS_M02', 43;
%        'PWS_M03', 30;
%        'PWS_M04', 21;
%        'PWS_M05', 33;
%        'PWS_M06', 37;
%        'PWS_F06', 30;
%        'PWS_M07', 14;
%        'PWS_M08', 28;
%        'PWS_M09', 22;
%        'PWS_M10', 19;
%        'PWS_M11', 27};  % OBSOLETE
   
if length(varargin) == 0
    opt = 'total';
else
    opt = varargin{1};
end

cols = {'sID', 'total', 'freq', 'dur', 'concom'};
colIdx = fsic(cols, opt);
if isempty(colIdx)
    error('Unrecognized option: %s', opt)
end
   
% {sID, total, frequency, duration, concomitants}
tab = {'PWS_M01', 23, 8, 8, 6;      %%%%xxxx%%%%
       'PWS_F01', 32, 14, 12, 6;
       'PWS_F02', 16, 8, 6, 2;  % Not done.
       'PWS_F03', 6, 4, 2, 0;
       'PWS_F04', 20, 12, 6, 2;
       'PWS_F05', 17, 8, 6, 3;
       'PWS_M02', 43, 16, 14, 13;
       'PWS_M03', 30, 15, 10, 8;    %%%%xxxx%%%%
       'PWS_M04', 23, 11, 8, 4;
       'PWS_M05', 34, 14, 12, 8;
       'PWS_M06', 37, 15, 14, 8;
       'PWS_F06', 30, 14, 6, 10;
       'PWS_M07', 16, 9, 4, 3;
       'PWS_M08', 28, 9, 10, 9;
       'PWS_M09', 23, 11, 10, 2;
       'PWS_M10', 20, 6, 10, 4;
       'PWS_M11', 27, 14, 8, 5;
       'PWS_F07', 16, 7, 8, 1;
       'PWS_M12', 25, 11, 8, 6;
       'PWS_M13', 29, 9, 8, 12;
       'PWS_M14', 13, 6, 4, 3;
       'PWS_M15', 17, 7, 6, 4;
       'PWS_M16', 36, 17, 8, 11};
   
idx = fsic(tab(:,1), sID);
if ~isempty(idx)
    SSI4 = tab{idx, colIdx};
else
    SSI4 = NaN;
end

return