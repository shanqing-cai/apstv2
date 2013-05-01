function draw_sgn_bar(tAxis, ps, pthresh1, pthresh2, y, width, color_edge, color1, color2)

rectangle('Position', [tAxis(1), y - 0.5 * width, tAxis(end) - tAxis(1) + tAxis(2) - tAxis(1), width], ...
                  'EdgeColor', color_edge, 'FaceColor', 'none');
for i0 = 1 : 2
    if i0 == 1
        bk1 = ps < pthresh1;
    else
        bk1 = ps < pthresh2;
    end
    
    bk1(1) = 0;
    if bk1(end) == 1
        bk1 = [bk1; 0];
    end
    [iBegin, iEnd] = getNonZeroChunks(bk1);
    for i1 = 1 : numel(iBegin)
        if i0 == 1
            fillColor = color1;
        else
            fillColor = color2;
        end
        
        if tAxis(iEnd(i1)) - tAxis(iBegin(i1)) > 0
            rectangle('Position', [tAxis(iBegin(i1)), y - 0.5 * width, tAxis(iEnd(i1)) - tAxis(iBegin(i1)), width], ...
                      'EdgeColor', color_edge, 'FaceColor', fillColor);
        else
            rectangle('Position', [tAxis(iBegin(i1) - 1), y - 0.5 * width, tAxis(2) - tAxis(1), width], ...
                      'EdgeColor', color_edge, 'FaceColor', fillColor);
%             plot([tAxis(iBegin(i1)), tAxis(iBegin(i1))], [y - 0.5 * width, y + 0.5 * width], ...
%                       'Color', color_edge);
        end
        hold on;
    end
end

set(gca, 'XLim', [tAxis(1), tAxis(end)]);

return
