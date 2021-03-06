function metaPlot_2grp(meas, fld1, fld2, measName, y_label, colors, varargin)
fontSize = 12;
if ~isempty(fsic(varargin, 'fontSize'))
    fontSize = varargin{fsic(varargin, 'fontSize') + 1};
end

%%
wndW = 500;
wndH = 360;
if ~isempty(fsic(varargin, 'wndW'))
    wndW = varargin{fsic(varargin, 'wndW') + 1};
end
if ~isempty(fsic(varargin, 'wndH'))
    wndH = varargin{fsic(varargin, 'wndH') + 1};
end

figure('name', measName, 'Position', [200, 200, wndW, wndH]);
set(gca, 'FontSize', fontSize)
if ~isempty(fsic(varargin, 'zeroLine'))
    plot([0, 3], [0, 0], '-', 'Color', [0.5, 0.5, 0.5]);
    hold on;
end

boxplot([meas.(fld1), meas.(fld2)], ...
        [1 * ones(size(meas.(fld1))), 2 * ones(size(meas.(fld2)))]);
hold on;
% boxplot([meas.(fld1), NaN NaN], [1 * ones(size(meas.(fld1))), 2, 2], 'Color', colors.(fld1));
% hold on;
% boxplot([NaN, NaN, meas.(fld1)], [1, 1, 2 * ones(size(meas.(fld1)))], 'Color', colors.(fld2));

if isempty(fsic(varargin, 'noIndPts'))
    plot(repmat(0.8, 1, numel(meas.(fld1))), meas.(fld1), 'o', 'Color', colors.PFS);
    plot(repmat(1.8, 1, numel(meas.(fld2))), meas.(fld2), 'o', 'Color', colors.PWS);
end

XLim = [0, 3];
if ~isempty(fsic(varargin, 'wndH'))
    XLim = varargin{fsic(varargin, 'XLim') + 1};
end

set(gca, 'XLim', XLim, 'XTick', [1, 2], 'XTickLabel', {'PFS', 'PWS'});

plot(1.2, mean(meas.(fld1)), 'o', 'Color', colors.PFS);
plot([1.2, 1.2], mean(meas.(fld1)) + [-1, 1] * ste(meas.(fld1)), '-', 'Color', colors.PFS);
plot(2.2, mean(meas.(fld2)), 'o', 'Color', colors.PWS);
plot([2.2, 2.2], mean(meas.(fld2)) + [-1, 1] * ste(meas.(fld2)), '-', 'Color', colors.PWS);

cohend = abs(cohen_d(meas.(fld2), meas.(fld1)));
xs = get(gca, 'XLim'); ys = get(gca, 'YLim');
text(xs(1) + 0.035 * range(xs), ys(2) - 0.05 * range(ys), ...
    sprintf('PFS: N = %d; PWS: N = %d', numel(meas.(fld1)), numel(meas.(fld2))), 'FontSize', fontSize);
[h_t, p_t] = ttest2(meas.(fld2), meas.(fld1));
text(xs(1) + 0.035 * range(xs), ys(2) - 0.10 * range(ys), sprintf('t-test: p = %.5f', p_t), 'FontSize', fontSize);
text(xs(1) + 0.035 * range(xs), ys(2) - 0.15 * range(ys), sprintf('Cohen''s d = %.5f', cohend), 'FontSize', fontSize);
[p_rs] = ranksum(meas.(fld2), meas.(fld1));
% text(xs(1) + 0.035 * range(xs), ys(2) - 0.15 * range(ys), sprintf('Rank-sum: p = %.5f', p_rs));

ylabel(y_label);

for i0 = 1 : 2
    if i0 == 1
        fld = fld1;
    else
        fld = fld2;
    end
    fprintf('%s: N = %d, mean = %f, SD = %f, SE = %f\n', fld, numel(meas.(fld)), mean(meas.(fld)), std(meas.(fld)), ste(meas.(fld)));
end
return