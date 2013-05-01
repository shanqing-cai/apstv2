function prop = calc_prop(x, lim, opt)
x = x(~isnan(x));
if isequal(opt, 'gt')
    prop = numel(find(x > lim)) / numel(x);
else
    prop = numel(find(x < lim)) / numel(x);
end
return