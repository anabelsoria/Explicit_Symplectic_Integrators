function subplots_same_limits(ax1,ax2)
xlim1 = get(ax1,'Xlim'); xlim2 = get(ax2,'Xlim');
ylim1 = get(ax1,'Ylim'); ylim2 = get(ax2,'Ylim');
xlim_union = [min([xlim1(1), xlim2(1)]), max([xlim1(2), xlim2(2)])];
ylim_union = [min([ylim1(1), ylim2(1)]), max([ylim1(2), ylim2(2)])];

% Apply to both axes (so they match)
xlim(ax1, xlim_union);
ylim(ax1, ylim_union);
xlim(ax2, xlim_union);
ylim(ax2, ylim_union);