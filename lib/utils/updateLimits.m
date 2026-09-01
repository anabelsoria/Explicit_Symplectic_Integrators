function updateLimits(ax, x, y,x_upper)

    if any(y < 10)
        xlim_curr = get(ax, 'XLim');
        new_xlim = [min(xlim_curr(1), min(x)), min(xlim_curr(2), 1)];
        xlim(ax, new_xlim);

        ylim_curr = get(ax, 'YLim');
        new_ylim = [max(ylim_curr(1), 1E1), max(ylim_curr(2), max(y))];
        ylim(ax, new_ylim);

    else
    % Expand X limits
    xlim_curr = get(ax, 'XLim');
    new_xlim = [min(xlim_curr(1), min(x)), max(xlim_curr(2), max(x))];

    % Apply special rule for X limits
    if exist('x_upper','var')
        if new_xlim(2) > x_upper
            new_xlim(2) = x_upper;
        end
    end
    xlim(ax, new_xlim);


    % Expand Y limits
    ylim_curr = get(ax, 'YLim');
    new_ylim = [min(ylim_curr(1), min(y)), max(ylim_curr(2), max(y))];
    ylim(ax, new_ylim);
    end
end
