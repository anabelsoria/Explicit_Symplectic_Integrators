function labels3d(units)
    if nargin < 1
        units = 'ND';
    end

    xlabel(sprintf('$x$ [%s]', units), 'interpreter', 'latex')
    ylabel(sprintf('$y$ [%s]', units),'interpreter', 'latex')
    zlabel(sprintf('$z$ [%s]', units),'interpreter', 'latex')
end