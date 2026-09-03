function c = plot_colors()
% Named RGB triples from MATLAB's default axes color order, so example
% scripts stop retyping them by hand.

co = get(groot,'defaultAxesColorOrder');

c.blue   = co(1,:);
c.orange = co(2,:);
c.yellow = co(3,:);
c.purple = co(4,:);
c.green  = co(5,:);
c.cyan   = co(6,:);
c.maroon = co(7,:);
end
