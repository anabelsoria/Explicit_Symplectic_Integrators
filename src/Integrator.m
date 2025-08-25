    classdef Integrator < handle
    
        properties
            name % Name of the integrator
            prob % Problem class containing parameters, initial conditions, etc.
            sol  % Propagated solution (trajectory and tspan)
            time_solver = NaN % Time taken by integrator
        end
    
        methods

            function obj = Integrator(prob,name, sol)
                arguments
                    prob = []
                    name = 'NaN'
                    sol  = []
                end
            % Constructor for the RK class
            obj.name = name;
            obj.prob = prob;
            obj.sol  = sol;
        end
    
            function plot_traj(obj,options)
                arguments
                    obj
                    options.title = []
                    options.fig = []
                    options.plot_2d_xy = false
                    options.plot_2d_xz = false
                    options.font_size = 12;
                    options.show_steps = false
                    options.color = []
                end
    
                if isempty(options.fig)
                    figure;
                    if isfield(obj.prob.DS,'r2')
                    scatter3(obj.prob.DS.r2,0,0,15,'k','filled',DisplayName='Moon')
                    end
                else
                    figure(options.fig);
                end
    
                hold on;
                sgtitle(options.title, 'Interpreter', 'latex');

                if options.show_steps
                    name_plot = obj.name + " $(n = " + num2str(obj.sol.Nsteps) + ")$";
                else
                    name_plot = obj.name;
                end

                if isempty(options.color)
                    plot3(obj.sol.x(1,:),obj.sol.x(2,:),obj.sol.x(3,:), ...
                        LineWidth=1,DisplayName=name_plot)

                else
                    plot3(obj.sol.x(1,:),obj.sol.x(2,:),obj.sol.x(3,:), ...
                        LineWidth=1,DisplayName=name_plot,Color=options.color)

                end

                if options.plot_2d_xy
                    view(2)
                elseif options.plot_2d_xz
                    view(2)
                else
                    view(3)
                end
    
                labels3d('LU'); %grid on;
                legend('Interpreter', 'latex', 'Location','best');
                set(findall(gcf, '-property', 'FontSize'), 'FontSize', options.font_size);
    
            end
    
            function plot_conserved_quantity(obj,options)
                arguments
                    obj
                    options.title = []
                    options.fig = []
                    options.font_size = 12
                    options.quantity = 'energy' % options: 'jacobi', 'energy', 'hamiltonian', 'kamiltonian'
                    options.show_steps = false
                    options.color = []
                end
    
                [dVal, label] = obj.get_conserved_fluctuation(options.quantity);
    
                if isempty(options.fig)
                    figure;
                else
                    figure(options.fig);
                end
                hold on;
                sgtitle(options.title, 'Interpreter', 'latex');

                if options.show_steps
                    name_plot = obj.name + " $(n = " + num2str(obj.sol.Nsteps) + ")$";
                else
                    name_plot = obj.name;
                end
                
                if isempty(options.color)
                    plot(obj.sol.t/obj.prob.Tp, dVal,LineWidth=2,DisplayName=name_plot);
                else
                    plot(obj.sol.t/obj.prob.Tp,dVal,LineWidth=2,DisplayName=name_plot,...
                    Color=options.color);
                end

                xlabel('Revolutions', 'Interpreter', 'latex');
                ylabel(label, 'Interpreter', 'latex');
                set(gca, 'yscale', 'log'); %grid on;
                legend('Interpreter', 'latex', 'Location','best');
                set(findall(gcf, '-property', 'FontSize'), 'FontSize', options.font_size);
    
            end
    
            function plot_traj_with_drift(obj,options)
                arguments
                    obj
                    options.colormap = cool
                    options.one_rev = true
                    options.title = obj.name 
                    options.fig = []
                    options.font_size = 12;
                    options.quantity = 'energy' % options: 'jacobi', 'energy', 'hamiltonian', 'kamiltonian'
                end
    
                [dVal, label] = obj.get_conserved_fluctuation(options.quantity);
    
                % Slice to one revolution if requested
                if options.one_rev
                    idx = obj.sol.t <= obj.prob.Tp;  % take only points in first rev
                else
                    idx = true(1, length(obj.sol.t));  % all points
                end
    
                % Extract data for plotting
                t_plot = obj.sol.t(idx);
                dVal_plot = dVal(idx);
                x_plot = obj.sol.x(:, idx);
    
                if isempty(options.fig)
                    figure;
                else
                    figure(options.fig);
                end
                sgtitle(options.title, 'Interpreter', 'latex');
    
                subplot(1,2,1)
                scatter(t_plot/obj.prob.Tp,dVal_plot,5,dVal_plot,'filled')
                set(gca,'yscale','log'); %grid on;
                set(gca, 'ColorScale', 'log');
                xlabel('Revolutions','Interpreter','latex')
                ylabel(label,'Interpreter','latex')
    
                subplot(1,2,2); hold on; %grid on;
                scatter3(x_plot(1,:), x_plot(2,:), x_plot(3,:), 5, dVal_plot, 'filled')
                scatter3(obj.prob.DS.r2,0,0,15,'k','filled')
                str = 'Moon'; text(obj.prob.DS.r2,-0.1,str)
                colormap(options.colormap)
                cb = colorbar;
                set(gca, 'ColorScale', 'log');
                % clim([min(dVal_plot) max(dVal_plot)])
                ylabel(cb,label,'Rotation',270,'Interpreter','latex')
                axis equal; view(2); labels3d('LU')
    
                set(findall(gcf, '-property', 'FontSize'), 'FontSize', options.font_size);
            end
    
            function [dVal, label] = get_conserved_fluctuation(obj, quantity)
                % Returns the fluctuation |Q - Q0| and the appropriate label string
                % based on the conserved quantity requested.

                switch lower(quantity)
                    case 'jacobi'
                        if ~isfield(obj.sol, 'C')
                            obj.sol.C = obj.prob.DS.jacobiconstant(obj.sol);
                        end
                        val = obj.sol.C;
                        label = '$\delta C$'; %'$|C - C_0|$';
                    case 'energy'
                        if ~isfield(obj.sol, 'E')
                            obj.sol.E = obj.prob.DS.total_energy(obj.sol.x);
                        end
                        val = obj.sol.E;
                        label = '$|E - E_0|$';
                    case 'hamiltonian'
                        if ~isfield(obj.sol, 'H')
                            obj.sol.H = obj.prob.DS.hamiltonian(obj.sol.x);
                        end
                        val = obj.sol.H;
                        label = '$|H - H_0|$';
                    case {'kamiltonian', 'kamilt'}
                        if ~isfield(obj.sol, 'K')
                            obj.sol.K = obj.prob.DS.kamiltonian(obj.sol.x);
                        end
                        val = obj.sol.K;
                        label = '$|K - K_0|$';
                    otherwise
                        error('Unknown quantity "%s". Choose "jacobi", "energy", or "kamiltonian".', quantity);
                end

                % dVal = abs(val - val(1));
                dVal = abs((val - val(1))/val(1));
            end

            function plot_state_error(obj,state_error,options)
                arguments
                    obj
                    state_error
                    options.title = []
                    options.fig = []
                    options.font_size = 12
                    options.show_steps = false
                    options.color = []
                end

                if isempty(options.fig)
                    figure;
                else
                    figure(options.fig);
                end
                hold on
                sgtitle(options.title, 'Interpreter', 'latex');

                if options.show_steps
                    name_plot = obj.name + " $(n = " + num2str(obj.sol.Nsteps) + ")$";
                else
                    name_plot = obj.name;
                end

                plot(obj.sol.t/obj.prob.Tp, state_error,LineWidth=2,DisplayName=name_plot);

                xlabel('Revolutions', 'Interpreter', 'latex');
                ylabel('State error', 'Interpreter', 'latex');
                set(gca, 'yscale', 'log'); %grid on;
                legend('Interpreter', 'latex', 'Location','best');
                set(findall(gcf, '-property', 'FontSize'), 'FontSize', options.font_size);
   
            end


        end
    end