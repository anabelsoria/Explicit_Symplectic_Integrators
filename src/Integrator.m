    classdef Integrator < handle
    
        properties
            sol % Propagated solution (trajectory and tspan)
        end
    
        methods
    
            function plot_traj(obj,options)
                arguments
                    obj
                    options.title = []
                    options.fig = []
                    options.plot_2d_xy = false
                    options.plot_2d_xz = false
                    options.font_size = 12;
                end
    
                if isempty(options.fig)
                    figure;
                    scatter3(obj.prob.DS.r2,0,0,15,'k','filled',DisplayName='Moon')
                else
                    figure(options.fig);
                end
    
                hold on;
                sgtitle(options.title, 'Interpreter', 'latex');
    
                plot3(obj.sol.x(1,:),obj.sol.x(2,:),obj.sol.x(3,:),LineWidth=1,DisplayName=obj.name)
                if options.plot_2d_xy
                    view(2)
                elseif options.plot_2d_xz
                    view(2)
                else
                    view(3)
                end
    
                labels3d('LU'); grid on;
                legend('Interpreter', 'latex', 'Location','best');
                set(findall(gcf, '-property', 'FontSize'), 'FontSize', options.font_size);
    
            end
    
            function plot_conserved_quantity(obj,options)
                arguments
                    obj
                    options.title = []
                    options.fig = []
                    options.font_size = 12;
                    options.quantity = 'energy' % options: 'jacobi', 'energy', 'hamiltonian', 'kamiltonian'
                end
    
                [dVal, label] = obj.get_conserved_fluctuation(options.quantity);
    
                if isempty(options.fig)
                    figure;
                else
                    figure(options.fig);
                end
                hold on;
                sgtitle(options.title, 'Interpreter', 'latex');
    
                plot(obj.sol.t/obj.prob.Tp,dVal,LineWidth=2,DisplayName=obj.name)
    
                xlabel('Revolutions', 'Interpreter', 'latex');
                ylabel(label, 'Interpreter', 'latex');
                set(gca, 'yscale', 'log'); grid on;
                legend('Interpreter', 'latex', 'Location','best');
                set(findall(gcf, '-property', 'FontSize'), 'FontSize', options.font_size);
    
            end
    
            function plot_traj_with_drift(obj,options)
                arguments
                    obj
                    options.colormap = cool
                    options.one_rev = true
                    options.title = []
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
                set(gca,'yscale','log'); grid on;
                xlabel('Revolutions','Interpreter','latex')
                ylabel(label,'Interpreter','latex')
    
                subplot(1,2,2); hold on; grid on;
                scatter3(x_plot(1,:), x_plot(2,:), x_plot(3,:), 5, dVal_plot, 'filled')
                scatter3(obj.prob.DS.r2,0,0,15,'k','filled')
                str = 'Moon'; text(obj.prob.DS.r2,-0.1,str)
                colormap(options.colormap)
                cb = colorbar;
                clim([min(dVal_plot) max(dVal_plot)])
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
                        label = '$|C - C_0|$';
                    case 'energy'
                        if ~isfield(obj.sol, 'E')
                            obj.sol.E = obj.prob.DS.energy(obj.sol.x);
                        end
                        val = obj.sol.E;
                        label = '$|E - E_0|$';
                    case {'kamiltonian', 'kamilt'}
                        if ~isfield(obj.sol, 'K')
                            obj.sol.K = obj.prob.DS.kamiltonian(obj.sol.x);
                        end
                        val = obj.sol.K;
                        label = '$|K - K_0|$';
                    otherwise
                        error('Unknown quantity "%s". Choose "jacobi", "energy", or "kamiltonian".', quantity);
                end

                dVal = abs(val - val(1));
            end


        end
    end