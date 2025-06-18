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

        end

        function plot_dC(obj,options)
            arguments
                obj
                options.title = []
                options.fig = []
            end

            if ~isfield(obj.sol,'dC')
                obj.sol.C = obj.prob.DS.jacobiconstant(obj.sol);
                obj.sol.dC = abs(obj.sol.C-obj.sol.C(1));
            end

            if isempty(options.fig)
                figure;
            else
                figure(options.fig);
            end
            hold on;
            sgtitle(options.title, 'Interpreter', 'latex');
            
            plot(obj.sol.t/obj.prob.Tp,obj.sol.dC,LineWidth=2,DisplayName=obj.name)

            xlabel('Revolutions', 'Interpreter', 'latex');                
            ylabel('$|C - C_0|$', 'Interpreter', 'latex');                
            set(gca, 'yscale', 'log'); grid on;                                    
            legend('Interpreter', 'latex', 'Location','best'); 

        end

        function plot_traj_dCfluct(obj,options)
            arguments
                obj
                options.colormap = cool
                options.one_rev = true
                options.title = []
                options.fig = []
            end

            if ~isfield(obj.sol,'dC')
                obj.sol.C = obj.prob.DS.jacobiconstant(obj.sol);
                obj.sol.dC = abs(obj.sol.C-obj.sol.C(1));
            end

            % Slice to one revolution if requested
            if options.one_rev
                idx = obj.sol.t <= obj.prob.Tp;  % take only points in first rev
            else
                idx = true(1, length(obj.sol.t));  % all points
            end

            % Extract data for plotting
            t_plot = obj.sol.t(idx);
            dC_plot = obj.sol.dC(idx);
            x_plot = obj.sol.x(:, idx);

            if isempty(options.fig)
                figure;
            else
                figure(options.fig);
            end
            sgtitle(options.title, 'Interpreter', 'latex'); 

            subplot(1,2,1)
            scatter(t_plot/obj.prob.Tp,dC_plot,5,dC_plot,'filled')
            set(gca,'yscale','log'); grid on;
            xlabel('Revolutions','Interpreter','latex')
            ylabel('$\delta C$','Interpreter','latex')

            subplot(1,2,2); hold on; grid on;
            scatter3(x_plot(1,:), x_plot(2,:), x_plot(3,:), 5, dC_plot, 'filled')
            scatter3(obj.prob.DS.r2,0,0,15,'k','filled')
            str = 'Moon'; text(obj.prob.DS.r2,-0.1,str)
            colormap(options.colormap)
            cb = colorbar;
            clim([min(dC_plot) max(dC_plot)])
            ylabel(cb,'$\delta C$','Rotation',270,'Interpreter','latex')
            axis equal; view(2); labels3d('LU')
        end
    end
end