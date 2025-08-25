classdef TimeRegularized < Integrator

    properties
        integrator     % handle to base integrator (RK, SI, etc.)
        reg_fun        % time regularization function (handle)
        params         % parameters needed 
    end

    methods
        function obj = TimeRegularized(integrator, method, reg_params)
            obj.integrator = integrator;
            obj.prob       = integrator.prob;
            obj.params     = reg_params;
            obj.name       = integrator.name;
            % obj.order       = integrator.order;
            % obj.scheme       = integrator.scheme;
            obj.reg_fun    = obj.select_time_regularization(method);
            
        end

        function reg_fun = select_time_regularization(obj,method)
            switch lower(method)
                case 'sundman'
                    alpha = obj.params.alpha;
                    reg_fun = @(q, p) obj.sundman_regularization(q, p, alpha);
                    if ~strcmp(obj.integrator.prob.center,'p2')
                        error("Sundman regularization requires the center to be at p2, but current center is " + obj.integrator.prob.center)
                    end

                case 'russell'
                    alpha = obj.params.alpha;
                    reg_fun = @(q, p) obj.sundman_regularization(q, p, alpha);
                otherwise
                    error('Unknown time regularization method.');
            end
        end

        function [z, G] = sundman_regularization(obj,q, p, alpha)
            z = 1 / (q(1:3)'*q(1:3))^(alpha/2);
            Hp = obj.integrator.prob.DS.Hp([q;p]);
            G = -alpha * Hp(1:3)' * q(1:3) / (q(1:3)'*q(1:3));
        end


        % --- Propagation method ---
        function [X, tspan] = propagate(obj, t0, tf, epsilon)
            if isa(obj.integrator, 'RK')  
                [X, tspan] = obj.propagateTimeRegRK(t0, tf, epsilon);
            elseif isa(obj.integrator, 'SI') 
                [X, tspan] = obj.propagateTimeRegSI(t0, tf, epsilon);
            else
                error('Unsupported integrator type.');
            end
        end


        function varargout = propagateTimeRegSI(obj, t0, tf, epsilon)
            ns = length(obj.prob.nu0);
            nq = ns / 2;
            
            q = obj.integrator.prob.nu0(1:nq);
            p = obj.integrator.prob.nu0(nq + 1: end);
            t = 0;
            [z0,~] = obj.reg_fun(q, p);
            z = z0;
            X = [q; p];
            tspan = t0;
            ii = 1;

            tic
            while t < tf
                
                [~,G] = obj.reg_fun(q, p);
                z = z + epsilon * G;
                dt = epsilon / z;

                for jj = 1:length(obj.integrator.gamma)
                    [q, p] = obj.integrator.prob.DS.SI_EOM(obj.integrator.gamma(jj)*dt, obj.integrator.scheme, [q; p]);
                    t = t + obj.integrator.gamma(jj)*dt;
                end

                ii = ii + 1;
                X(:, ii) = [q; p];
                tspan(ii) = t;
            end
            obj.time_solver = toc;
            
            % Post process
            
            obj.sol.x = X;
            obj.sol.t = tspan;
            obj.sol.coord = 'hamiltonian';
            obj.sol.Nsteps = round(length(tspan)/obj.prob.Nrevs);

            if nargout > 0
                varargout{1} = X;
                varargout{2} = tspan;
                %varargout{3} = taus;
            end
        end

        function varargout = propagateTimeRegRK(obj,t0, tf, epsilon)
            tauspan = t0:epsilon:tf;
            nt = length(tauspan);
            ns = length(obj.integrator.prob.nu0);

            X = zeros(ns, nt);
            X(1:ns, 1) = obj.integrator.prob.nu0;
            X(end, 1) = t0;

            i = 1;
            t = t0;
            tau(i) = 0;

            f = @(t,x)obj.regularized_EOM(t,x);

            tic
            while t < tf
                X(:, i+1) = obj.integrator.stepFun(epsilon, tau(i), X(:, i), f);
                t = X(end, i+1);
                tau(i+1) = tau(i) + epsilon;
                i = i + 1;
            end
            obj.time_solver = toc;
            
            % Post process
            tspan = X(end,:);
            X = X(1:6,:);

            obj.sol.x = X;
            obj.sol.t = tspan;
            obj.sol.coord = 'hamiltonian';
            obj.sol.Nsteps = round(length(tspan)/obj.prob.Nrevs);

            if nargout > 0
                varargout{1} = X;
                varargout{2} = tspan;
                varargout{3} = tau;
            end
        end

        function d_dtau = regularized_EOM(obj,t,s)

            ns = length(obj.prob.nu0);
            nq = ns / 2;

            q = s(1:nq);
            p = s(nq+1:end);

            dh = obj.integrator.prob.DS.Hamiltons_EOM(t, s);
            z  = obj.sundman_regularization(q, p, obj.params.alpha);

            dh_dtau  = 1/z*dh;
            dt_dtau = 1/z;

            d_dtau = [dh_dtau;dt_dtau];
  
        end


    end
end
