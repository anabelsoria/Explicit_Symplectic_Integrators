classdef TimeRegularized < Integrator

    properties
        integrator     % handle to base integrator (RK, SI, etc.)
        reg_fun        % time regularization function (handle)
        params         % parameters needed
        method         % time regularization method (name)
    end

    methods
        function obj = TimeRegularized(integrator, method, reg_params)
            obj.integrator = integrator;
            obj.prob       = integrator.prob;
            obj.params     = reg_params;
            obj.name       = integrator.name + " " + method;
            % obj.order       = integrator.order;
            % obj.scheme       = integrator.scheme;
            obj.method     = method;
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
                    reg_fun = @(q, p) obj.russell_regularization(q, p);
                    % if ~strcmp(obj.integrator.prob.center,'p2')
                    %     error("Russell regularization requires the center to be at p2, but current center is " + obj.integrator.prob.center)
                    % end

                case 'heggie'
                    reg_fun = @(q,p) obj.heggie_regularization(q,p);

                otherwise
                    error('Unknown time regularization method.');
            end
        end

        function [z, G] = sundman_regularization(obj,q, p, alpha)
            z = 1 / (q(1:3)'*q(1:3))^(alpha/2);
            Hp = obj.integrator.prob.DS.Hp([q;p]);
            G = -alpha * Hp(1:3)' * q(1:3) / (q(1:3)'*q(1:3));
        end

        function [z, G] = russell_regularization(obj,q, p)

            A_E = obj.params.A_E; B_E = obj.params.B_E; C_E = obj.params.C_E;
            rH_E = obj.params.rH_E; alpha_E = obj.params.alpha_E;

            A_M = obj.params.A_M; B_M = obj.params.B_M; C_M = obj.params.C_M;
            rH_M = obj.params.rH_M; alpha_M = obj.params.alpha_M;

            rE = norm(q(1:3)-[-1,0,0]');
            rM = norm(q(1:3));

            rho_E = obj.rho_fun(A_E,B_E,C_E,rE,rH_E);
            rho_M = obj.rho_fun(A_M,B_M,C_M,rM,rH_M);

            g = rho_E^alpha_E * rho_M^alpha_M;

            z = 1/g;

            G = obj.G_Russell(q,p);
        end

        function rho = rho_fun(~,A,B,C,r,rH)
            rho = (A^2 + A) / ( A + (A*(1-C)/(C+A))^(r/(B*rH)) ) - A;
        end

        function G = G_Russell(obj,q,p)
            A_E = obj.params.A_E; B_E = obj.params.B_E; C_E = obj.params.C_E;
            rH_E = obj.params.rH_E; alpha_E = obj.params.alpha_E;

            A_M = obj.params.A_M; B_M = obj.params.B_M; C_M = obj.params.C_M;
            rH_M = obj.params.rH_M; alpha_M = obj.params.alpha_M;

            mu = obj.prob.DS.mu;

            q1 = q(1); q2 = q(2); q3 = q(3);
            p1 = p(1); p2 = p(2); p3 = p(3);

            t2 = abs(q1);
            t3 = abs(q2);
            t4 = abs(q3);
            t5 = sign(q2);
            t6 = sign(q3);
            t7 = A_E.^2;
            t8 = A_M.^2;
            t9 = q1+1.0;
            t10 = A_E+C_E;
            t11 = A_M+C_M;
            t16 = -A_E;
            t17 = -A_M;
            t18 = 1.0./B_E;
            t19 = 1.0./B_M;
            t20 = C_E-1.0;
            t21 = C_M-1.0;
            t22 = -alpha_E;
            t23 = -alpha_M;
            t24 = alpha_E-1.0;
            t25 = alpha_M-1.0;
            t26 = 1.0./rH_E;
            t27 = 1.0./rH_M;
            t12 = abs(t9);
            t13 = t2.^2;
            t14 = t3.^2;
            t15 = t4.^2;
            t28 = A_E+t7;
            t29 = A_M+t8;
            t30 = 1.0./t10;
            t31 = 1.0./t11;
            t32 = t12.^2;
            t35 = t13+t14+t15;
            t36 = t16.*t20.*t30;
            t37 = t17.*t21.*t31;
            t38 = t14+t15+t32;
            t39 = log(t36);
            t40 = log(t37);
            t43 = sqrt(t35);
            t41 = conj(t39);
            t42 = conj(t40);
            t44 = 1.0./t43;
            t45 = sqrt(t38);
            t47 = t19.*t27.*t43;
            t46 = 1.0./t45;
            t48 = t18.*t26.*t45;
            t49 = t37.^t47;
            t50 = conj(t49);
            t51 = t36.^t48;
            t52 = A_M+t49;
            t53 = A_M+t50;
            t54 = conj(t51);
            t55 = A_E+t51;
            t56 = 1.0./t52;
            t57 = A_E+t54;
            t58 = 1.0./t53.^2;
            t59 = 1.0./t55;
            t61 = t29.*t56;
            t60 = 1.0./t57.^2;
            t62 = t28.*t59;
            t63 = t17+t61;
            t64 = t16+t62;
            t65 = t63.^alpha_M;
            t69 = t63.^t23;
            t70 = t63.^t25;
            t66 = conj(t65);
            t67 = t64.^alpha_E;
            t71 = conj(t70);
            t72 = t64.^t22;
            t73 = t64.^t24;
            t68 = conj(t67);
            t74 = conj(t73);
            G = t69.*t72.*(alpha_E.*t3.*t5.*t18.*t26.*t28.*t41.*t46.*t54.*t60.*t66.*t74+alpha_M.*t3.*t5.*t19.*t27.*t29.*t42.*t44.*t50.*t58.*t68.*t71).*(mu+p2-t9)+t69.*t72.*(p1+q2).*(alpha_M.*t2.*t19.*t27.*t29.*t42.*t44.*t50.*t58.*t68.*t71.*sign(q1)+alpha_E.*t12.*t18.*t26.*t28.*t41.*t46.*t54.*t60.*t66.*t74.*sign(t9))+p3.*t69.*t72.*(alpha_E.*t4.*t6.*t18.*t26.*t28.*t41.*t46.*t54.*t60.*t66.*t74+alpha_M.*t4.*t6.*t19.*t27.*t29.*t42.*t44.*t50.*t58.*t68.*t71);
        end

        function [z, G] = heggie_regularization(obj,q,p)
            r12 = norm(q(1:3)-[obj.r1,0,0]'); % s/c - primary
            r13 = norm(q(1:3)-[obj.r2,0,0]'); % s/c - secondary
            r23 = norm([obj.r1,0,0]' - [obj.r2,0,0]'); % primary - secondary

            g = r12*r13*r23 / (r12 + r13 + r23)^(3/2);

            z = 1/g;

            % TODO: diff_g_heggie(q,r2,r3) only ever worked in the removed
            % HR4BP branch -- its r2/r3 args (HR4BP primary position
            % vectors) were never defined here, and diff_g_heggie.m itself
            % has been removed. There's a commented-out hand-expanded
            % formula below (diff_g_q1/q2/q3) that looks like an attempt
            % at a non-HR4BP version, using the same r2/r3 naming -- worth
            % checking whether that's meant to use [obj.r1,0,0]'/[obj.r2,0,0]'
            % in place of the HR4BP r2/r3 vectors before this is usable.
            error('heggie_regularization:notImplemented', ...
                'Heggie regularization is not implemented for non-HR4BP systems.');

            % diff_g_q1 = ((2*q(1) - 2*r2(1))*r13*r23)/(2*r12*(r12 + r13 + r23)^(3/2)) - (3*((q(1) - r2(1))/r12 + (q(1) - r3(1))/r13)*r12*r13*r23)/(2*(r12 + r13 + r23)^(5/2)) + ((2*q(1) - 2*r3(1))*r12*r23)/(2*r13*(r12 + r13 + r23)^(3/2));
            % diff_g_q2 = ((2*q(1) - 2*r2(2))*r13*r23)/(2*r12*(r12 + r13 + r23)^(3/2)) - (3*((q(1) - r2(2))/r12 + (q(1) - r3(2))/r13)*r12*r13*r23)/(2*(r12 + r13 + r23)^(5/2)) + ((2*q(1) - 2*r3(2))*r12*r23)/(2*r13*(r12 + r13 + r23)^(3/2));
            % diff_g_q3 = ((2*q(3) - 2*r2(3))*r13*r23)/(2*r12*(r12 + r13 + r23)^(3/2)) - (3*((q(3) - r2(3))/r12 + (q(3) - r3(3))/r13)*r12*r13*r23)/(2*(r12 + r13 + r23)^(5/2)) + ((2*q(3) - 2*r3(3))*r12*r23)/(2*r13*(r12 + r13 + r23)^(3/2));

            % Hp = obj.prob.DS.Hp([q;p])';
            % G = -1/g * [diff_g_q1; diff_g_q2; diff_g_q3]' * Hp(1:3)';
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

        function varargout = propagateTimeRegSI(obj, t0, tf, h)

            ns = length(obj.prob.nu0);
            nq = ns / 2;

            q = obj.integrator.prob.nu0(1:nq);
            p = obj.integrator.prob.nu0(nq + 1: end);
            t = t0;
            z0= obj.reg_fun(q, p);
            z = z0;
            X = [q; p];
            tspan = t0;
            ii = 1;

            tic
            while (h > 0 && t < tf) || (h < 0 && t > tf)

                %,rho1(ii),rho2(ii)
                [~,G] = obj.reg_fun(q, p);
                z = z + h * G;
                dt = h / z;
                for jj = 1:length(obj.integrator.gamma)
                    [q, p] = obj.integrator.prob.DS.SI_EOM_nomp(obj.integrator.gamma(jj)*dt, obj.integrator.scheme, [q; p]);
                    t = t + obj.integrator.gamma(jj)*dt;
                end

                ii = ii + 1;
                X(:, ii) = [q; p];
                tspan(ii) = t;
             
                if abs(X(1,ii)) > 10
                    disp('Exploted')
                    continue
                end

                if (h > 0 && any(diff(tspan)<0) ) || (h < 0 && any(diff(tspan)>0) )
                    t
                    dt
                    disp('time exploted')
                    continue
                end
            end
            obj.time_solver = toc;

            % Post process

            obj.sol.x = X;
            obj.sol.t = tspan;
            obj.sol.coord = 'hamiltonian';
            obj.sol.Nsteps = round(length(tspan)/obj.prob.Nrevs);
            % obj.sol.rho1 = rho1;
            % obj.sol.rho2 = rho2;

            if nargout > 0
                varargout{1} = X;
                varargout{2} = tspan;
                %varargout{3} = taus;
            end
        end

        function varargout = propagateTimeRegRK(obj,t0, tf, h)
            X = obj.integrator.prob.nu0;
            X = [X;t0];

            ns = length(obj.integrator.prob.nu0);

            i = 1;
            t = t0;
            tau(i) = 0;

            f = @(t,x)obj.regularized_EOM(t,x);

            tic
            while (t < tf && h > 0) || (t>tf && h < 0)
                %,k1(:,i),k2(:,i),dh(:,i),g(i),G(i),rE(i),rM(i)
                X(:, i+1)  = obj.integrator.stepFun(h, tau(i), X(:, i), f);
                t = X(end, i+1);
                tau(i+1) = tau(i) + h;
                i = i + 1;
            end
            obj.time_solver = toc;

            % Post process
            tspan = X(end,:);
            X = X(1:ns,:);

            obj.sol.x = X;
            obj.sol.t = tspan;
            obj.sol.coord = 'hamiltonian';
            obj.sol.Nsteps = round(length(tspan)/obj.prob.Nrevs);
            % obj.sol.k1 = k1; obj.sol.k2 = k2;
            % obj.sol.dh = dh; obj.sol.g = g;
            % obj.sol.G = G; obj.sol.rE = rE; obj.sol.rM = rM;

            if nargout > 0
                varargout{1} = X;
                varargout{2} = tspan;
                varargout{3} = tau;
            end
        end

        function [d_dtau] = regularized_EOM(obj,t,s)

            ns = length(obj.prob.nu0);
            nq = ns / 2;

            q = s(1:nq);
            p = s(nq+1:end);

            dh = obj.integrator.prob.DS.Hamiltons_EOM(t, s);

            [z,~]  = obj.reg_fun(q, p);

            g = 1/z;
            dh_dtau  = g*dh;
            dt_dtau = g;

            d_dtau = [dh_dtau;dt_dtau];

        end

        % ---------------- MEX files ----------------
        function [X, tspan] = propagate_mex(obj, t0, tf, h,CompSum)
            arguments
                obj
                t0
                tf
                h
                CompSum = false
            end

            if isa(obj.integrator, 'RK')
                [X, tspan] = obj.propagateTimeRegRK_mex(t0, tf, h);

            elseif isa(obj.integrator, 'SI')

                if CompSum
                    [X, tspan] = obj.propagateTimeRegSI_CompSum_mex(t0, tf, h);
                else
                    [X, tspan] = obj.propagateTimeRegSI_mex(t0, tf, h);
                end

            else
                error('Unsupported integrator type.');
            end
        end

        function varargout = propagateTimeRegSI_mex(obj, t0, tf,h)
            p = obj.prob;
            if isa(obj.prob.DS, 'astro.CR3BP')
                tic
                [X, tspan,coord] = SI_TR_CR3BP_propagate_mex(t0, tf, h, p.nu0,...
                    int32(obj.integrator.scheme),int32(obj.integrator.order),...
                    p.DS.mu, p.DS.r1, p.DS.r2,obj.params,char(obj.method));
                obj.time_solver = toc;
            elseif isa(obj.prob.DS, 'astro.ER3BP')
                tic
                [X, tspan,coord] = SI_TR_ER3BP_propagate_mex(t0, tf, h, p.nu0,...
                    int32(obj.integrator.scheme),int32(obj.integrator.order),...
                    p.DS.mu, p.DS.r1, p.DS.r2,p.DS.e,obj.params,char(obj.method));
                obj.time_solver = toc;
            end

            obj.sol.x = X;
            obj.sol.t = tspan;
            obj.sol.coord = coord;
            obj.sol.Nsteps = round(length(tspan)/obj.prob.Nrevs);

            if nargout > 0
                varargout{1} = X;
                varargout{2} = tspan;
                %varargout{3} = taus;
            end
        end

        function varargout = propagateTimeRegSI_CompSum_mex(obj, t0, tf,h)
            p = obj.prob;
            GM1 = 1-p.DS.mu;
            GM2 = p.DS.mu;
            x1  = -p.DS.mu;
            x2  = 1-p.DS.mu;
            tic
            [X, tspan,coord] = SI_TR_CR3BP_CompSum_propagate_mex(t0, tf, h, p.nu0,...
                int32(obj.integrator.order),p.DS.mu, ...
                GM1,GM2,x1,x2,obj.params,char(obj.method));
            obj.time_solver = toc;

            obj.sol.x = X;
            obj.sol.t = tspan;
            obj.sol.coord = coord;
            obj.sol.Nsteps = round(length(tspan)/obj.prob.Nrevs);

            if nargout > 0
                varargout{1} = X;
                varargout{2} = tspan;
                %varargout{3} = taus;
            end
        end

        function varargout = propagateTimeRegRK_mex(obj, t0, tf,h)
            p = obj.prob;
            ns = length(p.nu0);
            if isa(obj.prob.DS, 'astro.CR3BP')
                tic
                [X,coord] = RK_TR_CR3BP_propagate_mex(t0, tf, h,p.nu0, ...
                    int32(obj.integrator.order),p.DS.mu, p.DS.r1, p.DS.r2, ...
                    obj.params,char(obj.method));
                obj.time_solver = toc;
            elseif isa(obj.prob.DS, 'astro.ER3BP')
                tic
                [X,coord] = RK_TR_ER3BP_propagate_mex(t0, tf, h,p.nu0, ...
                    int32(obj.integrator.order),p.DS.mu, p.DS.r1, p.DS.r2,p.DS.e,...
                    obj.params,char(obj.method));
                obj.time_solver = toc;
            end

            X = [X{:}];

            obj.sol.x = X(1:ns,:);
            obj.sol.t = X(end,:);
            obj.sol.coord = coord;
            obj.sol.Nsteps = round(length(obj.sol.t)/obj.prob.Nrevs);

            if nargout > 0
                varargout{1} = obj.sol.x;
                varargout{2} = obj.sol.t;
            end
        end
    end
end
