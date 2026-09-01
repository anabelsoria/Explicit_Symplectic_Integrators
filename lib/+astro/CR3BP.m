%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Author : Anabel Soria-Carro
% Date   : May 22, 2025
% Affiliation: The University of Texas at Austin
%              Controls Group for Distributed and Uncertain Systems (CDUS)
% Description:
%  This class defines the Circular Restricted Three-Body Problem (CR3BP)
%  dynamical system.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

classdef CR3BP < astro.DynamicalSystem

    % These properties will be set by the user
    properties
        center  string % Origin coordinates
    end

    % These properties are internally set by the class
    properties
        r1      % Position of larger primary
        r2      % Position of secondary primary
    end

    properties (Constant)
        % From xi = (r,v) to nu = (q,p)
        P_xi_nu = [eye(3), zeros(3);
            0 -1 0 1 0 0;
            1  0 0 0 1 0;
            0  0 0 0 0 1];

        P_nu_xi = [eye(3),  zeros(3);
            0  1 0   1 0 0;
            -1  0 0   0 1 0;
            0  0 0   0 0 1];
    end

    methods
        function obj = CR3BP(mu,center,LU,TU)
            obj@astro.DynamicalSystem();
            obj.mu = mu;
            obj.LU = LU;
            obj.TU = TU;
            obj.center = center;
        end

        function ds = EOM(obj, t, s)
            n = length(s);
            planar = (n <= 4);

            mu1 = 1 - obj.mu;
            mu2 = obj.mu;

            if planar
                x = s(1); y = s(2);
                xDot = s(3); yDot = s(4);

                r13 = ((x - obj.r1)^2 + y^2)^1.5;
                r23 = ((x - obj.r2)^2 + y^2)^1.5;

                Ux = x - mu1 * (x - obj.r1) / r13 - mu2 * (x - obj.r2) / r23;
                Uy = y - mu1 * y / r13 - mu2 * y / r23;

                xDDot = 2 * yDot + Ux;
                yDDot = -2 * xDot + Uy;

                ds = [xDot; yDot; xDDot; yDDot];
            else
                x = s(1); y = s(2); z = s(3);
                xDot = s(4); yDot = s(5); zDot = s(6);

                r13 = (x - obj.r1)^2 + y^2 + z^2;
                r23 = (x - obj.r2)^2 + y^2 + z^2;

                r13 = r13^1.5;
                r23 = r23^1.5;

                Ux = x - mu1 * (x - obj.r1) / r13 - mu2 * (x - obj.r2) / r23;
                Uy = y - mu1 * y / r13 - mu2 * y / r23;
                Uz = -mu1 * z / r13 - mu2 * z / r23;

                xDDot = 2 * yDot + Ux;
                yDDot = -2 * xDot + Uy;
                zDDot = Uz;

                ds = [xDot; yDot; zDot; xDDot; yDDot; zDDot];
            end
        end

        function dh = Hamiltons_EOM(obj, t, s)
            x = s(1); y = s(2); z = s(3);
            px = s(4); py = s(5); pz = s(6);

            mu1 = 1 - obj.mu;
            mu2 = obj.mu;

            r13 = sqrt((x - obj.r1)^2 + y^2 + z^2);
            r23 = sqrt((x - obj.r2)^2 + y^2 + z^2);

            % Ux = mu1 * (x - obj.r1) / r13^3 + mu2 * (x - obj.r2) / r23^3;
            % Uy = mu1 * y / r13^3       + mu2 * y / r23^3;
            % Uz = mu1 * z / r13^3       + mu2 * z / r23^3;
            dUdx = partialU(obj, s(1:3));
            Ux = dUdx(1); Uy = dUdx(2); Uz = dUdx(3);

            dh = zeros(6, 1);
            dh(1) = px + y;
            dh(2) = py - (x + (1 - obj.mu - obj.r2)); % Last part is zero when bary
            dh(3) = pz;
            dh(4) = -( -py + Ux );
            dh(5) = -(  px + Uy );
            dh(6) = -Uz;
        end

        function dUdx = partialU(obj, q)
            mu1 = 1 - obj.mu; mu2 = obj.mu;

            x = q(1); y = q(2);
            z = 0;
            if length(q) > 2
                z = q(3);
            end

            r13 = sqrt((x - obj.r1)^2 + y^2 + z^2);
            r23 = sqrt((x - obj.r2)^2 + y^2 + z^2);

            Ux = mu1 * (x - obj.r1) / r13^3 + mu2 * (x - obj.r2) / r23^3;
            Uy = mu1 * y / r13^3 + mu2 * y / r23^3;
            Uz = mu1 * z / r13^3 + mu2 * z / r23^3;

            dUdx = [Ux; Uy];
            if length(q) > 2
                dUdx = [dUdx; Uz];
            end
        end

        function Hp = Hp(obj,s)
            q = s(1:3); p = s(4:end);

            Hp = [p(1) + q(2);...
                p(2) - (1-obj.mu - obj.r2)- q(1);...
                p(3)];

            % Hp = [p1+q2,obj.mu+p2-q1-1.0,p3]';
        end

        function U = potentialU(obj,q)
            x = q(1,:); y = q(2,:);
            z = 0;
            if size(q,1) > 2
                z = q(3,:);
            end

            r13 = sqrt((x - obj.r1).^2 + y.^2 + z.^2);
            r23 = sqrt((x - obj.r2).^2 + y.^2 + z.^2);

            U = (1-obj.mu)./r13 + obj.mu./r23;

        end

        function [q_n1, p_n1, e_dt] = SI_EOM(obj, dt, scheme, X,e_dt)
            q = X(1:3);
            p = X(4:6);

            dt2 = (dt) / 2;
            den = 1 + (dt2^2);

            % [den_compsum, e_dt] = comp_sum(1,e_dt,vpa(dt2^2));
            % den = den_compsum;

            % T = (1/den) * [1, dt2, 0;
            %     -dt2, 1, 0;
            %     0,  0, den];
            %
            % D = [1, dt2, 0;
            %     -dt2, 1,  0;
            %     0,   0,  1];

            %T = double(T); D = double(D);

            x = dt2;
            absx = abs(x);
            if absx < 1e-3
                % Use Taylor series for a = 1/(1+x^2) and b = x/(1+x^2)
                % a = 1 - x^2 + x^4 - x^6 + ...
                % b = x*(1 - x^2 + x^4 - ...)
                x2 = x*x;
                x4 = x2*x2;
                x6 = x4*x2;
                % keep up to x^4 terms (usually enough for tiny x)
                a = 1 - x2 + x4 -x6;
                b = x*(1 - x2 + x4 -x6);
            else
                % Direct stable computation
                den = 1 + x*x;    % safe here because x not extremely tiny
                a = 1.0 / den;
                b = x * a;
            end

            % Build T and D
            T = [ a,  b, 0;
                -b,  a, 0;
                0,  0, 1];

            D = [1, x, 0;
                -x,1,  0;
                0, 0,  1];


            switch scheme
                case 1 % Stormer-Verlet A
                    if ~strcmp(obj.center,"bary")
                        error("Scheme 1 not implemented with center thats not bary")
                    end
                    p_n2 = T * (p - dt2 * obj.partialU(q));
                    q_n1 = T * (D * q + dt * p_n2);
                    p_n1 = D * p_n2 - dt2 * obj.partialU(q_n1);

                case 2 % Stormer-Verlet B
                    % q_n2 = T * (q + dt2 * (p - [0; 1-obj.mu - obj.r2;0]));
                    % p_n1 = T * (D * p - dt * obj.partialU(q_n2));
                    % q_n1 = D * q_n2 + dt2 * (p_n1 - [0; 1-obj.mu - obj.r2;0]);

                    x = dt / 2;                  % dt2
                    shift = [0; 1 - obj.mu - obj.r2; 0];
                    one_dd = obj.dd_from_double(1.0);
                    x_dd   = obj.dd_from_double(x);

                    % --- den = 1 + x^2 (dd)
                    x2_dd  = obj.dd_mul(x_dd, x_dd);
                    den_dd = obj.dd_add(one_dd, x2_dd);

                    % --- a = 1 / den (dd), b = x / (1 + x^2) (dd)
                    a_dd = obj.dd_recip(den_dd);
                    b_dd = obj.dd_mul(x_dd, a_dd);

                    % ==========================================================
                    % q_n2 = T * (q + dt2 * (p - shift))
                    % ==========================================================
                    tmp = q + x * (p - shift);
                    q_n2 = zeros(3,1);

                    % first two components via dd rotation
                    q_n2(1) = obj.dd2double(obj.dd_add(obj.dd_mul(a_dd, obj.dd_from_double(tmp(1))), ...
                        obj.dd_mul(b_dd, obj.dd_from_double(tmp(2)))));
                    q_n2(2) = obj.dd2double(obj.dd_add(obj.dd_mul(obj.dd_neg(b_dd), obj.dd_from_double(tmp(1))), ...
                        obj.dd_mul(a_dd, obj.dd_from_double(tmp(2)))));
                    q_n2(3) = tmp(3);  % third unchanged

                    % ==========================================================
                    % p_n1 = T * (D * p - dt * gradU(q_n2))
                    % ==========================================================
                    grad = obj.partialU(q_n2);

                    % D*p
                    Dp = [p(1) + x*p(2);
                        -x*p(1) + p(2);
                        p(3)];

                    v = Dp - dt * grad;

                    p_n1 = zeros(3,1);
                    p_n1(1) = obj.dd2double(obj.dd_add(obj.dd_mul(a_dd, obj.dd_from_double(v(1))), ...
                        obj.dd_mul(b_dd, obj.dd_from_double(v(2)))));
                    p_n1(2) = obj.dd2double(obj.dd_add(obj.dd_mul(obj.dd_neg(b_dd), obj.dd_from_double(v(1))), ...
                        obj.dd_mul(a_dd, obj.dd_from_double(v(2)))));
                    p_n1(3) = v(3);

                    % ==========================================================
                    % q_n1 = D * q_n2 + dt2 * (p_n1 - shift)
                    % ==========================================================
                    Dq_n2 = [q_n2(1) + x*q_n2(2);
                        -x*q_n2(1) + q_n2(2);
                        q_n2(3)];

                    q_n1 = Dq_n2 + x * (p_n1 - shift);

            end

            %q_n1 = double(q_n1); p_n1 = double(p_n1);

        end

        function C = jacobiconstant(obj, sol)
            s = sol.x;
            if strcmp(sol.coord,'hamiltonian')
                s = obj.nu2xi(s);
            end
            if size(s,1) > 4
                x = s(1,:); y = s(2,:); z = s(3,:);
                v = s(4:end,:);
            else
                x = s(1,:); y = s(2,:); z = 0;
                v = s(3:end,:);
            end

            Ustar = obj.potentialU(s(1:3,:)) + ...
                1/2 * ((x + (1 - obj.mu - obj.r2)).^2 + y.^2);
            C = 2*Ustar - sum(v.^2); %vecnorm(v).^2;
        end

        function H = hamiltonian(obj,sol)
            s = sol.x;
            if strcmp(sol.coord,'cartesian')
                s = obj.xi2nu(s);
            end
            q1 = s(1,:); q2 = s(2,:); q3 = s(3,:);
            p1 = s(4,:); p2 = s(5,:); p3 = s(6,:);

            Ustar = obj.potentialU(s(1:3,:)) + 1/2*((q1 + ...
                (1 - obj.mu - obj.r2)).^2 + q2.^2);
            H = 1/2 * ((p1 + q2).^2 + (p2 - (q1+(1 - obj.mu - obj.r2)) ).^2 + p3.^2)...
                -Ustar;
        end

        function xi = nu2xi(obj,nu)
            switch obj.center
                case 'bary'
                    xi = obj.P_nu_xi * nu;
                case 'p2'
                    xi = obj.P_nu_xi * (nu + [1-obj.mu;0;0;0;0;0]);
                    xi(1,:) = xi(1,:) - (1-obj.mu);
                case 'p1'
                    xi = obj.P_nu_xi * (nu + [-obj.mu;0;0;0;0;0]);
                    xi(1,:) = xi(1,:) - (-obj.mu);
            end
        end

        function nu = xi2nu(obj,xi)
            switch obj.center
                case 'bary'
                    nu = obj.P_xi_nu * xi;
                case 'p2'
                    nu = obj.P_xi_nu * (xi + [1-obj.mu;0;0;0;0;0]);
                    nu(1,:) = nu(1,:) - (1-obj.mu);
                case 'p1'
                    nu = obj.P_xi_nu * (xi + [-obj.mu;0;0;0;0;0]);
                    nu(1,:) = nu(1,:) - (-obj.mu);
            end
        end
    end

    methods (Static)
        % =========================
        % --- DOUBLE-DOUBLE ARITHMETIC ---
        % =========================
        function dd = dd_from_double(x)
            dd = [x, 0.0];
        end

        function d = dd2double(dd)
            d = dd(1) + dd(2);
        end

        function dd = dd_add(a, b)
            [s, e1] = astro.CR3BP.twoSum(a(1), b(1));
            e = a(2) + b(2) + e1;
            [hi, lo] = astro.CR3BP.quickTwoSum(s, e);
            dd = [hi, lo];
        end

        function dd = dd_sub(a, b)
            [s, e1] = astro.CR3BP.twoSum(a(1), -b(1));
            e = a(2) - b(2) + e1;
            [hi, lo] = astro.CR3BP.quickTwoSum(s, e);
            dd = [hi, lo];
        end

        function dd = dd_mul(a, b)
            [p1, e1] = astro.CR3BP.twoProd(a(1), b(1));
            t = a(1)*b(2) + a(2)*b(1);
            s = p1 + t;
            e = e1 + (t - (s - p1)) + a(2)*b(2);
            [hi, lo] = astro.CR3BP.quickTwoSum(s, e);
            dd = [hi, lo];
        end

        function dd = dd_neg(a)
            dd = [-a(1), -a(2)];
        end

        function dd = dd_recip(d)
            r0 = 1.0 / d(1);
            r = astro.CR3BP.dd_from_double(r0);
            one = astro.CR3BP.dd_from_double(1.0);
            dr = astro.CR3BP.dd_mul(d, r);
            tmp = astro.CR3BP.dd_sub(one, dr);
            corr = astro.CR3BP.dd_mul(r, tmp);
            r = astro.CR3BP.dd_add(r, corr);
            [hi, lo] = astro.CR3BP.quickTwoSum(r(1), r(2));
            dd = [hi, lo];
        end

        % --- ERROR-FREE TRANSFORMS ---
        function [s, e] = twoSum(a, b)
            s = a + b;
            bp = s - a;
            e = (a - (s - bp)) + (b - bp);
        end

        function [s, e] = quickTwoSum(a, b)
            s = a + b;
            e = b - (s - a);
        end

        function [p, e] = twoProd(a, b)
            c = 134217729 * a; % 2^27+1
            a_hi = c - (c - a);
            a_lo = a - a_hi;
            c = 134217729 * b;
            b_hi = c - (c - b);
            b_lo = b - b_hi;
            p = a * b;
            e = ((a_hi * b_hi - p) + a_hi*b_lo + a_lo*b_hi) + a_lo*b_lo;
        end
    end
end
