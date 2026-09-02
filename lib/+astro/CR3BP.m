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
            x = s(1); y = s(2);
            px = s(4); py = s(5); pz = s(6);

            dUdx = obj.partialU(s(1:3));
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

        function [q_n1, p_n1] = SI_EOM(obj, dt, scheme, X)
            q = X(1:3);
            p = X(4:6);

            dt2 = dt / 2;
            den = 1 + dt2^2;

            T = (1/den) * [1, dt2, 0;
                -dt2, 1, 0;
                0,  0, den];

            D = [1, dt2, 0;
                -dt2, 1,  0;
                0,   0,  1];

            switch scheme
                case 1 % Stormer-Verlet A
                    p_n2 = T * (p - dt2 * obj.partialU(q));
                    q_n1 = T * (D * q + dt * p_n2);
                    p_n1 = D * p_n2 - dt2 * obj.partialU(q_n1);

                case 2 % Stormer-Verlet B
                    q_n2 = T * (q + dt2 * (p - [0; 1-obj.mu - obj.r2;0]));
                    p_n1 = T * (D * p - dt * obj.partialU(q_n2));
                    q_n1 = D * q_n2 + dt2 * (p_n1 - [0; 1-obj.mu - obj.r2;0]);
            end

        end


        % Precision-ladder kernels for Scheme 2 (Stormer-Verlet B).

        function [x_n1,p_n1] = SI_EOM_Expanded(obj, dt, phi_l, x, p)
            % Uncompensated scalar form. Baseline for SI_EOM_CS.
            dt  = phi_l*dt;
            hdt = dt/2;
            dt2_4 = 4 + dt^2;

            %% x_n2(1) = (4*x(1) + 2*dt*p(1) + 2*dt*x(2) + dt^2*p(2)) / (dt^2+4)
            s1 = 4*x(1) + 2*dt*p(1);
            s2 = s1 + 2*dt*x(2);
            s3 = s2 + dt^2*p(2);
            x_n2_1 = s3/dt2_4;

            %% x_n2(2) = (-2*dt*x(1) - dt^2*p(1) + 4*x(2) + 2*dt*p(2)) / (dt^2+4)
            s1 = 4*x(2) + 2*dt*p(2);
            s2 = s1 - 2*dt*x(1);
            s3 = s2 - dt^2*p(1);
            x_n2_2 = s3/dt2_4;

            %% x_n2(3) = x(3) + hdt*p(3)
            x_n2_3 = x(3) + hdt*p(3);

            x_n2 = [x_n2_1; x_n2_2; x_n2_3];

            %% Force at x_n2
            dU = obj.partialU(x_n2);

            %% p_n1(1)
            s1 = dt^2*p(1) - 4*p(1);
            s2 = s1 - 4*dt*p(2);
            s3 = s2 + 4*dt*dU(1);
            s4 = s3 + 2*dt^2*dU(2);
            p_n1_1 = -s4/dt2_4;

            %% p_n1(2)
            s1 = 4*dt*p(1) - 4*p(2);
            s2 = s1 + dt^2*p(2);
            s3 = s2 - 2*dt^2*dU(1);
            s4 = s3 + 4*dt*dU(2);
            p_n1_2 = -s4/dt2_4;

            %% p_n1(3) = p(3) - dt*dU(3)
            p_n1_3 = p(3) - dt*dU(3);

            %% x_n1(1) = x_n2(1) + hdt*x_n2(2) + hdt*p_n1(1)
            s1 = x_n2_1 + hdt*x_n2_2;
            x_n1_1 = s1 + hdt*p_n1_1;

            %% x_n1(2) = x_n2(2) - hdt*x_n2(1) + hdt*p_n1(2)
            s1 = x_n2_2 - hdt*x_n2_1;
            x_n1_2 = s1 + hdt*p_n1_2;

            %% x_n1(3) = x_n2(3) + hdt*p_n1(3)
            x_n1_3 = x_n2_3 + hdt*p_n1_3;

            x_n1 = [x_n1_1; x_n1_2; x_n1_3];
            p_n1 = [p_n1_1; p_n1_2; p_n1_3];
        end

        function [x,p] = SI_EOM_Increment(obj, dt, phi_l, x, p)
            % Increment form of SI_EOM_ICS, uncompensated.
            dt    = phi_l*dt;
            hdt   = dt/2;
            dt2_4 = 4 + dt^2;

            %% Stage 1, x <- x_n2
            dx1 = dt*(2*p(1) + 2*x(2) + dt*p(2) - dt*x(1))/dt2_4;
            dx2 = dt*(2*p(2) - 2*x(1) - dt*p(1) - dt*x(2))/dt2_4;
            dx3 = hdt*p(3);

            x(1) = x(1) + dx1;
            x(2) = x(2) + dx2;
            x(3) = x(3) + dx3;

            %% Force at x_n2
            dU = obj.partialU(x);

            %% Stage 2, p <- p_n1
            dp1 = dt*(-2*dt*p(1) + 4*p(2) - 4*dU(1) - 2*dt*dU(2))/dt2_4;
            dp2 = dt*(-4*p(1) - 2*dt*p(2) + 2*dt*dU(1) - 4*dU(2))/dt2_4;
            dp3 = -dt*dU(3);

            p(1) = p(1) + dp1;
            p(2) = p(2) + dp2;
            p(3) = p(3) + dp3;

            %% Stage 3, x <- x_n1
            dxx1 = hdt*(x(2) + p(1));
            dxx2 = hdt*(p(2) - x(1));
            dxx3 = hdt*p(3);

            x(1) = x(1) + dxx1;
            x(2) = x(2) + dxx2;
            x(3) = x(3) + dxx3;
        end

        function [x,p,e_x,e_p] = SI_EOM_ICS(obj, dt, phi_l, x, p, e_x, e_p)
            % Increment + Kahan-compensated summation (9 adds/step).
            % Center-agnostic via obj.r2, same as SI_EOM: shift s=1-mu-r2
            % is 0 at bary, 1-mu at p2.
            dt  = phi_l*dt;
            hdt = dt/2;
            dt2_4 = 4 + dt^2;
            s = 1 - obj.mu - obj.r2;

            %% Stage 1, x <- x_n2
            dx1 = dt*(2*p(1) + 2*x(2) + dt*p(2) - dt*x(1) - dt*s)/dt2_4;
            dx2 = dt*(2*p(2) - 2*x(1) - dt*p(1) - dt*x(2) - 2*s)/dt2_4;
            dx3 = hdt*p(3);

            [x(1), e_x(1)] = comp_sum(x(1),e_x(1),dx1);
            [x(2), e_x(2)] = comp_sum(x(2),e_x(2),dx2);
            [x(3), e_x(3)] = comp_sum(x(3),e_x(3),dx3);

            %% Force at x_n2
            dU = obj.partialU(x);

            %% Stage 2, p <- p_n1
            dp1 = dt*(-2*dt*p(1) + 4*p(2) - 4*dU(1) - 2*dt*dU(2))/dt2_4;
            dp2 = dt*(-4*p(1) - 2*dt*p(2) + 2*dt*dU(1) - 4*dU(2))/dt2_4;
            dp3 = -dt*dU(3);

            [p(1), e_p(1)] = comp_sum(p(1),e_p(1),dp1);
            [p(2), e_p(2)] = comp_sum(p(2),e_p(2),dp2);
            [p(3), e_p(3)] = comp_sum(p(3),e_p(3),dp3);

            %% Stage 3, x <- x_n1
            dxx1 = hdt*(x(2) + p(1));
            dxx2 = hdt*(p(2) - x(1) - s);
            dxx3 = hdt*p(3);

            [x(1), e_x(1)] = comp_sum(x(1),e_x(1),dxx1);
            [x(2), e_x(2)] = comp_sum(x(2),e_x(2),dxx2);
            [x(3), e_x(3)] = comp_sum(x(3),e_x(3),dxx3);
        end

        function [x_n1,p_n1,e_x2,e_x,e_p,e_dt2_4] = SI_EOM_CS(obj, dt, phi_l, x, p, e_x2, e_x, e_p, e_dt2_4)
            % Closed-form update, Kahan-compensated on every partial sum
            % (24 adds/step). Center-agnostic like SI_EOM_ICS.
            % e_x2 is 9 slots, e_x is 6.
            dt = phi_l*dt;
            hdt = dt/2;
            s = 1 - obj.mu - obj.r2;

            %% dt^2 + 4
            [dt2_4, e_dt2_4] = comp_sum(4,e_dt2_4,dt^2);

            %% x_n2(1) = (4*x(1) + 2*dt*p(1) + 2*dt*x(2) + dt^2*p(2) - dt^2*s) / (dt^2+4)
            [s1, e_x2(1)] = comp_sum(4*x(1),e_x2(1),2*dt*p(1));
            [s2, e_x2(2)] = comp_sum(s1,e_x2(2),2*dt*x(2));
            [s3, e_x2(3)] = comp_sum(s2,e_x2(3),dt^2*p(2));
            [s4, e_x2(4)] = comp_sum(s3,e_x2(4),-dt^2*s);
            x_n2_1 = s4/dt2_4;

            %% x_n2(2) = (-2*dt*x(1) - dt^2*p(1) + 4*x(2) + 2*dt*p(2) - 2*dt*s) / (dt^2+4)
            [s1, e_x2(5)] = comp_sum(4*x(2),e_x2(5),2*dt*p(2));
            [s2, e_x2(6)] = comp_sum(s1,e_x2(6),-2*dt*x(1));
            [s3, e_x2(7)] = comp_sum(s2,e_x2(7),-dt^2*p(1));
            [s4, e_x2(8)] = comp_sum(s3,e_x2(8),-2*dt*s);
            x_n2_2 = s4/dt2_4;

            %% x_n2(3) = x(3) + hdt*p(3)
            [x_n2_3, e_x2(9)] = comp_sum(x(3),e_x2(9),hdt*p(3));

            x_n2 = [x_n2_1; x_n2_2; x_n2_3];

            %% Force at x_n2
            dU = obj.partialU(x_n2);

            %% p_n1(1) = -(dt^2*p(1) - 4*p(1) - 4*dt*p(2) + 4*dt*dU(1) + 2*dt^2*dU(2)) / (dt^2+4)
            [s1, e_p(1)] = comp_sum(dt^2*p(1),e_p(1),-4*p(1));
            [s2, e_p(2)] = comp_sum(s1,e_p(2),-4*dt*p(2));
            [s3, e_p(3)] = comp_sum(s2,e_p(3),4*dt*dU(1));
            [s4, e_p(4)] = comp_sum(s3,e_p(4),2*dt^2*dU(2));
            p_n1_1 = -s4/dt2_4;

            %% p_n1(2) = -(4*dt*p(1) - 4*p(2) + dt^2*p(2) - 2*dt^2*dU(1) + 4*dt*dU(2)) / (dt^2+4)
            [s1, e_p(5)] = comp_sum(4*dt*p(1),e_p(5),-4*p(2));
            [s2, e_p(6)] = comp_sum(s1,e_p(6),dt^2*p(2));
            [s3, e_p(7)] = comp_sum(s2,e_p(7),-2*dt^2*dU(1));
            [s4, e_p(8)] = comp_sum(s3,e_p(8),4*dt*dU(2));
            p_n1_2 = -s4/dt2_4;

            %% p_n1(3) = p(3) - dt*dU(3)
            [p_n1_3, e_p(9)] = comp_sum(p(3),e_p(9),-dt*dU(3));

            %% x_n1(1) = x_n2(1) + hdt*x_n2(2) + hdt*p_n1(1)
            [s1, e_x(1)] = comp_sum(x_n2_1,e_x(1),hdt*x_n2_2);
            [x_n1_1, e_x(2)] = comp_sum(s1,e_x(2),hdt*p_n1_1);

            %% x_n1(2) = x_n2(2) - hdt*x_n2(1) + hdt*p_n1(2) - hdt*s
            [s1, e_x(3)] = comp_sum(x_n2_2,e_x(3),-hdt*x_n2_1);
            [s2, e_x(4)] = comp_sum(s1,e_x(4),hdt*p_n1_2);
            [x_n1_2, e_x(5)] = comp_sum(s2,e_x(5),-hdt*s);

            %% x_n1(3) = x_n2(3) + hdt*p_n1(3)
            [x_n1_3, e_x(6)] = comp_sum(x_n2_3,e_x(6),hdt*p_n1_3);

            x_n1 = [x_n1_1; x_n1_2; x_n1_3];
            p_n1 = [p_n1_1; p_n1_2; p_n1_3];
        end

        % Double-double (dd) precision kernels, state-only (no STM).
        % Arithmetic toolkit lives in lib/utils/.

        function [xh,xl,ph,pl] = SI_EOM_dd(obj, dt, phi_l, xh, xl, ph, pl)
            % Scheme 2 (Stormer-Verlet B) with T/D linear algebra in dd;
            % force still evaluated at double precision.
            dt  = phi_l*dt;
            hdt = dt/2;

            %% den = 1 + hdt^2, inv_den = 1/den, hdt/den
            hdt2   = dd_mul([hdt,0], [hdt,0]);
            den    = dd_add([1,0], hdt2);
            invden = dd_recip(den);
            hdtden = dd_mul([hdt,0], invden);

            %% v = x + hdt*p
            v1 = dd_add([xh(1),xl(1)], dd_mul([hdt,0],[ph(1),pl(1)]));
            v2 = dd_add([xh(2),xl(2)], dd_mul([hdt,0],[ph(2),pl(2)]));
            v3 = dd_add([xh(3),xl(3)], dd_mul([hdt,0],[ph(3),pl(3)]));

            %% x_n2 = Tinv * v
            xn2_1 = dd_add(dd_mul(invden, v1), dd_mul(hdtden, v2));
            xn2_2 = dd_add(dd_mul(dd_neg(hdtden), v1), dd_mul(invden, v2));
            xn2_3 = v3;

            %% Force at x_n2, evaluated at double precision (high word only)
            x_n2 = [xn2_1(1); xn2_2(1); xn2_3(1)];
            dU = obj.partialU(x_n2);

            %% w = D*p - dt*dU
            w1 = dd_add([ph(1),pl(1)], dd_mul([hdt,0],[ph(2),pl(2)]));
            w2 = dd_add(dd_mul([-hdt,0],[ph(1),pl(1)]), [ph(2),pl(2)]);
            w3 = [ph(3),pl(3)];

            w1 = dd_add(w1, dd_mul([-dt,0],[dU(1),0]));
            w2 = dd_add(w2, dd_mul([-dt,0],[dU(2),0]));
            w3 = dd_add(w3, dd_mul([-dt,0],[dU(3),0]));

            %% p_n1 = Tinv * w
            pn1_1 = dd_add(dd_mul(invden, w1), dd_mul(hdtden, w2));
            pn1_2 = dd_add(dd_mul(dd_neg(hdtden), w1), dd_mul(invden, w2));
            pn1_3 = w3;

            %% x_n1 = D*x_n2 + hdt*p_n1
            g1 = dd_add(xn2_1, dd_mul([hdt,0], xn2_2));
            g2 = dd_add(dd_mul([-hdt,0], xn2_1), xn2_2);

            xn1_1 = dd_add(g1, dd_mul([hdt,0], pn1_1));
            xn1_2 = dd_add(g2, dd_mul([hdt,0], pn1_2));
            xn1_3 = dd_add(xn2_3, dd_mul([hdt,0], pn1_3));

            %% Output
            xh = [xn1_1(1); xn1_2(1); xn1_3(1)];
            xl = [xn1_1(2); xn1_2(2); xn1_3(2)];
            ph = [pn1_1(1); pn1_2(1); pn1_3(1)];
            pl = [pn1_1(2); pn1_2(2); pn1_3(2)];
        end

        function [xh,xl,ph,pl] = SI_EOM_ddInc(obj, dt, phi_l, xh, xl, ph, pl)
            % dd analogue of SI_EOM_ICS: same increments as
            % SI_EOM_Increment, folded in via dd_accum instead of Kahan.
            dt    = phi_l*dt;
            hdt   = dt/2;
            dt2_4 = 4 + dt^2;

            x = xh + xl;
            p = ph + pl;

            %% Stage 1, x <- x_n2
            dx1 = dt*(2*p(1) + 2*x(2) + dt*p(2) - dt*x(1))/dt2_4;
            dx2 = dt*(2*p(2) - 2*x(1) - dt*p(1) - dt*x(2))/dt2_4;
            dx3 = hdt*p(3);

            [xh(1),xl(1)] = dd_accum(xh(1),xl(1),dx1);
            [xh(2),xl(2)] = dd_accum(xh(2),xl(2),dx2);
            [xh(3),xl(3)] = dd_accum(xh(3),xl(3),dx3);

            x = xh + xl;

            %% Force at x_n2
            dU = obj.partialU(x);

            %% Stage 2, p <- p_n1
            dp1 = dt*(-2*dt*p(1) + 4*p(2) - 4*dU(1) - 2*dt*dU(2))/dt2_4;
            dp2 = dt*(-4*p(1) - 2*dt*p(2) + 2*dt*dU(1) - 4*dU(2))/dt2_4;
            dp3 = -dt*dU(3);

            [ph(1),pl(1)] = dd_accum(ph(1),pl(1),dp1);
            [ph(2),pl(2)] = dd_accum(ph(2),pl(2),dp2);
            [ph(3),pl(3)] = dd_accum(ph(3),pl(3),dp3);

            p = ph + pl;

            %% Stage 3, x <- x_n1
            dxx1 = hdt*(x(2) + p(1));
            dxx2 = hdt*(p(2) - x(1));
            dxx3 = hdt*p(3);

            [xh(1),xl(1)] = dd_accum(xh(1),xl(1),dxx1);
            [xh(2),xl(2)] = dd_accum(xh(2),xl(2),dxx2);
            [xh(3),xl(3)] = dd_accum(xh(3),xl(3),dxx3);
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

end
