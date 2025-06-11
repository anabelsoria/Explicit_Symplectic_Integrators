%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Author : Anabel Soria-Carro
% Date   : May 22, 2025
% Affiliation: The University of Texas at Austin
%              Controls Group for Distributed and Uncertain Systems (CDUS)
% Description:
%  This class defines the Circular Restricted Three-Body Problem (CR3BP)
%  dynamical system.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%   schemes for propagation.
classdef CR3BP < astro.DynamicalSystem
    properties
        mu % Mass parameter (μ = m2 / (m1 + m2))
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
        function obj = CR3BP(mu)
            obj@astro.DynamicalSystem();
            obj.mu = mu;
        end

        function ds = EOM(obj, t, s)
            n = length(s);
            planar = (n <= 4);

            mu1 = 1 - obj.mu;
            mu2 = obj.mu;

            if planar
                x = s(1); y = s(2);
                xDot = s(3); yDot = s(4);

                r1 = ((x + mu2)^2 + y^2)^1.5;
                r2 = ((x - mu1)^2 + y^2)^1.5;

                Ux = x - mu1 * (x + mu2) / r1 - mu2 * (x - mu1) / r2;
                Uy = y - mu1 * y / r1 - mu2 * y / r2;

                xDDot = 2 * yDot + Ux;
                yDDot = -2 * xDot + Uy;

                ds = [xDot; yDot; xDDot; yDDot];
            else
                x = s(1); y = s(2); z = s(3);
                xDot = s(4); yDot = s(5); zDot = s(6);

                r1 = (x + mu2)^2 + y^2 + z^2;
                r2 = (x - mu1)^2 + y^2 + z^2;

                r13 = r1^1.5;
                r23 = r2^1.5;

                Ux = x - mu1 * (x + mu2) / r13 - mu2 * (x - mu1) / r23;
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

            x1 = -obj.mu; x2 = 1 - obj.mu;
            GM1 = 1 - obj.mu;
            GM2 = obj.mu;

            R1 = sqrt((x - x1)^2 + y^2 + z^2);
            R2 = sqrt((x - x2)^2 + y^2 + z^2);

            Ux = GM1 * (x - x1) / R1^3 + GM2 * (x - x2) / R2^3;
            Uy = GM1 * y / R1^3       + GM2 * y / R2^3;
            Uz = GM1 * z / R1^3       + GM2 * z / R2^3;

            dh = zeros(6, 1);
            dh(1) = px + y;
            dh(2) = py - x;
            dh(3) = pz;
            dh(4) = -( -py + Ux );
            dh(5) = -(  px + Uy );
            dh(6) = -Uz;
        end

        function dUdx = partialU(obj, q)
            x1 = -obj.mu; x2 = 1 - obj.mu;
            GM1 = 1 - obj.mu; GM2 = obj.mu;

            x = q(1); y = q(2);
            z = 0;
            if length(q) > 2
                z = q(3);
            end

            R1 = sqrt((x - x1)^2 + y^2 + z^2);
            R2 = sqrt((x - x2)^2 + y^2 + z^2);

            Ux = GM1 * (x - x1) / R1^3 + GM2 * (x - x2) / R2^3;
            Uy = GM1 * y / R1^3 + GM2 * y / R2^3;
            Uz = GM1 * z / R1^3 + GM2 * z / R2^3;

            dUdx = [Ux; Uy];
            if length(q) > 2
                dUdx = [dUdx; Uz];
            end
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
                    q_n2 = T * (q + dt2 * p);
                    p_n1 = T * (D * p - dt * obj.partialU(q_n2));
                    q_n1 = D * q_n2 + dt2 * p_n1;
            end

        end

        function C = jacobiconstant(obj, s)
            if size(s,1) > 4
                x = s(1,:); y = s(2,:); z = s(3,:);
                v = s(4:end,:);
            else
                x = s(1,:); y = s(2,:); z = 0;
                v = s(3:end,:);
            end

            r1 = sqrt((x + obj.mu).^2 + y.^2 + z.^2);
            r2 = sqrt((x - 1 + obj.mu).^2 + y.^2 + z.^2);
            U  = (1 - obj.mu) ./ r1 + obj.mu ./ r2;

            C = x.^2 + y.^2 + 2 * U - vecnorm(v).^2;
        end

        function Hp = Hp(obj,s)
            q1 = s(1);
            q2 = s(2);
            q3 = s(3);
            p1 = s(4);
            p2 = s(5);
            p3 = s(6);

            t2 = p1+q2;
            t3 = -p2;
            t4 = obj.mu+q1+t3-1.0;
            Hp = [abs(t2).*sign(t2),-abs(t4).*sign(t4),abs(p3).*sign(p3)];
        end
    end
end
