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
        function obj = CR3BP(mu,center)
            obj@astro.DynamicalSystem();
            obj.mu = mu;
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

            Ux = mu1 * (x - obj.r1) / r13^3 + mu2 * (x - obj.r2) / r23^3;
            Uy = mu1 * y / r13^3       + mu2 * y / r23^3;
            Uz = mu1 * z / r13^3       + mu2 * z / r23^3;

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

            r13 = sqrt((x - obj.r1).^2 + y.^2 + z.^2);
            r23 = sqrt((x - obj.r2).^2 + y.^2 + z.^2);
            U  = (1 - obj.mu) ./ r13 + obj.mu ./ r23;

            C = (x + (1 - obj.mu - obj.r2)).^2 + y.^2 + 2 * U - vecnorm(v).^2;
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

    end
end
