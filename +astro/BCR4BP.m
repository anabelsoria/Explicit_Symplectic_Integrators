%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Author : Anabel Soria-Carro
% Date   : June 19, 2025
% Affiliation: The University of Texas at Austin
%              Controls Group for Distributed and Uncertain Systems (CDUS)
% Description:
%  This class defines the Elliptic Restricted Three-Body Problem (ER3BP)
%  dynamical system.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

classdef BCR4BP < astro.DynamicalSystem

    % These properties will be set by the user
    properties
        center  string % Origin coordinates
    end

    % These properties are internally set by the class
    properties
        r1        % Position of larger primary
        r2        % Position of secondary primary
        theta_dot %
        theta_0   % Initial Sun-B1 angle
        mu_em     % mu Earh-Moon
        mu_s      % mu Sun-B1
        a_s       % radius of the Sun wrt B1
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
        function obj = BCR4BP(mu_em,mu_s,center,LU,m_sun,theta_0)

            import astro.Constants;

            obj@astro.DynamicalSystem();
            obj.mu = mu_em;
            obj.mu_em = mu_em;
            obj.mu_s  = mu_s;
            obj.center = center;
            obj.a_s = (Constants.AU/LU) - obj.mu_em;

            obj.theta_0   = theta_0;
            obj.theta_dot = sqrt((1+obj.mu_s)/obj.a_s^3) - 1;

            obj.LU = LU; % [km]
            obj.TU = sqrt(LU^3/(Constants.G*m_sun)); % [s]

            obj.theta_dot = sqrt((1+obj.mu_s)/obj.a_s^3) - 1;
        end

        function ds = EOM(obj, t, s)
            % Unpack the state variables
            x = s(1);
            y = s(2);
            z = s(3);

            % Find P4 position
            theta = s(end);

            x4 = obj.a_s*cos(theta);
            y4 = obj.a_s*sin(theta);
            z4 = 0;

            mu1 = 1 - obj.mu_em;
            mu2 = obj.mu_em;

            % The Distances from the larger and Smaller Primary
            r13 = (x - obj.r1)^2 + y^2 + z^2;      % r13: distance to m1, LARGER MASS
            r23 = (x - obj.r2)^2 + y^2 + z^2;      % r23: distance to m2, smaller mass
            r43 = (x - x4)^2 + (y - y4)^2 + (z - z4)^2;

            % Partial Derivative of the Pseudo Potential Function
            psi_x =  x - mu1*(x - obj.r1)/r13^(3/2) - mu2*(x - obj.r2)/r23^(3/2)...
                - ( (obj.mu_s*x4)/obj.a_s^3 + obj.mu_s*(x - x4)/r43^(3/2) );

            psi_y =  y - (mu1*y)/r13^(3/2)  - (mu2*y)/r23^(3/2)...
                - ((obj.mu_s*y4)/obj.a_s^3 + (obj.mu_s*(y - y4))/r43^(3/2));

            psi_z = - mu1*z/r13^(3/2) - (mu2*z)/r23^(3/2)...
                - ((obj.mu_s*z4)/obj.a_s^3 + (obj.mu_s*(z - z4))/r43^(3/2));

            xDot  = s(4);
            yDot  = s(5);
            zDot  = s(6);
            xDDot = 2*yDot + psi_x;
            yDDot = -2*xDot + psi_y;
            zDDot = psi_z;
            thetaDot = obj.theta_dot;

            ds = [xDot;yDot;zDot;xDDot;yDDot;zDDot;thetaDot];
        end

        function dh = Hamiltons_EOM(obj, t, s)

            dq =   obj.Hp(s);
            dp = - obj.Hq(s);

            dh = [dq;dp];
        end

        function Hp = Hp(obj,s)

            q = s(1:4); p = s(5:end);

            theta = obj.theta_0 + obj.theta_dot*p(4);

            xs = obj.a_s*cos(theta);
            ys = obj.a_s*sin(theta);
            zs = 0;
            r43 = sqrt( (q(1) - xs)^2 + (q(2) - ys)^2 + (q(3) - zs)^2) ;

            Hp = [p(1) + q(2);...
                p(2) - q(1);...
                p(3)       ;...
                (obj.mu_s*obj.theta_dot*(q(2)*cos(theta) - q(1)*sin(theta)))/obj.a_s^2 ...
                - (obj.a_s*obj.mu_s*obj.theta_dot*(q(2)*cos(theta) ...
                - q(1)*sin(theta)))/r43^3];
        end

        function Hq = Hq(obj,s)

            q = s(1:4); p = s(5:end); 

            % Partial Derivative of the Pseudo Potential Function (not psi_star)
            psi = obj.partialU(p(4),q);

            Hq = [-p(2) - psi(1);...
                   p(1) - psi(2);...
                        - psi(3);...
                  -1];
        end

        function dU = partialU(obj,p4,q)

            x = q(1); y = q(2);
            z = 0;
            if length(q) > 2
                z = q(3);
            end

            % Find P4 position
            theta = obj.theta_0 + obj.theta_dot*p4;

            x4 = obj.a_s*cos(theta);
            y4 = obj.a_s*sin(theta);
            z4 = 0;

            mu1 = 1 - obj.mu_em;
            mu2 = obj.mu_em;

            % The Distances from the larger and Smaller Primary
            r13 = (x - obj.r1)^2 + y^2 + z^2;      % r13: distance to m1, LARGER MASS
            r23 = (x - obj.r2)^2 + y^2 + z^2;      % r23: distance to m2, smaller mass
            r43 = (x - x4)^2 + (y - y4)^2 + (z - z4)^2;

            % Partial Derivative of the Pseudo Potential Function
            psi_x = - mu1*(q(1) - obj.r1)/r13^(3/2) - mu2*(q(1) - obj.r2)/r23^(3/2)...
                - ( (obj.mu_s*x4)/obj.a_s^3 + obj.mu_s*(q(1) - x4)/r43^(3/2) );

            psi_y = - mu1*q(2)/r13^(3/2)  - (mu2*q(2))/r23^(3/2)...
                - ((obj.mu_s*y4)/obj.a_s^3 + (obj.mu_s*(q(2) - y4))/r43^(3/2));

            psi_z = - mu1*q(3)/r13^(3/2) - (mu2*q(3))/r23^(3/2)...
                - ((obj.mu_s*z4)/obj.a_s^3 + (obj.mu_s*(q(3) - z4))/r43^(3/2));

            dU = [psi_x; psi_y; psi_z];

        end

        function [q_n1, p_n1] = SI_EOM(obj, dt, scheme, X)
            q = X(1:4);
            p = X(5:end);

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
                    error('Scheme 1 not implemented yet')

                case 2 % Stormer-Verlet B
                    q_n2      = zeros(4,1);
                    q_n2(1:3) = T * (q(1:3) + dt2 * (p(1:3) - [0; 1-obj.mu - obj.r2;0]));
                    Hp_qn2 = obj.Hp([q_n2;p]);
                    q_n2(4)   = q(4) + dt/2 * Hp_qn2(end);

                    p_n1    = zeros(4,1);
                    p_n1(4) = p(4) + dt;

                    dU_n  = obj.partialU(p(4),q_n2(1:3)+ [1-obj.mu - obj.r2;0;0]);
                    dU_n1 = obj.partialU(p_n1(4),q_n2(1:3)+ [1-obj.mu - obj.r2;0;0]);

                    p_n1(1:3) = T*( D*p(1:3,1) - dt/2 *  (-dU_n - dU_n1) );

                    q_n1 = q_n2 + dt/2*obj.Hp([q_n2; p_n1-[0;1-obj.mu - obj.r2;0;0]]);
            end

        end

        function H = hamiltonian(obj,s)

            q = s(1:4,:); p = s(5:end,:);

            qx = q(1,:); qy = q(2,:); qz = q(3,:);
            qt = q(4,:);  % Unused but extracted for completeness
            px = p(1,:); py = p(2,:); pz = p(3,:);
            pt = p(4,:);  

            mu1 = 1 - obj.mu_em;
            mu2 = obj.mu_em;

            t = pt;

            theta = obj.theta_0 + obj.theta_dot * t;

            x4 = obj.a_s * cos(theta);
            y4 = obj.a_s * sin(theta);
            z4 = zeros(size(theta));

            % Distances
            r13 = sqrt((qx - obj.r1).^2 + qy.^2 + qz.^2);
            r23 = sqrt((qx - obj.r2).^2 + qy.^2 + qz.^2);
            r43 = sqrt((qx - x4).^2 + (qy - y4).^2 + (qz - z4).^2);

            % Psi potential
            psi = 0.5*(qx.^2 + qy.^2) + mu1./r13 + mu2./r23 + ...
                (obj.mu_s ./ r43 ...
                - obj.mu_s ./ obj.a_s.^3 .* (x4 .* qx + y4 .* qy + z4 .* qz));

            % Hamiltonian
            % H = 2*psi - ((px + qy).^2 + (py - qx).^2 + pz.^2); % Stephen
            H = 1/2 * ((px + qy).^2 + (py - qx).^2 + pz.^2) - psi - qt; 
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
