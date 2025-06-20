%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Author : Anabel Soria-Carro
% Date   : June 19, 2025
% Affiliation: The University of Texas at Austin
%              Controls Group for Distributed and Uncertain Systems (CDUS)
% Description:
%  This class defines the Elliptic Restricted Three-Body Problem (ER3BP)
%  dynamical system.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

classdef ER3BP < astro.DynamicalSystem
    
    % These properties will be set by the user
    properties 
        center  string % Origin coordinates
    end

    % These properties are internally set by the class
    properties
        r1      % Position of larger primary
        r2      % Position of secondary primary
        e       % Eccentricity

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
        function obj = ER3BP(mu,center,e)
            obj@astro.DynamicalSystem();
            obj.mu = mu;
            obj.center = center;
            obj.e = e;
        end

        function ds = EOM(obj, t, s)
            ds = NaN;
        end

        function dh = Hamiltons_EOM(obj, t, s)
            q = s(1:4);
            p = s(5:end);

            f = 1/(1+obj.e*cos(p(4)));

            dq = [p(1) + q(2);...
                  p(2) - q(1);...
                  p(3)       ;...
                  -0.5 * f^2 * obj.e * sin(p(4)) * ( norm(q(1:3))^2 - 2*obj.U(q) )];

            dU = obj.partialU(s);

            dp = -[-p(2) + f*( obj.e*cos(p(4))* q(1) + dU(1));...
                    p(1) + f*( obj.e*cos(p(4))* q(2) + dU(2));...
                    f*( obj.e*cos(p(4))* q(3) + dU(3));...
                   -1];

            dh = [dq;dp];
        end

        function U = U(obj,q)

            if isrow(q)
                q = reshape(q, [], 1); % Ensure q is a 2D matrix
            end

            % Extract components
            x = q(1,:);
            y = q(2,:);
            z = q(3,:);

            % Mass parameters
            mu1 = 1 - obj.mu; % Mass of larger primary
            mu2 = obj.mu;     % Mass of smaller primary

            % Distances (squared)
            d = (x - obj.r1).^2 + y.^2 + z.^2; % Distance to larger primary
            r = (x - obj.r2).^2 + y.^2 + z.^2; % Distance to smaller primary

            % Compute U
            U = -mu1 ./ sqrt(d) - obj.mu ./ sqrt(r);
        end

        function dUdx = partialU(obj, q) % Same as CR3BP
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
                    % p_n2 = T * (p - dt2 * obj.partialU(q));
                    % q_n1 = T * (D * q + dt * p_n2);
                    % p_n1 = D * p_n2 - dt2 * obj.partialU(q_n1);
                    error('Scheme 1 not implemented yet')

                case 2 % Stormer-Verlet B
                    q_n2    = T * (q(1:3,1) + dt2 * (p(1:3,1) - [0; 1-obj.mu - obj.r2;0]) );
                    q_n2(4) = q(4) - dt*obj.e*sin(p(4))/(4*(1+obj.e*cos(p(4)))^2)*...
                              (norm(q_n2(1:3))^2 - 2*obj.U(q_n2));

                    p_n1 = zeros(4,1);
                    p_n1(4) = p(4) + dt;
                    p_n1(1:3) = T*( D*p(1:3,1) ...
                        - dt2 * 1/(1+obj.e*cos(p(4)))   * (obj.e*cos(p(4)) *q_n2(1:3,1) + obj.partialU(q_n2)) +...
                        - dt2 * 1/(1+obj.e*cos(p_n1(4)))* (obj.e*cos(p_n1(4))*q_n2(1:3,1) + obj.partialU(q_n2)));


                    q_n1 = zeros(4,1);
                    q_n1(1:3) = D*q_n2(1:3,1) + dt2 * (p_n1(1:3,1) - [0; 1-obj.mu - obj.r2;0]);
                    q_n1(4)   = q_n2(4) - dt/2 * obj.e/2*sin(p_n1(4))/(1+obj.e*cos(p_n1(4)))^2 * (norm(q_n2(1:3))^2 - 2*obj.U(q_n2));
            end

        end

        function K = kamiltonian(obj,X)

            q = X(1:4,:);
            p = X(5:end,:);
            
            f = p(end,:);

            K = 0.5 * vecnorm(p(1:3,:),2,1).^2 + p(1,:).*q(2,:) - p(2,:).*q(1,:) + ...
            1./(1+obj.e.*cos(f)) .* ( obj.e/2 .*cos(f).*vecnorm(q(1:3,:),2,1).^2 + obj.U(q) ) - q(4,:);

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
