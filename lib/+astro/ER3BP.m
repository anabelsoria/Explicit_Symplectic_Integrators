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
        function obj = ER3BP(mu,center,e,LU,TU)
            obj@astro.DynamicalSystem();
            obj.mu = mu;
            obj.LU = LU;
            obj.TU = TU;
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

            dq = obj.Hp(s);

            dU = obj.partialU(s);

            dp = -[-p(2) + f*( obj.e*cos(p(4))* (q(1)+1-obj.mu - obj.r2) + dU(1));...
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

        function Hp = Hp(obj,s)
            q = s(1:4); p = s(5:end);

            f = 1/(1+obj.e*cos(p(4)));
            Hp = [p(1) + q(2);...
                  p(2) - (q(1)+ 1-obj.mu - obj.r2);... 
                  p(3);...
                  -0.5 * f^2 * obj.e * sin(p(4)) * ( norm(q(1:3)+[1-obj.mu - obj.r2;0;0])^2 - 2*obj.U(q) )];
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
                              (norm(q_n2(1:3)+[1-obj.mu - obj.r2;0;0])^2 - 2*obj.U(q_n2));

                    p_n1 = zeros(4,1);
                    p_n1(4) = p(4) + dt;
                    p_n1(1:3) = T*( D* p(1:3,1) ...
                        - dt2 * 1/(1+obj.e*cos(p(4)))   * (obj.e*cos(p(4)) * (q_n2(1:3,1)+[1-obj.mu - obj.r2;0;0]) +...
                        obj.partialU(q_n2)) +...
                        - dt2 * 1/(1+obj.e*cos(p_n1(4)))* (obj.e*cos(p_n1(4))*(q_n2(1:3,1)+[1-obj.mu - obj.r2;0;0]) +...
                        obj.partialU(q_n2)));


                    q_n1 = zeros(4,1);
                    q_n1(1:3) = D*q_n2(1:3,1) + dt2 * (p_n1(1:3,1) - [0; 1-obj.mu - obj.r2;0]);
                    q_n1(4)   = q_n2(4) - dt/2 * obj.e/2*sin(p_n1(4))/(1+obj.e*cos(p_n1(4)))^2 * ...
                    (norm(q_n2(1:3)+[1-obj.mu - obj.r2;0;0])^2 - 2*obj.U(q_n2));
            end

        end


        % Precision-ladder kernels. Bary-centered (no p2 shift), unlike
        % SI_EOM above which assumes center='p2'.

        function [q,p,e_q,e_p] = SI_EOM_ICS(obj, dt, phi_l, q, p, e_q, e_p)
            % Increment + Kahan-compensated summation. Spatial increments
            % match the CR3BP update; force term h*dU is replaced by G,
            % evaluated before/after the true anomaly advances.
            mu = obj.mu;
            ec = obj.e;

            h   = phi_l*dt;
            hd  = h/2;
            d4  = h*h + 4;

            %% Stage 1, q(1:3) advances to q_n2. Increments formed before any update.
            dq1 = h*( 2*p(1) + 2*q(2) + h*p(2) - h*q(1) )/d4;
            dq2 = h*( 2*p(2) - 2*q(1) - h*p(1) - h*q(2) )/d4;
            dq3 = hd*p(3);

            [q(1), e_q(1)] = comp_sum(q(1),e_q(1),dq1);
            [q(2), e_q(2)] = comp_sum(q(2),e_q(2),dq2);
            [q(3), e_q(3)] = comp_sum(q(3),e_q(3),dq3);

            % Values at q_n2, needed again in stage 3 before q(1:3) advances further.
            n1 = q(1);  n2 = q(2);  n3 = q(3);
            qn2 = [n1;n2;n3];

            W  = n1*n1 + n2*n2 + n3*n3 - 2*obj.U(qn2);
            dU = obj.partialU(qn2);

            %% Stage 1b, the coordinate conjugate to the true anomaly
            nu  = p(4);
            c   = cos(nu);
            dq4 = -h*ec*sin(nu)/(4*(1 + ec*c)^2) * W;
            [q(4), e_q(4)] = comp_sum(q(4),e_q(4),dq4);

            %% Stage 2, the true anomaly advances, then the momenta
            [p(4), e_p(4)] = comp_sum(p(4),e_p(4),h);
            nu2 = p(4);
            c2  = cos(nu2);

            G = hd*( (ec*c *qn2 + dU)/(1 + ec*c ) ...
                   + (ec*c2*qn2 + dU)/(1 + ec*c2) );

            den = 1 + hd*hd;
            dp1 = ( -h*hd*p(1) + h*p(2) - G(1) - hd*G(2) )/den;
            dp2 = ( -h*p(1) - h*hd*p(2) + hd*G(1) - G(2) )/den;
            dp3 = -G(3);

            [p(1), e_p(1)] = comp_sum(p(1),e_p(1),dp1);
            [p(2), e_p(2)] = comp_sum(p(2),e_p(2),dp2);
            [p(3), e_p(3)] = comp_sum(p(3),e_p(3),dp3);

            %% Stage 3, q(1:3) advances to q_n1 using q_n2 and the updated momenta
            dr1 = hd*( n2 + p(1) );
            dr2 = hd*( p(2) - n1 );
            dr3 = hd*p(3);

            [q(1), e_q(1)] = comp_sum(q(1),e_q(1),dr1);
            [q(2), e_q(2)] = comp_sum(q(2),e_q(2),dr2);
            [q(3), e_q(3)] = comp_sum(q(3),e_q(3),dr3);

            %% Stage 3b
            dq4b = -hd*(ec/2)*sin(nu2)/(1 + ec*c2)^2 * W;
            [q(4), e_q(4)] = comp_sum(q(4),e_q(4),dq4b);
        end

        function [q_n1,p_n1,e_q2,e_q,e_dt] = SI_EOM_CS(obj, dt, phi_l, q, p, e_q2, e_q, e_dt)
            % Compensated summation on the spatial q update only.
            mu = obj.mu;
            e  = obj.e;

            dt = phi_l*dt;

            [den, e_dt] = comp_sum(1,e_dt,(dt/2)^2);

            T = 1/den * [1         dt/2   0;...
                         -dt/2       1      0;...
                         0           0      den];

            D = [1         dt/2    0;...
                -dt/2      1       0;...
                0            0       1];

            % q_n2 = T * (q(1:3) + dt/2 * p(1:3));
            delta_q2 = dt/2 * p(1:3,1);
            [q_n2, e_q2] = comp_sum(T*q(1:3,1),e_q2,T*delta_q2);
            q_n2(4) = q(4) - dt/2 * e/2*sin(p(4))/(1+e*cos(p(4)))^2 * ( norm(q_n2(1:3))^2 - 2*obj.U(q_n2));

            p_n1 = zeros(4,1);
            p_n1(4) = p(4) + dt;
            p_n1(1:3) = T*( D*p(1:3,1) ...
                - dt/2 * 1/(1+e*cos(p(4)))   * (e*cos(p(4)) *q_n2(1:3,1) + obj.partialU(q_n2)) +...
                - dt/2 * 1/(1+e*cos(p_n1(4)))* (e*cos(p_n1(4))*q_n2(1:3,1) + obj.partialU(q_n2)));

            q_n1 = zeros(4,1);
            % q_n1(1:3) = D*q_n2(1:3) + dt/2 * p_n1(1:3);
            delta_q = dt/2 * p_n1(1:3,1);
            [q_n1(1:3), e_q] = comp_sum(D*q_n2(1:3,1),e_q,delta_q);

            q_n1(4)   = q_n2(4) - dt/2 * e/2*sin(p_n1(4))/(1+e*cos(p_n1(4)))^2 * (norm(q_n2(1:3))^2 - 2*obj.U(q_n2));
        end

        function K = kamiltonian(obj,sol)

            q = sol.x(1:4,:);
            p = sol.x(5:end,:);

            f = p(end,:);

            K = 0.5 * vecnorm(p(1:3,:),2,1).^2 + p(1,:).*q(2,:) - p(2,:).* (q(1,:) + 1-obj.mu - obj.r2) + ...
            1./(1+obj.e.*cos(f)) .* ( obj.e/2 .*cos(f).*vecnorm((q(1:3,:)+[1-obj.mu - obj.r2;0;0]),2,1).^2 + ...
            obj.U(q) ) - q(4,:);

        end

        function xi = nu2xi(obj,nu)
            switch obj.center
                case 'bary'
                    xi = obj.P_nu_xi * nu;
                case 'p2'
                    xi = nu;
                    xi([1:3,5:7],:) = obj.P_nu_xi * (nu([1:3,5:7],:) + [1-obj.mu;0;0;0;0;0]);
                    xi(1,:) = xi(1,:) - (1-obj.mu);
                case 'p1'
                    xi = obj.P_nu_xi * (nu + [-obj.mu;0;0;0;0;0;0;0]);
                    xi(1,:) = xi(1,:) - (-obj.mu);
            end
        end

    end
end
