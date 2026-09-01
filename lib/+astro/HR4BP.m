%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Author : Anabel Soria-Carro
% Date   : June 19, 2025
% Affiliation: The University of Texas at Austin
%              Controls Group for Distributed and Uncertain Systems (CDUS)
% Description:
%  This class defines the Elliptic Restricted Three-Body Problem (ER3BP)
%  dynamical system.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

classdef HR4BP < astro.DynamicalSystem

    % These properties will be set by the user
    properties
        center  string % Origin coordinates
    end

    % These properties are internally set by the class
    properties
        r1        % Position of larger primary
        r2        % Position of secondary primary
        m 
    end

    properties (Constant)
        n_ord = 4;
        
        % Note: Expansion is to order 9 in m
        A_a0_sub = [1, -2/3, 7/18, -4/81, 19565/62208, -47161/93312, -2284055/3359232, -1152145/1259712, -65557603/1934917632, 4005971079828870083/13850679916489605120]; % storage index is (p + 1)
        A_anoa0 = [0, 0, 3/16, 1/2, 7/12, 11/36, -30749/110592, -1010521/829440, -18445871/6220800, -2114557853/373248000;...
            0, 0, 0, 0, 25/256, 803/1920, 6109/7200, 897599/864000, 237203647/368640000, -11098919887/14515200000;...
            0, 0, 0, 0, 0, 0, 833/12288, 27943/71680, 12275527/11289600, 27409853579/14224896000;...
            0, 0, 0, 0, 0, 0, 0, 0, 3537/65536, 18638507/48168960]; % storage index is (n,p + 1)
        A_amnoa0 = [0, 0, -19/16, -5/3, -43/36, -14/27, -7381/82944, 3574153/2488320, 55218889/9331200, 13620153029/1119744000;...
            0, 0, 0, 0, 0, 23/640, 299/2400, 56339/288000, 238200053/1105920000, 146886277/537600000;...
            0, 0, 0, 0, 0, 0, 1/192, 7477/215040, 65239/627200, 2674679587/14224896000;...
            0, 0, 0, 0, 0, 0, 0, 0, 23/6144, 795829/28901376]; % storage index is (-n,p + 1)
    end

    properties (Dependent)
        n_vals
        m2oa03
        mn
        a_0
        anoa0
        amnoa0
        P_xi_nu
        P_nu_xi
    end

    methods
        function obj = HR4BP(center,LU)

            import astro.Constants;

            obj@astro.DynamicalSystem();
            obj.center = center;

            obj.LU = LU; % [km]
            obj.TU = 29.53*24*60*60/(2*pi); % [s]
        end

        % Dependent property getters
        function val = get.P_xi_nu(obj)
            % From xi = (r,v) to nu = (q,p)
        val = [eye(3), zeros(3);
            0 -(1+obj.m) 0 1 0 0;
            (1+obj.m)  0 0 0 1 0;
            0  0         0 0 0 1];
        end

        function val = get.P_nu_xi(obj)
            [eye(3),  zeros(3);
            0    (1+obj.m) 0 1 0 0;
         -(1+obj.m)  0     0 0 1 0;
            0        0     0 0 0 1];
        end
        function val = get.n_vals(obj)
            val = (1:obj.n_ord)';
        end

        function val = get.mn(obj)
            val = transpose(obj.m.^(0:(length(obj.A_a0_sub)-1)));
        end

        function val = get.a_0(obj)
            val = (obj.m^(2/3)) * (obj.A_a0_sub * obj.mn(1:length(obj.A_a0_sub)));
        end

        function val = get.m2oa03(obj)
            if obj.m ~= 0
                val = (obj.m^2)/(obj.a_0^3);
            else
                val = 1;
            end
        end

        function val = get.anoa0(obj)
            val = obj.A_anoa0 * obj.mn;
        end

        function val = get.amnoa0(obj)
            val = obj.A_amnoa0 * obj.mn;
        end

        function ds = EOM(obj, tau, s)
            % Unpack the state variables
            x = s(1);
            y = s(2);
            z = s(3);

            dU = obj.partialU(tau,s);

            xDot  = s(4);
            yDot  = s(5);
            zDot  = s(6);
            xDDot = 2*(1+obj.m)*yDot + dU(1);
            yDDot = -2*(1+obj.m)*xDot + dU(2);
            zDDot = dU(3);

            ds = [xDot;yDot;zDot;xDDot;yDDot;zDDot];
        end

        function dh = Hamiltons_EOM(obj, t, s)

            dq =   obj.Hp(s);
            dp = - obj.Hq(s);

            dh = [dq;dp];
        end

        function Hp = Hp(obj,s)

            q = s(1:4); p = s(5:end);

            Hp = [p(1) + (1+obj.m) * q(2);...
                  p(2) - (1+obj.m) * q(1);...
                  p(3)       ;...
                  1];
        end

        function Hq = Hq(obj,s)

            q = s(1:4); 
            tau = q(end);
            p = s(5:end); 

            psi = obj.partialU(tau,s);

            Hq = [-(1+obj.m) * p(2) + (1+obj.m)^2 * q(1) - psi(1);...
                   (1+obj.m) * p(1) + (1+obj.m)^2 * q(2) - psi(2);...
                                                         - psi(3);...
                  -obj.partialU_tau(tau,s)];
        end

        function V = potential(obj,tau,s)
            x = s(1,:);
            y = s(2,:);
            z = s(3,:);

            m = obj.m;
            mu = obj.mu;

            xi_bar = sum((obj.anoa0 + obj.amnoa0).*cos(2*obj.n_vals*tau));
            eta_bar = sum((obj.anoa0 - obj.amnoa0).*sin(2*obj.n_vals*tau));

            R_a = sqrt((x + mu*(1 + xi_bar)).^2 + (y + mu*eta_bar).^2 + z.^2);
            R_b = sqrt((x - (1 - mu)*(1 + xi_bar)).^2 + (y - (1 - mu)*eta_bar).^2 + z.^2);
            
            V = 1/2*(1 + 2*m + 3/2*m^2)*(x.^2+y.^2) - 1/2*m^2*z.^2 ...
                + 3/4*m^2*((x.^2-y.^2).*cos(2*tau) - 2*x.*y.*sin(2*tau))...
                + obj.m2oa03 * ( (1-mu)./R_a + mu./R_b );

        end

        function dU = partialU(obj,tau,s)

            x = s(1);
            y = s(2);
            z = s(3);

            m = obj.m;
            mu = obj.mu;

            xi_bar = sum((obj.anoa0 + obj.amnoa0).*cos(2*obj.n_vals*tau));
            eta_bar = sum((obj.anoa0 - obj.amnoa0).*sin(2*obj.n_vals*tau));

            R_a = sqrt((x + mu*(1 + xi_bar))^2 + (y + mu*eta_bar)^2 + z^2);
            R_b = sqrt((x - (1 - mu)*(1 + xi_bar))^2 + (y - (1 - mu)*eta_bar)^2 + z^2);

            % Partial Derivative of the Pseudo Potential Function
            V_x = (1 + 2*m + (3/2)*m^2)*x + 0 + ...
                (3/4)*(m^2)*(2*x*cos(2*tau) - 2*y*sin(2*tau)) - ...
                obj.m2oa03*((1 - mu)*(x + mu*(1 + xi_bar))/R_a^3 + mu*(x - (1 - mu)*(1 + xi_bar))/R_b^3);
            V_y = (1 + 2*m + (3/2)*m^2)*y + 0 + ...
                (3/4)*(m^2)*(-2*y*cos(2*tau) - 2*x*sin(2*tau)) - ...
                obj.m2oa03*((1 - mu)*(y + mu*eta_bar)/R_a^3 + mu*(y - (1 - mu)*eta_bar)/R_b^3);
            V_z = 0 - (m^2)*z + ...
                0 - ...
                obj.m2oa03*((1 - mu)*z/R_a^3 + mu*z/R_b^3);
            
            dU = [V_x; V_y; V_z];

        end

        function dV_dtau = partialU_tau(obj,tau,s)

            x = s(1);
            y = s(2);
            z = s(3);

            m = obj.m;
            mu = obj.mu;
            
            xi_bar = sum((obj.anoa0 + obj.amnoa0).*cos(2*obj.n_vals*tau));
            eta_bar = sum((obj.anoa0 - obj.amnoa0).*sin(2*obj.n_vals*tau));

            R_a = sqrt((x + mu*(1 + xi_bar))^2 + (y + mu*eta_bar)^2 + z^2);
            R_b = sqrt((x - (1 - mu)*(1 + xi_bar))^2 + (y - (1 - mu)*eta_bar)^2 + z^2);

            dV_dRa = obj.m2oa03 * (mu - 1)/R_a^2;
            dV_dRb = -obj.m2oa03* mu/R_b^2;

            dRa_deta = (mu*(y + eta_bar*mu))/R_a;
            dRa_dxi = (mu*(x + mu*(xi_bar + 1)))/R_a;
            dRb_deta = ((y + eta_bar*(mu - 1))*(mu - 1))/R_b;
            dRb_dxi = ((mu - 1)*(x + (mu - 1)*(xi_bar + 1)))/R_b;

            deta_dtau = sum( -2*(obj.amnoa0 - obj.anoa0) .* (obj.n_vals.*cos(2*obj.n_vals*tau)) );
            dxi_dtau  = sum( -2*(obj.amnoa0 + obj.anoa0) .* (obj.n_vals.*sin(2*obj.n_vals*tau)) );

            dV_dtau = -3/2*m^2*(sin(2*tau)*(x^2 - y^2) + 2*x*y*cos(2*tau)) +...
                dV_dRa*(dRa_dxi*dxi_dtau + dRa_deta*deta_dtau) + ...
                dV_dRb*(dRb_dxi*dxi_dtau + dRb_deta*deta_dtau);

            % n_vals1 = obj.n_vals(1); n_vals2 = obj.n_vals(2);
            % n_vals3 = obj.n_vals(3); n_vals4 = obj.n_vals(4);
            % a0 = obj.a_0;
            % dV_dtau_v2 = (m^2*((mu*(2*(y - (mu - 1)*(sin(2*n_vals1*tau)*(obj.amnoa0(1) - obj.anoa0(1)) + sin(2*n_vals2*tau)*(obj.amnoa0(2) - obj.anoa0(2)) + sin(2*n_vals3*tau)*(obj.amnoa0(3) - obj.anoa0(3)) + sin(2*n_vals4*tau)*(obj.amnoa0(4) - obj.anoa0(4))))*(mu - 1)*(2*n_vals1*cos(2*n_vals1*tau)*(obj.amnoa0(1) - obj.anoa0(1)) + 2*n_vals2*cos(2*n_vals2*tau)*(obj.amnoa0(2) - obj.anoa0(2)) + 2*n_vals3*cos(2*n_vals3*tau)*(obj.amnoa0(3) - obj.anoa0(3)) + 2*n_vals4*cos(2*n_vals4*tau)*(obj.amnoa0(4) - obj.anoa0(4))) + 2*(x + (mu - 1)*(cos(2*n_vals1*tau)*(obj.amnoa0(1) + obj.anoa0(1)) + cos(2*n_vals2*tau)*(obj.amnoa0(2) + obj.anoa0(2)) + cos(2*n_vals3*tau)*(obj.amnoa0(3) + obj.anoa0(3)) + cos(2*n_vals4*tau)*(obj.amnoa0(4) + obj.anoa0(4)) + 1))*(mu - 1)*(2*n_vals1*sin(2*n_vals1*tau)*(obj.amnoa0(1) + obj.anoa0(1)) + 2*n_vals2*sin(2*n_vals2*tau)*(obj.amnoa0(2) + obj.anoa0(2)) + 2*n_vals3*sin(2*n_vals3*tau)*(obj.amnoa0(3) + obj.anoa0(3)) + 2*n_vals4*sin(2*n_vals4*tau)*(obj.amnoa0(4) + obj.anoa0(4)))))/(2*((y - (mu - 1)*(sin(2*n_vals1*tau)*(obj.amnoa0(1) - obj.anoa0(1)) + sin(2*n_vals2*tau)*(obj.amnoa0(2) - obj.anoa0(2)) + sin(2*n_vals3*tau)*(obj.amnoa0(3) - obj.anoa0(3)) + sin(2*n_vals4*tau)*(obj.amnoa0(4) - obj.anoa0(4))))^2 + (x + (mu - 1)*(cos(2*n_vals1*tau)*(obj.amnoa0(1) + obj.anoa0(1)) + cos(2*n_vals2*tau)*(obj.amnoa0(2) + obj.anoa0(2)) + cos(2*n_vals3*tau)*(obj.amnoa0(3) + obj.anoa0(3)) + cos(2*n_vals4*tau)*(obj.amnoa0(4) + obj.anoa0(4)) + 1))^2 + z^2)^(3/2)) - ((2*mu*(y - mu*(sin(2*n_vals1*tau)*(obj.amnoa0(1) - obj.anoa0(1)) + sin(2*n_vals2*tau)*(obj.amnoa0(2) - obj.anoa0(2)) + sin(2*n_vals3*tau)*(obj.amnoa0(3) - obj.anoa0(3)) + sin(2*n_vals4*tau)*(obj.amnoa0(4) - obj.anoa0(4))))*(2*n_vals1*cos(2*n_vals1*tau)*(obj.amnoa0(1) - obj.anoa0(1)) + 2*n_vals2*cos(2*n_vals2*tau)*(obj.amnoa0(2) - obj.anoa0(2)) + 2*n_vals3*cos(2*n_vals3*tau)*(obj.amnoa0(3) - obj.anoa0(3)) + 2*n_vals4*cos(2*n_vals4*tau)*(obj.amnoa0(4) - obj.anoa0(4))) + 2*mu*(x + mu*(cos(2*n_vals1*tau)*(obj.amnoa0(1) + obj.anoa0(1)) + cos(2*n_vals2*tau)*(obj.amnoa0(2) + obj.anoa0(2)) + cos(2*n_vals3*tau)*(obj.amnoa0(3) + obj.anoa0(3)) + cos(2*n_vals4*tau)*(obj.amnoa0(4) + obj.anoa0(4)) + 1))*(2*n_vals1*sin(2*n_vals1*tau)*(obj.amnoa0(1) + obj.anoa0(1)) + 2*n_vals2*sin(2*n_vals2*tau)*(obj.amnoa0(2) + obj.anoa0(2)) + 2*n_vals3*sin(2*n_vals3*tau)*(obj.amnoa0(3) + obj.anoa0(3)) + 2*n_vals4*sin(2*n_vals4*tau)*(obj.amnoa0(4) + obj.anoa0(4))))*(mu - 1))/(2*((x + mu*(cos(2*n_vals1*tau)*(obj.amnoa0(1) + obj.anoa0(1)) + cos(2*n_vals2*tau)*(obj.amnoa0(2) + obj.anoa0(2)) + cos(2*n_vals3*tau)*(obj.amnoa0(3) + obj.anoa0(3)) + cos(2*n_vals4*tau)*(obj.amnoa0(4) + obj.anoa0(4)) + 1))^2 + (y - mu*(sin(2*n_vals1*tau)*(obj.amnoa0(1) - obj.anoa0(1)) + sin(2*n_vals2*tau)*(obj.amnoa0(2) - obj.anoa0(2)) + sin(2*n_vals3*tau)*(obj.amnoa0(3) - obj.anoa0(3)) + sin(2*n_vals4*tau)*(obj.amnoa0(4) - obj.anoa0(4))))^2 + z^2)^(3/2))))/a0^3 - (3*m^2*(2*sin(2*tau)*(x^2 - y^2) + 4*x*y*cos(2*tau)))/4;
        end

        function [q_n1, p_n1] = SI_EOM(obj, dt, scheme, X)
            q = X(1:4);
            p = X(5:end);

            dt2 = dt / 2;

            den = (dt^2*obj.m^2 + 2*dt^2*obj.m + dt^2 + 4);

            T = (4/den) * [1, dt2*(1+obj.m), 0;
                          -dt2*(1+obj.m), 1, 0;
                           0,             0, den/4];

            D = [1,  dt2*(1+obj.m), 0;
                -dt2*(1+obj.m), 1,  0;
                 0,             0,  1];

            switch scheme
                case 1 % Stormer-Verlet A
                    error('Scheme 1 not implemented yet')

                case 2 % Stormer-Verlet B
                    q_n2      = vpa(zeros(4,1));
                    q_n2(1:3) = T * (q(1:3) + dt2 * (p(1:3)));
                    q_n2(4)   = q(4) + dt/2;

                    p_n1    = vpa(zeros(4,1));
                    %dU_n  = obj.partialU(p(4),q_n2(1:3));
                    %dU_n1 = obj.partialU(p_n1(4),q_n2(1:3));
                    dU = obj.partialU(q_n2(4),q_n2);
                    p_n1(1:3) = T*( D*p(1:3,1) ...
                        - dt* ( (1+obj.m)^2*[q_n2(1);q_n2(2);0] -dU) );
                    p_n1(4) = p(4) - dt* (-obj.partialU_tau(q_n2(4),q_n2));

                    q_n1 = q_n2 + dt/2*obj.Hp([q_n2;p_n1]);
            end

        end

        function [q_n1, p_n1] = SI_EOM_nomp(obj, dt, scheme, X)
            q = X(1:4);
            p = X(5:end);

            dt2 = dt / 2;

            den = (dt^2*obj.m^2 + 2*dt^2*obj.m + dt^2 + 4);

            T = (4/den) * [1, dt2*(1+obj.m), 0;
                          -dt2*(1+obj.m), 1, 0;
                           0,             0, den/4];

            D = [1,  dt2*(1+obj.m), 0;
                -dt2*(1+obj.m), 1,  0;
                 0,             0,  1];

            switch scheme
                case 1 % Stormer-Verlet A
                    error('Scheme 1 not implemented yet')

                case 2 % Stormer-Verlet B
                    q_n2      = (zeros(4,1));
                    q_n2(1:3) = T * (q(1:3) + dt2 * (p(1:3)));
                    q_n2(4)   = q(4) + dt/2;

                    p_n1    = (zeros(4,1));
                    %dU_n  = obj.partialU(p(4),q_n2(1:3));
                    %dU_n1 = obj.partialU(p_n1(4),q_n2(1:3));
                    dU = obj.partialU(q_n2(4),q_n2);
                    p_n1(1:3) = T*( D*p(1:3,1) ...
                        - dt* ( (1+obj.m)^2*[q_n2(1);q_n2(2);0] -dU) );
                    p_n1(4) = p(4) - dt* (-obj.partialU_tau(q_n2(4),q_n2));

                    q_n1 = q_n2 + dt/2*obj.Hp([q_n2;p_n1]);
            end

        end

        function [rM,rE] = state_Earth_Moon(obj,tau)
            xi_bar = sum((obj.anoa0 + obj.amnoa0).*cos(2*obj.n_vals*tau));
            eta_bar = sum((obj.anoa0 - obj.amnoa0).*sin(2*obj.n_vals*tau));

            rM = (1 - obj.mu)*[1 + xi_bar', eta_bar', zeros(length(eta_bar),1)]';
            rE = -obj.mu*[1 + xi_bar', eta_bar', zeros(length(eta_bar),1)]';
        end

        function H = hamiltonian(obj,sol)

            s = sol.x;
            if strcmp(sol.coord,'cartesian')
                s = obj.xi2nu(s);
                q = s(1:3,:); p = s(4:end,:);
                qtau = sol.t;
                ptau = zeros(1,length(sol.t));
            else
                q = s(1:4,:); p = s(5:end,:);
                qtau = q(4,:);
                ptau = p(4,:);
            end

            px = p(1,:); py = p(2,:); pz = p(3,:);
            qx = q(1,:); qy = q(2,:); qz = q(3,:);
            
            %

            H = 1/2 * ((px + (1+obj.m)*qy).^2 + (py - (1+obj.m)*qx).^2 + pz.^2)... 
            - obj.potential(qtau,s) ;%+ ptau; 
        end

        function H = kamiltonian(obj,sol)

            s = sol.x;
            if strcmp(sol.coord,'cartesian')
                s = obj.xi2nu(s);
                q = s(1:3,:); p = s(4:end,:);
                qtau = sol.t;
                ptau = zeros(1,length(sol.t));
            else
                q = s(1:4,:); p = s(5:end,:);
                qtau = q(4,:);
                ptau = p(4,:);
            end

            px = p(1,:); py = p(2,:); pz = p(3,:);
            qx = q(1,:); qy = q(2,:); qz = q(3,:);
            
            %

            H = 1/2 * ((px + (1+obj.m)*qy).^2 + (py - (1+obj.m)*qx).^2 + pz.^2)... 
            - obj.potential(qtau,s) + ptau; 
        end

        function xi = nu2xi(obj,nu)
            switch obj.center
                case 'bary'
                    xi = obj.P_nu_xi * nu;
                case 'p2'
                    error("'Not implemented")
                case 'p1'
                    error("'Not implemented")
            end
        end

        function nu = xi2nu(obj,xi)
            switch obj.center
                case 'bary'
                    nu = obj.P_xi_nu * xi;
                case 'p2'
                    error("'Not implemented")
                case 'p1'
                    error("'Not implemented")
            end
        end

    end
end
