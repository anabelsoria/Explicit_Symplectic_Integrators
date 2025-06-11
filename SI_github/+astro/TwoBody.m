%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Author : Anabel Soria-Carro
% Date   : May 22, 2025
% Affiliation: The University of Texas at Austin
%              Controls Group for Distributed and Uncertain Systems (CDUS)
% Description:
%  This class defines the Two-Body Problem (TBP) dynamical system.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

classdef TwoBody < astro.DynamicalSystem
    properties
        mu % Gravitational parameter
    end

    methods
        function obj = TwoBody(mu)
            % Constructor: initializes TwoBody object
            obj.mu = mu;
        end

        function dxdt = EOM(obj, t, x)
            % Equations of motion (EOM) for the two-body problem.
            % Input:
            %   t - time (not used since system is autonomous)
            %   x - state vector [position(3); velocity(3)]
            % Output:
            %   dxdt - time derivative of state [velocity; acceleration]

            r = x(1:3);
            dxdt = [x(4:6); -obj.mu/norm(r)^3*r];
        end

        function A = dfdx(obj, t, x)
            % Computes the Jacobian matrix of the system dynamics w.r.t. state.
            % Input:
            %   t - time
            %   x - state vector [position; velocity]
            % Output:
            %   A - 6x6 Jacobian matrix of partial derivatives

            r = x(1:3);
            A = [zeros(3)                                       eye(3);
                obj.mu/norm(r)^3 * (3/norm(r)^2*(r*r')-eye(3))  zeros(3)];
        end

        function dUdx = partialU(obj,q)
            % Gradient of gravitational potential with respect to position.
            % Input:
            %   q - position vector
            % Output:
            %   dUdx - gradient vector

            x = q(1); y = q(2); z = q(3);

            Ux = (obj.mu*x)/(x^2 + y^2 + z^2)^(3/2) ;
            Uy = (obj.mu*y)/(x^2 + y^2 + z^2)^(3/2);
            Uz = (obj.mu*z)/(x^2 + y^2 + z^2)^(3/2);

            dUdx = [Ux; Uy; Uz];
        end

        function [q_n1,p_n1] = SI_EOM(obj,dt,scheme,X)
            % Symplectic integrator step for the two-body problem.
            % Inputs:
            %   dt     - time step size
            %   scheme - integrator scheme identifier (1 or 2)
            %   X      - current state vector [position; momentum]
            % Outputs:
            %   q_n1   - updated position
            %   p_n1   - updated momentum

            q = X(1:3); p = X(4:6);

            switch scheme
                case 1
                    p_n2 = p    - dt/2*obj.partialU(q);
                    q_n1 = q    + dt*p_n2;
                    p_n1 = p_n2 - dt/2*obj.partialU(q_n1);
                case 2
                    q_n2 = q    + dt/2 * p;
                    p_n1 = p    - dt * obj.partialU(q_n2);
                    q_n1 = q_n2 + dt/2 * p_n1;
            end

        end

        function E = total_energy(obj,y)
            r = y(1:3,:); v = y(4:6,:);
            E = 1/2*vecnorm(v).^2 - obj.mu./vecnorm(r);
            
        end
    end
end
