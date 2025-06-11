%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Author : Anabel Soria-Carro
% Date   : May 22, 2025
% Affiliation: The University of Texas at Austin
%              Controls Group for Distributed and Uncertain Systems (CDUS)
% Description:
%  This class defines a symplectic integrator framework with configurable 
%  order and scheme
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

classdef SI

    properties
        prob    % Problem definition object (contains initial conditions, dynamics, etc.)
        order   % Order of the symplectic integrator (e.g., 2, 3, or 4)
        scheme  % Scheme selector for equations of motion
        gamma   % Gamma coefficients specific to the chosen integrator order
    end

    methods

        function obj = SI(prob,order,scheme)
            % Constructor for the SI class
            obj.prob     = prob;
            obj.order = order;
            obj.scheme = scheme;
            obj.gamma = obj.generateGamma(order);
        end

        function gamma = generateGamma(~, order)
            % Generates the gamma coefficients based on integrator order
            % These coefficients define the splitting method for time integration
            switch order
                case 2
                    gamma = 1;
                case 3
                    s = 1 / (2 - 2^(1/3));
                    gamma = [s; 1 - 2*s; s];
                case 4
                    s = 1 / (4 - 4^(1/3));
                    gamma = [s; s; 1 - 4*s; s; s];
                otherwise
                    error('Gamma coefficients for order %d not implemented.', order);
            end
        end

        function [X,tspan] = propagate(obj,t0,tf,dt)
            % Propagate the system dynamics from t0 to tf using symplectic
            % integration
            % Inputs:
            %   t0   - Initial time
            %   tf   - Final time
            %   dt   - Time step
            % Outputs:
            %   X     - State history [q; p] over time
            %   tspan - Time vector

            tspan  = t0:dt:tf;
            nt = length(tspan);
            ns = length(obj.prob.S0);

            X = zeros(ns, nt);
            X(:,1) = obj.prob.S0_qp;
            q = X(1:3,1); p = X(4:6,1);

            for ii = 2:nt
                for jj = 1:length(obj.gamma)
                    [q,p] = obj.prob.DS.SI_EOM(obj.gamma(jj)*dt,obj.scheme,[q;p]);
                end
                X(1:3,ii) = q;
                X(4:6,ii) = p;
            end
        end
    end
end
