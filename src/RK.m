%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Author : Anabel Soria-Carro
% Date   : May 22, 2025
% Affiliation: The University of Texas at Austin
%              Controls Group for Distributed and Uncertain Systems (CDUS)
% Description:
%  This class implements explicit Runge-Kutta integrators for solving 
%  Hamiltonian systems.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

classdef RK

    properties
        prob      % Problem class containing parameters, initial conditions, etc.
        order     % Order of Runge-Kutta method 
        stepFun   % Function handle to specific RK method
    end

    methods
        function obj = RK(prob, order)
            % Constructor for the RK class
            obj.prob  = prob;
            obj.order = order;
            obj.stepFun = RK.selectStepFunction(order);
        end

        function [X, tspan] = propagate(obj, S0, t0, tf, dt, f)
            % Propagate the system using the selected RK method
            %
            % Inputs:
            %   S0  - Initial state vector
            %   t0  - Initial time
            %   tf  - Final time
            %   dt  - Time step
            %   f   - Function handle for the dynamics
            %
            % Outputs:
            %   X     - State history (each column is the state at a time step)
            %   tspan - Time vector

            tspan = t0:dt:tf;
            nt = length(tspan);
            ns = length(S0);

            X = zeros(ns, nt);
            X(:, 1) = S0;

            for i = 2:nt
                X(:, i) = obj.stepFun(dt, tspan(i-1), X(:, i-1), f);
            end
        end
    end

    methods (Static, Access = private)
        function stepFun = selectStepFunction(order)
            % Return appropriate RK method function handle
            switch order
                case 2
                    stepFun = @RK.RK2; % Midpoint method
                case 3
                    stepFun = @RK.RK3; 
                case 4
                    stepFun = @RK.RK4; % Classical RK4
                otherwise
                    error('Runge-Kutta method of order %d not implemented.', order);
            end
        end

        function X_k1 = RK2(dt, t, X_k, f)
            % 2nd-order Runge-Kutta (Midpoint Method)
            k1 = f(t, X_k);
            k2 = f(t + dt, X_k + dt * k1);
            X_k1 = X_k + (dt/2) * (k1 + k2);
        end

        function X_k1 = RK3(dt, t, X_k, f)
            % 3rd-order Runge-Kutta 
            k1 = f(t, X_k);
            k2 = f(t + dt/2, X_k + dt/2 * k1);
            k3 = f(t + dt, X_k -dt*k1 + 2*dt * k2);
            X_k1 = X_k + (dt/6) * (k1 + 4*k2 + k3);
        end

        function X_k1 = RK4(dt, t, X_k, f)
            % 4th-order Runge-Kutta
            k1 = f(t,         X_k);
            k2 = f(t + dt/2,  X_k + (dt/2) * k1);
            k3 = f(t + dt/2,  X_k + (dt/2) * k2);
            k4 = f(t + dt,    X_k + dt * k3);
            X_k1 = X_k + (dt/6) * (k1 + 2*k2 + 2*k3 + k4);
        end
    end
end
