%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Author : Anabel Soria-Carro
% Date   : May 22, 2025
% Affiliation: The University of Texas at Austin
%              Controls Group for Distributed and Uncertain Systems (CDUS)
% Description:
%  This class implements explicit Runge-Kutta integrators for solving 
%  Hamiltonian systems.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

classdef RK < Integrator

    properties
        order     % Order of Runge-Kutta method 
        stepFun   % Function handle to specific RK method
    end

    methods
        function obj = RK(prob, order)
            % Constructor for the RK class
            obj.name = 'RK';
            obj.prob  = prob;
            obj.order = order;
            obj.stepFun = RK.selectStepFunction(order);
        end

        function varargout = propagate(obj, S0, t0, tf, dt, f)
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

            obj.sol.x = X;
            obj.sol.t = tspan;
            obj.sol.coord = 'hamiltonian'; % Modify later to option cart

            if nargout > 0
                varargout{1} = X;
                varargout{2} = tspan;
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
                case 6
                    stepFun = @RK.RK6;
                case 8
                    stepFun = @RK.RK8;
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

        function X_k1 = RK6(dt, t, X_k, f)
            % 6th-order Runge-Kutta
            % A-57 from https://ntrs.nasa.gov/api/citations/19730024029/downloads/19730024029.pdf
            k_1 = f(t        ,X_k) ;
            k_2 = f(t+(1/3)*dt,X_k+(1/3)*k_1*dt) ;
            k_3 = f(t+(2/3)*dt,X_k+ dt*(2/3*k_2)) ;
            k_4 = f(t+(1/3)*dt,X_k+ dt*(1/12*k_1 + 1/3*k_2  - 1/12*k_3)) ;
            k_5 = f(t+1/2*dt,  X_k+ dt*(-1/16*k_1 + 9/8*k_2 -3/16*k_3 -3/8*k_4)) ;
            k_6 = f(t+1/2*dt,  X_k+dt*(9/8*k_2 -3/8*k_3 -3/4*k_4 + 1/2*k_5)) ;
            k_7 = f(t+dt,  X_k+dt*(9/44*k_1 - 9/11*k_2 + 63/44*k_3 + 18/11*k_4 -16/11*k_6)) ;
            X_k1 = X_k +  dt*(11/120*k_1 + 27/40*k_3 + 27/40*k_4 -4/15*k_5 -4/15*k_6 + 11/120*k_7);
        end

        function X_k1 = RK8(dt, t, X_k, f)
            % 8th-order Runge-Kutta
            % Table 6-1 from Goddard Trajectory Determination System (GTDS): Mathematical Theory,
            % Goddard Space Flight Center, 1989.
            k_1 = f(t         ,X_k                                                                           );
            k_2 = f(t+dt*(4/27),X_k+(dt*4/27)*k_1                                                              );
            k_3 = f(t+dt*(2/9) ,X_k+  (dt/18)*(k_1+3*k_2)                                                      );
            k_4 = f(t+dt*(1/3) ,X_k+  (dt/12)*(k_1+3*k_3)                                                      );
            k_5 = f(t+dt*(1/2) ,X_k+   (dt/8)*(k_1+3*k_4)                                                      );
            k_6 = f(t+dt*(2/3) ,X_k+  (dt/54)*(13*k_1-27*k_3+42*k_4+8*k_5)                                     );
            k_7 = f(t+dt*(1/6) ,X_k+(dt/4320)*(389*k_1-54*k_3+966*k_4-824*k_5+243*k_6)                         );
            k_8 = f(t+dt       ,X_k+  (dt/20)*(-231*k_1+81*k_3-1164*k_4+656*k_5-122*k_6+800*k_7)               );
            k_9 = f(t+dt*(5/6) ,X_k+ (dt/288)*(-127*k_1+18*k_3-678*k_4+456*k_5-9*k_6+576*k_7+4*k_8)            );
            k_10= f(t+dt       ,X_k+(dt/820)*(1481*k_1-81*k_3+7104*k_4-3376*k_5+72*k_6-5040*k_7-60*k_8+720*k_9));
            X_k1 = X_k + dt/840*(41*k_1+27*k_4+272*k_5+27*k_6+216*k_7+216*k_9+41*k_10);
        end

    end
end
