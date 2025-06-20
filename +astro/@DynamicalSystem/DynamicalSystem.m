%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DynamicalSystem is a class that defines the interface for dynamical systems.
% It is an abstract class that must be subclassed to be used.
% Child classes must implement the EOM, dfdx, dfdu, and add_impulse methods.
% The class provides methods for propagating the system, computing the
% state transition matrix, and discretizing the system for optimization.
% The current implementation assumes that the system is control affine.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef DynamicalSystem < dynamicprops 
    
    properties
        odeopts struct = odeset('RelTol', 1e-13, 'AbsTol', 1e-13);
        integrator function_handle = @ode78;

        mu % Mass parameter 
    end

    methods 

        function varargout = propagate(obj, x0, tspan, odeopts)
            arguments
                obj
                x0 
                tspan
                odeopts = obj.odeopts
            end

            [varargout{1:nargout}] = obj.integrator(@obj.EOM, tspan, x0, odeopts);
        end

        function varargout = propagate_with_STM(obj, x0, tspan, odeopts)
            arguments
                obj
                x0 
                tspan
                odeopts = obj.odeopts
            end
            y0 = obj.get_initial_state(x0);
            [varargout{1:nargout}] = obj.integrator(@obj.EOM_with_STM, tspan, y0, odeopts);
        end

        function [dydt, A] = EOM_with_STM(obj, t, y)
            x = y(1:obj.nx);
            STM = reshape(y(obj.nx+1:end), obj.nx, obj.nx);
            A = obj.dfdx(t, x);
            dydt = [obj.EOM(t, x); reshape(A*STM, [], 1)];
        end
        
        function [xi0, nu0] = shiftOrigin(obj, xi0, nu0, center)
            switch lower(center)
                case 'bary'
                    obj.r1 = -obj.mu;
                    obj.r2 = 1 - obj.mu;
                case 'p2'
                    xi0(1) = xi0(1) - (1 - obj.mu);
                    nu0(1) = nu0(1) - (1 - obj.mu);

                    obj.r1 = -1;
                    obj.r2 = 0;

                case 'p1'
                    xi0(1) = xi0(1)    - (- obj.mu);
                    nu0(1) = nu0(1) - (- obj.mu);

                    obj.r1 = 0;
                    obj.r2 = 1;
                otherwise
                    error('Unknown center option: %s', center);
            end
        end

    end

    methods (Abstract)
        dxdt = EOM(obj, t, x)
    end

end