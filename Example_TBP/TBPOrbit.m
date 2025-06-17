%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Author : Anabel Soria-Carro
% Date   : May 22, 2025
% Affiliation: The University of Texas at Austin
%              Controls Group for Distributed and Uncertain Systems (CDUS)
% Description:
%  This class defines predefined periodic orbits in the Two-Body Problem
%  (TBP). Provides initial state vectors, orbital
%  periods, and dynamical system setup for given orbit types.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

classdef TBPOrbit
    properties
        S0      % Initial state vector
        S0_qp
        Tp      % Orbit period
        DS      % Dynamical system object (e.g. CR3BP)
        oe
        body      % Body properties
        mu
    end

    methods
        function obj = TBPOrbit(body,type)
            % Constructor: initialize orbit parameters

            import astro.Constants;
            import astro.TwoBody;

            obj.body = Constants.getBodyConstants(body);

            % TBP characteristic properties
            obj.mu = 1;

            % Initialize TBP dynamical system object
            obj.DS = TwoBody(obj.mu);

            % Get initial conditions and period
            [obj.S0, obj.oe, obj.Tp] = obj.IC_po(type);

            obj.S0_qp = obj.S0;

        end

        function [S0, oe, Tp] = IC_po(obj, type)
            % Returns initial state vector and orbital period for known orbits
            %
            % Inputs:
            %   type - Orbit type (string)
            %
            % Outputs:
            %   S0 - Initial state vector [x; y; z; vx; vy; vz]
            %   Tp - Orbital period

            switch type
                case 1 % Planar Orbit, ecc = 0.1
                    e = 0.1;
                    x0 = 1-e; y0 = 0; z0 = 0;
                    vx0 = 0; vy0 = sqrt((1+e)/(1-e)); vz0 = 0;
                    rvec = [x0;y0;z0]; vvec = [vx0;vy0;vz0];
                    S0 = [x0, y0, z0, vx0, vy0, vz0]';
                    oe = astro.conics.cart2coe(S0,obj.mu,'MA');
                    Tp = 2*pi*sqrt(oe.sma^3/obj.mu); % [TU]

                case 2 % Planar Orbit, ecc = 0.7
                    e = 0.7;
                    x0 = 1-e; y0 = 0; z0 = 0;
                    vx0 = 0; vy0 = sqrt((1+e)/(1-e)); vz0 = 0;
                    rvec = [x0;y0;z0]; vvec = [vx0;vy0;vz0];
                    S0 = [x0, y0, z0, vx0, vy0, vz0]';
                    oe = astro.conics.cart2coe(S0,obj.mu,'MA');
                    Tp = 2*pi*sqrt(oe.sma^3/obj.mu); % [TU]

                otherwise
                    error('Unknown orbit type: %s', type);
            end

        end
    end
end
