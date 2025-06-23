%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Author : Anabel Soria-Carro
% Date   : June 20, 2025
% Affiliation: The University of Texas at Austin
%              Controls Group for Distributed and Uncertain Systems (CDUS)
% Description:
%  This class defines predefined periodic orbits in the Bicircular Restricted
%  Four-Body Problem (BCR4BP). Provides initial state vectors, orbital 
%  periods, and dynamical system setup for given orbit types.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

classdef BCR4BPOrbit < handle
    properties
        type    % Orbit type string
        center  % Origin coordinates
        xi0     % Cartesian initial conditions
        nu0     % Hamiltonian initial conditions 
        theta_0 % Initial Sun-B1 angle
        Tp      % Orbit period
        DS      % Dynamical system object (i.e. BCR4BP)
    end
    
    methods
        function obj = BCR4BPOrbit(type,center,theta_0)
            arguments
                type
                center 
                theta_0 = 0
            end

            % Constructor: initialize orbit parameters based on orbit type
            
            import astro.Constants;
            import astro.BCR4BP;
            
            obj.type = type;
            obj.center = center;
            
            % Bodies 
            moon  = Constants.getBodyConstants('moon');
            earth = Constants.getBodyConstants('earth');
            sun   = Constants.getBodyConstants('sun');

            mu_em = moon.m/(earth.m + moon.m);
            mu_s  = sun.m/(earth.m + moon.m);
            
            % BCR4BP characteristic properties 
            LU = Constants.LU_EM;
            
            % Initialize BCR4BP dynamical system object
            obj.DS = BCR4BP(mu_em,mu_s,obj.center,LU,sun.m, theta_0);

            obj.theta_0   = theta_0;

            % Get initial conditions and period
            [obj.xi0, obj.Tp] = obj.IC_po(type,Constants);
            
            % Transform initial conditions to canonical coordinates
            obj.nu0 = obj.DS.P_xi_nu * obj.xi0;

            % Switch center
            [obj.xi0, obj.nu0] = obj.DS.shiftOrigin(obj.xi0, obj.nu0, center);

            % Add pt and qt as variables
            obj.nu0 = [obj.nu0(1:3);0;obj.nu0(4:end);0];

            % Add theta_0 as variable to propagate
            obj.xi0 = [obj.xi0; obj.theta_0];
            % obj.nu0 = [obj.nu0; obj.theta_0];
        end
        
        function [xi0, Tp] = IC_po(obj, type, Constants)
            % Returns initial state vector and orbital period for known orbits
            %
            % Inputs:
            %   type - Orbit type (string)
            %
            % Outputs:
            %   xi0 - Initial state vector [x; y; z; vx; vy; vz]
            %   Tp  - Orbital period

            switch type
                case 'NRHO_3_1'
                    x0  = 1.086763230129224;             
                    y0  = -0.009790900673699;
                    z0  = -0.174577385347935;
                    vx0 = 0.012341101802055;
                    vy0 = -0.234896886697873;
                    vz0 = -0.076421195994642;
                    obj.theta_0 = -0.7854;
                    
                    % Orbital resonance (M:N)
                    M   = 3; % Earth-Moons revs
                    N   = 1; % Sun-B1 revs
                    Tp_3bp = 2.2635;
                    Tp  = M*Tp_3bp;

                    obj.DS.LU = 384400; % [km]
                    obj.DS.a_s = (Constants.AU/obj.DS.LU) - obj.DS.mu_em;
                    obj.DS.theta_dot = sqrt((1+obj.DS.mu_s)/obj.DS.a_s^3) - 1; 
                
                case 'NRHO_9_2'
                    x0  = 1.002805174989106;               
                    y0  = 0.019209365447395;
                    z0  = -0.168464000804656;
                    vx0 = 0.023898460222743;
                    vy0 = -0.081252704013115;
                    vz0 = -0.118918741044589;
                    obj.theta_0 = 0.1904;
                    obj.DS.theta_0 = 0.1904;

                    % Orbital resonance (M:N)
                    M   = 9; % Earth-Moons revs
                    N   = 2; % Sun-B1 revs
                    Tp_3bp = 1.5090;
                    Tp  = M*Tp_3bp;

                    obj.DS.LU = 384400; % [km]
                    obj.DS.a_s = (Constants.AU/obj.DS.LU) - obj.DS.mu_em;
                    obj.DS.theta_dot = sqrt((1+obj.DS.mu_s)/obj.DS.a_s^3) - 1;

                otherwise
                    error('Unknown orbit type: %s', type);
            end
            xi0 = [x0, y0, z0, vx0, vy0, vz0]';
        end
    end
end
