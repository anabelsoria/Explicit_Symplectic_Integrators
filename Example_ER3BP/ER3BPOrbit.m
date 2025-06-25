%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Author : Anabel Soria-Carro
% Date   : May 22, 2025
% Affiliation: The University of Texas at Austin
%              Controls Group for Distributed and Uncertain Systems (CDUS)
% Description:
%  This class defines predefined periodic orbits in the Elliptic Restricted
%  Three-Body Problem (ER3BP). Provides initial state vectors, orbital 
%  periods, and dynamical system setup for given orbit types.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

classdef ER3BPOrbit
    properties
        type    % Orbit type string
        center  % Origin coordinates
        xi0     % Cartesian initial conditions
        nu0     % Hamiltonian initial conditions  
        Tp      % Orbit period
        DS      % Dynamical system object (i.e. ER3BP)
    end
    
    methods
        function obj = ER3BPOrbit(type,center)
            % Constructor: initialize orbit parameters based on orbit type
            
            import astro.Constants;
            import astro.ER3BP;
            
            obj.type = type;
            obj.center = center;
            
            % CR3BP characteristic properties (Earth-Moon system)
            LU = Constants.LU_EM;
            TU = Constants.TU_EM;
            mu = Constants.MU_EM;
            e  = Constants.e_EM;
            
            % Initialize ER3BP dynamical system object
            obj.DS = ER3BP(mu,obj.center,e);
            
            % Get initial conditions and period
            [obj.xi0, obj.Tp] = obj.IC_po(type);
            
            % Transform initial conditions to canonical coordinates
            obj.nu0 = obj.DS.P_xi_nu * obj.xi0;

            % Switch center
            [obj.xi0, obj.nu0] = obj.DS.shiftOrigin(obj.xi0, obj.nu0, center);

            % Add pt and qt as variables
            obj.nu0 = [obj.nu0(1:3);0;obj.nu0(4:end);0];
        end
        
        function [xi0, Tp] = IC_po(~, type)
            % Returns initial state vector and orbital period for known orbits
            %
            % Inputs:
            %   type - Orbit type (string)
            %
            % Outputs:
            %   xi0 - Initial state vector [x; y; z; vx; vy; vz]
            %   Tp  - Orbital period

            switch type
                case 'DRO_1_1'
                    x0  = 0.138804934901161;
                    y0  = 0;
                    z0  = 0;
                    vx0 = 0;
                    vy0 = 3.24515236755642;
                    vz0 = 0;
                    % Orbital resonance (M:N)
                    M   = 1; % Orbit revs
                    N   = 1; % Primaries revs
                    Tp  = 2*pi*N;
                case 'DRO_3_2'
                    x0  = 0.744738588535288;
                    y0  = 0;
                    z0  = 0;
                    vx0 = 0;
                    vy0 = 0.600143162130891;
                    vz0 = 0;
                    % Orbital resonance (M:N)
                    M   = 3; % Orbit revs
                    N   = 2; % Primaries revs
                    Tp  = 2*pi*N;
                case 'DRO_5_2'
                    x0  = 0.830003941607134;
                    y0  = 0;
                    z0  = 0;
                    vx0 = 0;
                    vy0 = 0.483496169760793;
                    vz0 = 0;
                    % Orbital resonance (M:N)
                    M   = 5; % Orbit revs
                    N   = 2; % Primaries revs
                    Tp  = 2*pi*N;
                case 'Halo_L2_2_1'
                    x0  = 1.14521893547366;
                    y0  = 0;
                    z0  = -0.160895069179558;
                    vx0 = 0;
                    vy0 = -0.220898747857751;
                    vz0 = 0;
                    % Orbital resonance (M:N)
                    M   = 2; % Orbit revs
                    N   = 1; % Primaries revs
                    Tp  = 2*pi*N;
                case 'Halo_L2_3_1'
                    x0  = 1.06370588609306;
                    y0  = 0;
                    z0  = -0.212556037429961;
                    vx0 = 0;
                    vy0 = -0.163003551068639;
                    vz0 = 0;
                    % Orbital resonance (M:N)
                    M   = 3; % Orbit revs
                    N   = 1; % Primaries revs
                    Tp  = 2*pi*N;

                otherwise
                    error('Unknown orbit type: %s', type);
            end
            xi0 = [x0, y0, z0, vx0, vy0, vz0]';
        end
    end
end
