%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Author : Anabel Soria-Carro
% Date   : May 22, 2025
% Affiliation: The University of Texas at Austin
%              Controls Group for Distributed and Uncertain Systems (CDUS)
% Description:
%  This class defines predefined periodic orbits in the Circular Restricted
%  Three-Body Problem (CR3BP). Provides initial state vectors, orbital 
%  periods, and dynamical system setup for given orbit types.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

classdef CR3BPOrbit
    properties
        type    % Orbit type string
        center  % Origin coordinates
        xi0     % Initial state vector
        nu0     % Transformed initial conditions (optional)
        Tp      % Orbit period
        DS      % Dynamical system object (e.g. CR3BP)
    end
    
    methods
        function obj = CR3BPOrbit(type,center)
            % Constructor: initialize orbit parameters based on orbit type
            
            import astro.Constants;
            import astro.CR3BP;
            
            obj.type = type;
            obj.center = center;
            
            % CR3BP characteristic properties (Earth-Moon system)
            LU = Constants.LU_EM;
            TU = Constants.TU_EM;
            mu = Constants.MU_EM;
            
            % Initialize CR3BP dynamical system object
            obj.DS = CR3BP(mu,obj.center);
            
            % Get initial conditions and period
            [obj.xi0, obj.Tp] = obj.IC_po(type);
            
            % Transform initial conditions to canonical coordinates
            obj.nu0 = CR3BP.P_xi_nu * obj.xi0;

            % Switch center
            [obj.xi0, obj.nu0] = obj.DS.shiftOrigin(obj.xi0, obj.nu0, center);
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
                case 'DRO'
                    x0  = 3.8680865381232921E-1;
                    y0  = 6.8407409893399016E-24;
                    z0  = -5.6003346392265899E-24;
                    vx0 = 1.0915552878170719E-12;
                    vy0 = 1.6044630974809673E+0;
                    vz0 = 1.3293950271456153E-23;
                    Tp = 6.1543653112844199E+0;
                case 'NRHO_L2_S'
                    x0  = 1.0137267853351086E+0;
                    y0  = 6.6143443622012310E-28;
                    z0  = -1.7564926055839319E-1;
                    vx0 = -3.1488714679618939E-13;
                    vy0 = -8.4451795580222144E-2;
                    vz0 = -4.9619055939541527E-12;
                    Tp = 1.4004896487302059E+0;
                case 'Halo_L1_N'
                    x0  = 8.7606564511706009E-1;
                    y0  = 1.0491502986863478E-26;
                    z0  = 1.9181353568854601E-1;
                    vx0 = -3.9282267868172261E-14;
                    vy0 = 2.3055753889250671E-1	;
                    vz0 = 1.0601974505638561E-13;
                    Tp = 2.1764730139006754E+0;
                otherwise
                    error('Unknown orbit type: %s', type);
            end
            xi0 = [x0, y0, z0, vx0, vy0, vz0]';
        end
    end
end
