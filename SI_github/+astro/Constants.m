
classdef Constants
    properties (Constant)

        % ------------------- Time Constants
        HOUR2SEC = 3600;
        DAY2SEC = 86400;
        YEAR2DAY = 365.25;
        YEAR2SEC = 86400 * 365.25;
        SEC2HOUR = 1/3600;
        SEC2DAY = 1/86400;
        SEC2YEAR = 1/(86400 * 365.25);

        % ------------------- Astronomical Unit
        AU = 1.49597870700E+8; % [km]

        % ------------------- CR3BP Characteristic Properties: 
        % Earth-Moon System
        LU_EM = 389703; % [km]
        TU_EM = 382981; % [s]
        MU_EM = 1.215058560962404E-2;
    end

    methods (Static)
        function body = getBodyConstants(name)
            % Returns a struct with gravitational parameter and radius for 
            % a given celestial body
            % Radius: Data from JPL/NAIF's PCK00010.TPC Planetary Constants Kernel.
            % GM: Data from JPL/NAIF's DE440.BSP Kernel, summarized here:
            % https://ssd.jpl.nasa.gov/astro_par.html
            name = lower(name); % Case-insensitive

            switch name
                case 'sun'
                    body.name = 'Sun';
                    body.radius_km = 696000.0;
                    body.mu_km3ps2 = 132712440041.9394;

                case 'mercury'
                    body.name = 'Mercury';
                    body.radius_km = 2439.7;
                    body.mu_km3ps2 = 22031.868551;

                case 'venus'
                    body.name = 'Venus';
                    body.radius_km = 6051.8;
                    body.mu_km3ps2 = 324858.592;

                case 'earth'
                    body.name = 'Earth';
                    body.radius_km = 6378.1366;
                    body.mu_km3ps2 = 398600.435507;

                case 'moon'
                    body.name = 'Moon';
                    body.radius_km = 1737.4;
                    body.mu_km3ps2 = 4902.800118;

                case 'mars'
                    body.name = 'Mars';
                    body.mu_km3ps2 = 42828.375816;

                case 'jupiter'
                    body.name = 'Jupiter';
                    body.mu_km3ps2 = 126712764.1;

                case 'saturn'
                    body.name = 'Saturn';
                    body.mu_km3ps2 = 37940584.8418;

                case 'uranus'
                    body.name = 'Uranus';
                    body.mu_km3ps2 = 5794556.4;

                case 'neptune'
                    body.name = 'Neptune';
                    body.mu_km3ps2 = 6836527.10058;

                case 'pluto'
                    body.name = 'Pluto';
                    body.mu_km3ps2 = 975.5;

                otherwise
                    error('Unknown celestial body name: %s', name);
            end
        end
    end
end
