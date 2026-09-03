function sv = coe2cart(coe,MU,anomaly)
    % assumes all angles are in radians! 

    import astro.conics.*; % import astro package

    % a = coe(1);
    % e = coe(2);
    % i = coe(3);
    % W = coe(4);
    % w = coe(5);
    a = coe.sma;
    e = coe.ecc;
    i = coe.inc_r;
    W = coe.RAAN_r;
    w = coe.omega_r;

    switch nargin
        case 2
            TA = coe(6);
        case 3
            if strcmp(anomaly,'MA')
                %MA = coe(6);
                MA = coe.nu_r;
                TA = MAtoTA(MA, e);
            elseif strcmp(anomaly, 'EA')
                % EA = coe(6);
                EA = coe.nu_r;
                TA = EAtoTA(EA, e);
            elseif strcmp(anomaly, 'TA')
                %TA = coe(6);
                TA = coe.nu_r;
            else
                error('Invalid anomaly type given.')
            end
        otherwise
            error('Number of input arguments must be 2 or 3.')
    end

    h = sqrt(a*MU*(1-e^2));
    r_pf = h^2/MU / (1+e*cos(TA)) * [cos(TA); sin(TA); 0];
    v_pf = MU/h * [-sin(TA); e+cos(TA); 0];

    r = Rmat(3, -W)*Rmat(1, -i)*Rmat(3, -w) * r_pf;
    v = Rmat(3, -W)*Rmat(1, -i)*Rmat(3, -w) * v_pf;
    sv = [r;v];
end
