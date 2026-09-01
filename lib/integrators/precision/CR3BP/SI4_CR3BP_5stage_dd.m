function [xout,pout] = SI4_CR3BP_5stage_dd(t0,tf,dt,S0,vargin)

tspan  = t0:dt:tf;
ntspan = length(tspan);

xout = nan(3,ntspan);
pout = nan(3,ntspan);

xout(:,1) = S0(1:3)';
pout(:,1) = S0(4:6)';

xh = S0(1:3);
xl = zeros(3,1);
ph = S0(4:6);
pl = zeros(3,1);

s3 = 1/(4-4^(1/3));

phi_l(1) = s3;
phi_l(2) = s3;
phi_l(3) = 1- 4*s3;
phi_l(4) = s3;
phi_l(5) = s3;

for ii = 2:ntspan

    for jj = 1:length(phi_l)
        [xh,xl,ph,pl] = SI2_CR3BP_Scheme2_dd_onetimestep(dt,xh,xl,ph,pl,...
            phi_l(jj),vargin);
    end

    pout(:,ii) = ph;
    xout(:,ii) = xh;

end

end
