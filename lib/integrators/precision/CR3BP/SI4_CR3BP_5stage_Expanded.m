function [xout,pout] = SI4_CR3BP_5stage_Expanded(t0,tf,dt,S0,...
    w,partial_U,vargin)
% Five-stage SI4 using the hand-expanded scalar update WITHOUT compensation.
% Control case for isolating the cost of compensation against
% SI4_CR3BP_5stage_ExtCompSum, which evaluates the identical expressions.

if ~exist('vargin','var')
    vargin = [];
end

tspan  = t0:dt:tf;
ntspan = length(tspan);

xout = nan(3,ntspan);
pout = nan(3,ntspan);

xout(:,1) = S0(1:3)';
pout(:,1) = S0(4:6)';

x = xout(:,1);
p = pout(:,1);

s3    = 1/(4-4^(1/3));
phi_l = [s3, s3, 1-4*s3, s3, s3];

for ii = 2:ntspan
    for jj = 1:length(phi_l)
        [x,p] = SI2_CR3BP_onetimestep_Scheme2_Expanded(dt,x,p, ...
            partial_U,phi_l(jj),vargin);
    end
    xout(:,ii) = x;
    pout(:,ii) = p;
end
end
