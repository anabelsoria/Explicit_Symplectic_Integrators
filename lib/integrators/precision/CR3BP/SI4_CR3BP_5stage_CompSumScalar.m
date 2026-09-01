function [xout,pout] = SI4_CR3BP_5stage_CompSumScalar(t0,tf,dt,S0,...
    w,partial_U,vargin)
% Five-stage SI4 with compensated summation on the state accumulations only,
% written in the same hand-expanded scalar form as SI4_CR3BP_5stage_Expanded
% and SI4_CR3BP_5stage_ExtCompSum so that the three are timed on equal terms.

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

e_x = zeros(3,1);
e_p = zeros(3,1);

for ii = 2:ntspan
    for jj = 1:length(phi_l)
        [x,p,e_x,e_p] = SI2_CR3BP_onetimestep_Scheme2_CompSumScalar( ...
            dt,x,p,partial_U,phi_l(jj),vargin,e_x,e_p);
    end
    xout(:,ii) = x;
    pout(:,ii) = p;
end
end
