function [xout,pout,tout] = SI4_ER3BP_5stage_CompSum(t0,tf,dt,S0,U,partial_U,vargin)

if ~exist('vargin','var')
    vargin = [];
end

tspan  = t0:dt:tf;
ntspan = length(tspan);

xout(:,1) = S0(1:4)';
pout(:,1) = S0(5:end)';

x = xout;
p = pout;

tout(1) = 0;

s3 = 1/(4-4^(1/3));

phi_l(1) = s3;
phi_l(2) = s3; 
phi_l(3) = 1- 4*s3; 
phi_l(4) = s3;
phi_l(5) = s3;

e_dt = 0;
e_q2 = zeros(3,1); e_q = zeros(3,1);
for ii = 2:ntspan

     for jj = 1:length(phi_l)
        [x,p,e_q2,e_q,e_dt] =   SI2_ER3BP_onetimestep_CompSum(dt,x,p,...
                        U,partial_U,phi_l(jj),vargin,e_q2,e_q,e_dt);
     end

% Save variables
    pout(:,ii) = p;
    xout(:,ii) = x;
    tout(ii) = tspan(ii);
end

end