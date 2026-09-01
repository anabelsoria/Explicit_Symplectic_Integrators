function [xout,pout] = SI4_CR3BP_5stage_ExtCompSum(t0,tf,dt,S0,...
    w,partial_U,vargin)

x = S0(1:3);
p = S0(4:6);

nsteps = round((tf-t0)/dt);

xout = nan*ones(3,nsteps+1);
pout = nan*ones(3,nsteps+1);

xout(:,1) = x;
pout(:,1) = p;

s3 = 1/(4-4^(1/3));
phi_l = [s3, s3, 1-4*s3, s3, s3];

e_x2 = zeros(7,1);
e_x = zeros(5,1);
e_p = zeros(9,1);
e_dt2_4 = 0;

for ii = 1:nsteps
    for jj = 1:5
        [x,p,e_x2,e_x,e_p,e_dt2_4] = SI2_CR3BP_onetimestep_Scheme2_ExtCompSum(dt,x,p,...
            partial_U,phi_l(jj),vargin,e_x2,e_x,e_p,e_dt2_4);
    end
    xout(:,ii+1) = x;
    pout(:,ii+1) = p;
end
