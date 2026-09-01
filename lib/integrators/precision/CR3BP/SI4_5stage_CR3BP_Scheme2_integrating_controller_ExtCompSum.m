function [qout,pout,tout] = SI4_5stage_CR3BP_Scheme2_integrating_controller_ExtCompSum(tf,q,p,epsilon,...
    partial_U,w,alpha,vargin,size_vecs)

t = 0;
mu = vargin(2);

z_0 = 1/(q'*q)^(alpha/2);

qout = nan*ones(size_vecs,3); pout = nan*ones(size_vecs,3);
tout = nan*ones(size_vecs,1);

qout(1,:) = q;
pout(1,:) = p;
tout(1) = t;

s3 = 1/(4-4^(1/3));
phi_l = [s3, s3, 1-4*s3, s3, s3];

G = -alpha * Hp_CR3BP_centerAtP2(mu,p,q) * q / (q'*q);
z_n2 = z_0 + epsilon/2 * G;

e_q2 = zeros(9,1);
e_q = zeros(6,1);
e_p = zeros(9,1);
e_dt2_4 = 0;
e_t = 0;

for jj = 1:5
    [q,p,t,e_q2,e_q,e_p,e_dt2_4,e_t] = SI2_CR3BP_onetimestep_Scheme2_integrating_controller_ExtCompSum(epsilon,...
        q,p,z_n2,t,w,partial_U,vargin,phi_l(jj),e_q2,e_q,e_p,e_dt2_4,e_t);
end

ii = 2;
qout(ii,:) = q; pout(ii,:) = p; tout(ii) = t;

while t < tf

    G = -alpha * Hp_CR3BP_centerAtP2(mu,p,q) * q / (q'*q);
    z_n2 = z_n2 + epsilon * G;

    for jj = 1:5
        [q,p,t,e_q2,e_q,e_p,e_dt2_4,e_t] = SI2_CR3BP_onetimestep_Scheme2_integrating_controller_ExtCompSum(epsilon,...
            q,p,z_n2,t,w,partial_U,vargin,phi_l(jj),e_q2,e_q,e_p,e_dt2_4,e_t);
    end

    ii = ii + 1;
    qout(ii,:) = q; pout(ii,:) = p; tout(ii) = t;
end

NaNIndices = find(isnan(qout));
if ~isempty(NaNIndices)
    qout = qout(1:NaNIndices(1)-1,:);
    pout = pout(1:NaNIndices(1)-1,:);
    tout = tout(1:NaNIndices(1)-1);
end
