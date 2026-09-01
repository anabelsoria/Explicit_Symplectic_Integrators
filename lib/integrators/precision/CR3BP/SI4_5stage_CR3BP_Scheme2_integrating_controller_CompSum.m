function [qout,pout,tout] = SI4_5stage_CR3BP_Scheme2_integrating_controller_CompSum(tf,q,p,epsilon,...
    partial_U,w,alpha,vargin,size_vecs)


t = 0;
%z_nm2_v2 = (q'*q)^(alpha/2);

mu = vargin(2);

z_0 = 1/(q'*q)^(alpha/2);

qout = nan*ones(size_vecs,3); pout = nan*ones(size_vecs,3); 
tout = nan*ones(size_vecs,1); zout = nan*ones(size_vecs,1);

qout(1,:) = q;
pout(1,:) = p;
tout(1) = t;
zout(1) = z_0;

s3 = 1/(4-4^(1/3));

phi_l(1) = s3;
phi_l(2) = s3; 
phi_l(3) = 1- 4*s3; 
phi_l(4) = s3;
phi_l(5) = s3;


G = -alpha * Hp_CR3BP_centerAtP2(mu,p,q) * q / (q'*q); % Verified
z_n2 = z_0 + epsilon/2 * G;

e_q = zeros(3,1);  
e_p = zeros(3,1);  
e_q2 = zeros(3,1); 
e_dt = 0; e_t = 0;

for jj = 1:length(phi_l)
        [q,p,t,e_q2,e_q,e_p,e_dt,e_t] =   SI2_CR3BP_onetimestep_Scheme2_integrating_controller_CompSum(epsilon,...
    q,p,z_n2,t,w,partial_U,vargin,phi_l(jj),e_q2,e_q,e_p,e_dt,e_t);
end

ii = 2;
qout(ii,:) = q;
pout(ii,:) = p;
tout(ii,:) = t;
zout(ii) = z_n2;

while t < tf

    G = -alpha * Hp_CR3BP_centerAtP2(mu,p,q) * q / (q'*q);
    z_n2 = z_n2 + epsilon * G;

    for jj = 1:length(phi_l)
        [q,p,t,e_q2,e_q,e_p,e_dt,e_t] =   SI2_CR3BP_onetimestep_Scheme2_integrating_controller_CompSum(epsilon,...
                q,p,z_n2,t,w,partial_U,vargin,phi_l(jj),e_q2,e_q,e_p,e_dt,e_t);
    end


    ii = ii + 1;
     % Save variables
    qout(ii,:) = q;
    pout(ii,:) = p;
    tout(ii,:) = t;
    zout(ii) = z_n2;
end

% Trim the Nan Values 
NaNIndices = find(isnan(qout));
if ~isempty(NaNIndices) 
qout = qout(1:NaNIndices(1)-1,:);
pout = pout(1:NaNIndices(1)-1,:);
tout = tout(1:NaNIndices(1)-1);
zout = zout(1:NaNIndices(1)-1);
end
