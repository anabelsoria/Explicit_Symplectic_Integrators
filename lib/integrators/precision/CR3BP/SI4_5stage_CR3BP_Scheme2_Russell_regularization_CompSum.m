function [qout,pout,tout] = SI4_5stage_CR3BP_Scheme2_Russell_regularization_CompSum(tf,q,p,epsilon,...
    partial_U,w,vargin,vargin2,size_vecs)

t = 0;
mu = vargin(2);

A_E = vargin2(1).A; B_E = vargin2(1).B; C_E = vargin2(1).C;
rH_E = vargin2(1).rH; alpha_E = vargin2(1).alpha;
A_M = vargin2(2).A; B_M = vargin2(2).B; C_M = vargin2(2).C;
rH_M = vargin2(2).rH; alpha_M = vargin2(2).alpha;

rE = norm(q-[-1,0,0]');
rM = norm(q);
rho_1 = rho_fun(A_E,B_E,C_E,rE,rH_E);
rho_2 = rho_fun(A_M,B_M,C_M,rM,rH_M);
g = rho_1^alpha_E * rho_2^alpha_M;
z_0 = 1/g;

qout = nan*ones(size_vecs,3); pout = nan*ones(size_vecs,3);
tout = nan*ones(size_vecs,1);

qout(1,:) = q;
pout(1,:) = p;
tout(1) = t;

s3 = 1/(4-4^(1/3));
phi_l = [s3, s3, 1-4*s3, s3, s3];

G = G_Russell(q,p,vargin2);
z_n2 = z_0 + epsilon/2 * G;

e_q = zeros(3,1);
e_p = zeros(3,1);
e_q2 = zeros(3,1);
e_dt = 0; e_t = 0;

for jj = 1:5
    [q,p,t,e_q2,e_q,e_p,e_dt,e_t] = SI2_CR3BP_onetimestep_Scheme2_integrating_controller_CompSum(epsilon,...
        q,p,z_n2,t,w,partial_U,vargin,phi_l(jj),e_q2,e_q,e_p,e_dt,e_t);
end

ii = 2;
qout(ii,:) = q; pout(ii,:) = p; tout(ii) = t;

while t < tf

    G = G_Russell(q,p,vargin2);
    z_n2 = z_n2 + epsilon * G;

    for jj = 1:5
        [q,p,t,e_q2,e_q,e_p,e_dt,e_t] = SI2_CR3BP_onetimestep_Scheme2_integrating_controller_CompSum(epsilon,...
            q,p,z_n2,t,w,partial_U,vargin,phi_l(jj),e_q2,e_q,e_p,e_dt,e_t);
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
