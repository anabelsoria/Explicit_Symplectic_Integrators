function [q_n1,p_n1,t_n1,e_q2,e_q,e_p,e_dt,e_t] = SI2_CR3BP_onetimestep_Scheme2_integrating_controller_CompSum(epsilon,...
    q,p,z_n2,tn,w,partial_U,vargin,phi_l,e_q2,e_q,e_p,e_dt,e_t)

dt = epsilon/z_n2 *phi_l;

mu = vargin(2);

% den = 1 + (dt/2*w)^2;
[den, e_dt] = comp_sum(1,e_dt,(dt/2*w)^2);

T = 1/den * [1         dt/2*w   0;...
    -dt/2*w      1      0;...
    0            0      den];

D = [1         dt/2*w    0;...
    -dt/2*w      1       0;...
    0            0       1];

% q_n2 = T * (q + dt/2 * p - dt/2 * [0;1-mu;0]);
delta_q2 = dt/2 * (p - [0;1-mu;0]);
[q_n2, e_q2] = comp_sum(q,e_q2,delta_q2);
q_n2 = T*q_n2;


% p_n1 = T*( D*p - dt * partial_U(q_n2,vargin));
delta_pn1 = - dt * partial_U(q_n2,vargin);
[p_n1, e_p] = comp_sum(D*p,e_p,delta_pn1);
p_n1 = T * p_n1;

% q_n1 = D*q_n2 + dt/2 * p_n1 - dt/2 * [0;1-mu;0];
delta_qn1 = dt/2 * (p_n1 - [0;1-mu;0]);
[q_n1, e_q] = comp_sum(D*q_n2,e_q,delta_qn1);

% t_n1 = tn + dt;
[t_n1, e_t] = comp_sum(tn,e_t,dt);

end