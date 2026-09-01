function [q_n1,p_n1,e_q2,e_q,e_dt] = SI2_ER3BP_onetimestep_CompSum(dt,q,p,...
    U,partial_U,phi_l,vargin,e_q2,e_q,e_dt)

mu = vargin(2);
e  = vargin(end);

dt = phi_l*dt;

[den, e_dt] = comp_sum(1,e_dt,(dt/2)^2);

T = 1/den * [1         dt/2   0;...
             -dt/2       1      0;...
             0           0      den];

D = [1         dt/2    0;...
    -dt/2      1       0;...
    0            0       1];


% q_n2    = T * (q(1:3,1) + dt/2 * p(1:3,1));
delta_q2 = dt/2 * p(1:3,1);
[q_n2, e_q2] = comp_sum(T*q(1:3,1),e_q2,T*delta_q2);
q_n2(4) = q(4) - dt/2 * e/2*sin(p(4))/(1+e*cos(p(4)))^2 * ( norm(q_n2(1:3))^2 - 2*U(q_n2,mu));

p_n1 = zeros(4,1);
p_n1(4) = p(4) + dt;
p_n1(1:3) = T*( D*p(1:3,1) ...
    - dt/2 * 1/(1+e*cos(p(4)))   * (e*cos(p(4)) *q_n2(1:3,1) + partial_U(q_n2,vargin)) +...
    - dt/2 * 1/(1+e*cos(p_n1(4)))* (e*cos(p_n1(4))*q_n2(1:3,1) + partial_U(q_n2,vargin)));

q_n1 = zeros(4,1);
% q_n1(1:3) = D*q_n2(1:3,1) + dt/2 * p_n1(1:3,1);
delta_q = dt/2 * p_n1(1:3,1);
[q_n1(1:3), e_q] = comp_sum(D*q_n2(1:3,1),e_q,delta_q);

q_n1(4)   = q_n2(4) - dt/2 * e/2*sin(p_n1(4))/(1+e*cos(p_n1(4)))^2 * (norm(q_n2(1:3))^2 - 2*U(q_n2,mu));



end