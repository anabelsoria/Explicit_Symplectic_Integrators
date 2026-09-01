function [q_n1,p_n1,t_n1,e_q2,e_q,e_p,e_dt2_4,e_t] = SI2_CR3BP_onetimestep_Scheme2_integrating_controller_ExtCompSum(epsilon,...
    q,p,z_n2,tn,w,partial_U,vargin,phi_l,e_q2,e_q,e_p,e_dt2_4,e_t)

dt = epsilon/z_n2 * phi_l;
hdt = dt/2;
mu = vargin(2);
off2 = hdt*(1-mu);

%% dt^2 + 4
[dt2_4, e_dt2_4] = comp_sum(4,e_dt2_4,dt^2);

%% q_n2(1) = (4*q(1) + 2*dt*p(1) + 2*dt*q(2) + dt^2*p(2) - dt^2*(1-mu)) / (dt^2+4)
[s1, e_q2(1)] = comp_sum(4*q(1),e_q2(1),2*dt*p(1));
[s2, e_q2(2)] = comp_sum(s1,e_q2(2),2*dt*q(2));
[s3, e_q2(3)] = comp_sum(s2,e_q2(3),dt^2*p(2));
[s4, e_q2(4)] = comp_sum(s3,e_q2(4),-dt^2*(1-mu));
q_n2_1 = s4/dt2_4;

%% q_n2(2) = (-2*dt*q(1) - dt^2*p(1) + 4*q(2) + 2*dt*p(2) - 2*dt*(1-mu)) / (dt^2+4)
[s1, e_q2(5)] = comp_sum(4*q(2),e_q2(5),2*dt*p(2));
[s2, e_q2(6)] = comp_sum(s1,e_q2(6),-2*dt*q(1));
[s3, e_q2(7)] = comp_sum(s2,e_q2(7),-dt^2*p(1));
[s4, e_q2(8)] = comp_sum(s3,e_q2(8),-2*dt*(1-mu));
q_n2_2 = s4/dt2_4;

%% q_n2(3) = q(3) + hdt*p(3)
[q_n2_3, e_q2(9)] = comp_sum(q(3),e_q2(9),hdt*p(3));

q_n2 = [q_n2_1; q_n2_2; q_n2_3];

%% Force at q_n2
dU = partial_U(q_n2,vargin);

%% p_n1(1) = -(dt^2*p(1) - 4*p(1) - 4*dt*p(2) + 4*dt*dU(1) + 2*dt^2*dU(2)) / (dt^2+4)
[s1, e_p(1)] = comp_sum(dt^2*p(1),e_p(1),-4*p(1));
[s2, e_p(2)] = comp_sum(s1,e_p(2),-4*dt*p(2));
[s3, e_p(3)] = comp_sum(s2,e_p(3),4*dt*dU(1));
[s4, e_p(4)] = comp_sum(s3,e_p(4),2*dt^2*dU(2));
p_n1_1 = -s4/dt2_4;

%% p_n1(2) = -(4*dt*p(1) - 4*p(2) + dt^2*p(2) - 2*dt^2*dU(1) + 4*dt*dU(2)) / (dt^2+4)
[s1, e_p(5)] = comp_sum(4*dt*p(1),e_p(5),-4*p(2));
[s2, e_p(6)] = comp_sum(s1,e_p(6),dt^2*p(2));
[s3, e_p(7)] = comp_sum(s2,e_p(7),-2*dt^2*dU(1));
[s4, e_p(8)] = comp_sum(s3,e_p(8),4*dt*dU(2));
p_n1_2 = -s4/dt2_4;

%% p_n1(3) = p(3) - dt*dU(3)
[p_n1_3, e_p(9)] = comp_sum(p(3),e_p(9),-dt*dU(3));

%% q_n1(1) = q_n2(1) + hdt*q_n2(2) + hdt*p_n1(1)
[s1, e_q(1)] = comp_sum(q_n2_1,e_q(1),hdt*q_n2_2);
[q_n1_1, e_q(2)] = comp_sum(s1,e_q(2),hdt*p_n1_1);

%% q_n1(2) = q_n2(2) - hdt*q_n2(1) + hdt*p_n1(2) - off2
[s1, e_q(3)] = comp_sum(q_n2_2,e_q(3),-hdt*q_n2_1);
[s2, e_q(4)] = comp_sum(s1,e_q(4),hdt*p_n1_2);
[q_n1_2, e_q(5)] = comp_sum(s2,e_q(5),-off2);

%% q_n1(3) = q_n2(3) + hdt*p_n1(3)
[q_n1_3, e_q(6)] = comp_sum(q_n2_3,e_q(6),hdt*p_n1_3);

%% time
[t_n1, e_t] = comp_sum(tn,e_t,dt);

q_n1 = [q_n1_1; q_n1_2; q_n1_3];
p_n1 = [p_n1_1; p_n1_2; p_n1_3];

end
