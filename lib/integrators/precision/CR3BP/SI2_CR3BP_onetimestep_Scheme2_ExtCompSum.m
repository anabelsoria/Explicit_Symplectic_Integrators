function [x_n1,p_n1,e_x2,e_x,e_p,e_dt2_4] = SI2_CR3BP_onetimestep_Scheme2_ExtCompSum(dt,x,p,...
    partial_U,phi_l,vargin,e_x2,e_x,e_p,e_dt2_4)

dt = phi_l*dt;
hdt = dt/2;

%% dt^2 + 4
[dt2_4, e_dt2_4] = comp_sum(4,e_dt2_4,dt^2);

%% x_n2(1) = (4*x(1) + 2*dt*p(1) + 2*dt*x(2) + dt^2*p(2)) / (dt^2+4)
[s1, e_x2(1)] = comp_sum(4*x(1),e_x2(1),2*dt*p(1));
[s2, e_x2(2)] = comp_sum(s1,e_x2(2),2*dt*x(2));
[s3, e_x2(3)] = comp_sum(s2,e_x2(3),dt^2*p(2));
x_n2_1 = s3/dt2_4;

%% x_n2(2) = (-2*dt*x(1) - dt^2*p(1) + 4*x(2) + 2*dt*p(2)) / (dt^2+4)
[s1, e_x2(4)] = comp_sum(4*x(2),e_x2(4),2*dt*p(2));
[s2, e_x2(5)] = comp_sum(s1,e_x2(5),-2*dt*x(1));
[s3, e_x2(6)] = comp_sum(s2,e_x2(6),-dt^2*p(1));
x_n2_2 = s3/dt2_4;

%% x_n2(3) = x(3) + hdt*p(3)
[x_n2_3, e_x2(7)] = comp_sum(x(3),e_x2(7),hdt*p(3));

x_n2 = [x_n2_1; x_n2_2; x_n2_3];

%% Force at x_n2
dU = partial_U(x_n2,vargin);

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

%% x_n1(1) = x_n2(1) + hdt*x_n2(2) + hdt*p_n1(1)
[s1, e_x(1)] = comp_sum(x_n2_1,e_x(1),hdt*x_n2_2);
[x_n1_1, e_x(2)] = comp_sum(s1,e_x(2),hdt*p_n1_1);

%% x_n1(2) = x_n2(2) - hdt*x_n2(1) + hdt*p_n1(2)
[s1, e_x(3)] = comp_sum(x_n2_2,e_x(3),-hdt*x_n2_1);
[x_n1_2, e_x(4)] = comp_sum(s1,e_x(4),hdt*p_n1_2);

%% x_n1(3) = x_n2(3) + hdt*p_n1(3)
[x_n1_3, e_x(5)] = comp_sum(x_n2_3,e_x(5),hdt*p_n1_3);

x_n1 = [x_n1_1; x_n1_2; x_n1_3];
p_n1 = [p_n1_1; p_n1_2; p_n1_3];

end
