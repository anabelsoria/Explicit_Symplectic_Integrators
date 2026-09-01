function [x_n1,p_n1] = SI2_CR3BP_onetimestep_Scheme2_Expanded(dt,x,p,...
    partial_U,phi_l,vargin)
% Uncompensated twin of SI2_CR3BP_onetimestep_Scheme2_ExtCompSum.
%
% Exactly the same hand-expanded scalar expressions, evaluated in the same
% order, with every comp_sum(a,e,b) replaced by the plain addition a + b.
%
% This is the control needed to isolate the cost of the compensation. Timing
% the extended scheme against SI2_CR3BP_onetimestep_Scheme2 instead compares
% two different formulations of the update, since that one builds two 3x3
% matrices and forms four matrix-vector products at every substep.

dt  = phi_l*dt;
hdt = dt/2;

dt2_4 = 4 + dt^2;

%% x_n2(1) = (4*x(1) + 2*dt*p(1) + 2*dt*x(2) + dt^2*p(2)) / (dt^2+4)
s1 = 4*x(1) + 2*dt*p(1);
s2 = s1 + 2*dt*x(2);
s3 = s2 + dt^2*p(2);
x_n2_1 = s3/dt2_4;

%% x_n2(2) = (-2*dt*x(1) - dt^2*p(1) + 4*x(2) + 2*dt*p(2)) / (dt^2+4)
s1 = 4*x(2) + 2*dt*p(2);
s2 = s1 - 2*dt*x(1);
s3 = s2 - dt^2*p(1);
x_n2_2 = s3/dt2_4;

%% x_n2(3) = x(3) + hdt*p(3)
x_n2_3 = x(3) + hdt*p(3);

x_n2 = [x_n2_1; x_n2_2; x_n2_3];

%% Force at x_n2
dU = partial_U(x_n2,vargin);

%% p_n1(1)
s1 = dt^2*p(1) - 4*p(1);
s2 = s1 - 4*dt*p(2);
s3 = s2 + 4*dt*dU(1);
s4 = s3 + 2*dt^2*dU(2);
p_n1_1 = -s4/dt2_4;

%% p_n1(2)
s1 = 4*dt*p(1) - 4*p(2);
s2 = s1 + dt^2*p(2);
s3 = s2 - 2*dt^2*dU(1);
s4 = s3 + 4*dt*dU(2);
p_n1_2 = -s4/dt2_4;

%% p_n1(3) = p(3) - dt*dU(3)
p_n1_3 = p(3) - dt*dU(3);

%% x_n1(1) = x_n2(1) + hdt*x_n2(2) + hdt*p_n1(1)
s1 = x_n2_1 + hdt*x_n2_2;
x_n1_1 = s1 + hdt*p_n1_1;

%% x_n1(2) = x_n2(2) - hdt*x_n2(1) + hdt*p_n1(2)
s1 = x_n2_2 - hdt*x_n2_1;
x_n1_2 = s1 + hdt*p_n1_2;

%% x_n1(3) = x_n2(3) + hdt*p_n1(3)
x_n1_3 = x_n2_3 + hdt*p_n1_3;

x_n1 = [x_n1_1; x_n1_2; x_n1_3];
p_n1 = [p_n1_1; p_n1_2; p_n1_3];

end
