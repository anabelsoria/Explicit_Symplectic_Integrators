function [x,p,e_x,e_p] = SI2_CR3BP_onetimestep_Scheme2_CompSumScalar(dt,x,p,...
    partial_U,phi_l,vargin,e_x,e_p)
% Compensated summation on the state accumulations only, in scalar form.
%
% Shares the formulation of SI2_CR3BP_onetimestep_Scheme2_Expanded and
% ..._ExtCompSum, so the three differ only in how much is compensated:
%
%   Expanded        0 compensated additions per substep
%   CompSumScalar   9, one per state component per stage        <-- this
%   ExtCompSum     22, one per partial sum
%
% Each increment is derived ALGEBRAICALLY rather than by forming the full
% updated value and subtracting the current state. Subtracting would be an
% exact operation by Sterbenz, so adding the result back would recover the
% same value and leave the compensation register empty, making the scheme a
% no-op. The precision is lost while forming the full value, not in the final
% addition, so the increment has to avoid forming it. Writing
%
%   x_n2(1) - x(1) = dt*(2*p(1) + 2*x(2) + dt*p(2) - dt*x(1)) / (dt^2+4)
%
% keeps every term of order dt and no large quantities cancel.

dt  = phi_l*dt;
hdt = dt/2;

dt2_4 = 4 + dt^2;

%% Stage 1, x <- x_n2. Increments formed before any component is updated.
dx1 = dt*(2*p(1) + 2*x(2) + dt*p(2) - dt*x(1))/dt2_4;
dx2 = dt*(2*p(2) - 2*x(1) - dt*p(1) - dt*x(2))/dt2_4;
dx3 = hdt*p(3);

[x(1), e_x(1)] = comp_sum(x(1),e_x(1),dx1);
[x(2), e_x(2)] = comp_sum(x(2),e_x(2),dx2);
[x(3), e_x(3)] = comp_sum(x(3),e_x(3),dx3);

%% Force at x_n2
dU = partial_U(x,vargin);

%% Stage 2, p <- p_n1
dp1 = dt*(-2*dt*p(1) + 4*p(2) - 4*dU(1) - 2*dt*dU(2))/dt2_4;
dp2 = dt*(-4*p(1) - 2*dt*p(2) + 2*dt*dU(1) - 4*dU(2))/dt2_4;
dp3 = -dt*dU(3);

[p(1), e_p(1)] = comp_sum(p(1),e_p(1),dp1);
[p(2), e_p(2)] = comp_sum(p(2),e_p(2),dp2);
[p(3), e_p(3)] = comp_sum(p(3),e_p(3),dp3);

%% Stage 3, x <- x_n1
dxx1 = hdt*(x(2) + p(1));
dxx2 = hdt*(p(2) - x(1));
dxx3 = hdt*p(3);

[x(1), e_x(1)] = comp_sum(x(1),e_x(1),dxx1);
[x(2), e_x(2)] = comp_sum(x(2),e_x(2),dxx2);
[x(3), e_x(3)] = comp_sum(x(3),e_x(3),dxx3);

end
