function [xout,pout] = SI4_CR3BP_5stage_Increment(t0,tf,dt,S0,...
    w,partial_U,vargin)
% Five-stage SI4 with the update written as an analytic increment, and the
% increment accumulated by ordinary addition.
%
% Each stage is obtained by subtracting the current state from the closed
% form solution SYMBOLICALLY, so that the terms of order one cancel on paper
% and are never formed. For the first component,
%
%   x_n2(1)      = (4*x(1) + 2*dt*p(1) + 2*dt*x(2) + dt^2*p(2)) / (dt^2+4)
%   x_n2(1)-x(1) =  dt*(2*p(1) + 2*x(2) + dt*p(2) - dt*x(1)) / (dt^2+4)
%
% where the 4*x(1) terms cancel exactly. Every surviving term carries a
% factor dt, so the quantity computed is of order dt rather than of order
% one, and it is therefore resolved to full relative precision.
%
% This isolates the effect of the increment formulation alone. The companion
% files add compensated summation (CompSumScalar) or a double-double
% accumulation (ddInc) on top of the same increments.

if ~exist('vargin','var')
    vargin = [];
end

tspan  = t0:dt:tf;
ntspan = length(tspan);

xout = nan(3,ntspan);
pout = nan(3,ntspan);

xout(:,1) = S0(1:3)';
pout(:,1) = S0(4:6)';

x = xout(:,1);
p = pout(:,1);

s3    = 1/(4-4^(1/3));
phi_l = [s3, s3, 1-4*s3, s3, s3];

for ii = 2:ntspan
    for jj = 1:length(phi_l)
        [x,p] = step(dt,x,p,partial_U,phi_l(jj),vargin);
    end
    xout(:,ii) = x;
    pout(:,ii) = p;
end
end

%% ---------------------------------------------------------------- one step
function [x,p] = step(dt,x,p,partial_U,phi_l,vargin)
dt    = phi_l*dt;
hdt   = dt/2;
dt2_4 = 4 + dt^2;

%% Stage 1, x <- x_n2
dx1 = dt*(2*p(1) + 2*x(2) + dt*p(2) - dt*x(1))/dt2_4;
dx2 = dt*(2*p(2) - 2*x(1) - dt*p(1) - dt*x(2))/dt2_4;
dx3 = hdt*p(3);

x(1) = x(1) + dx1;
x(2) = x(2) + dx2;
x(3) = x(3) + dx3;

%% Force at x_n2
dU = partial_U(x,vargin);

%% Stage 2, p <- p_n1
dp1 = dt*(-2*dt*p(1) + 4*p(2) - 4*dU(1) - 2*dt*dU(2))/dt2_4;
dp2 = dt*(-4*p(1) - 2*dt*p(2) + 2*dt*dU(1) - 4*dU(2))/dt2_4;
dp3 = -dt*dU(3);

p(1) = p(1) + dp1;
p(2) = p(2) + dp2;
p(3) = p(3) + dp3;

%% Stage 3, x <- x_n1
dxx1 = hdt*(x(2) + p(1));
dxx2 = hdt*(p(2) - x(1));
dxx3 = hdt*p(3);

x(1) = x(1) + dxx1;
x(2) = x(2) + dxx2;
x(3) = x(3) + dxx3;
end
