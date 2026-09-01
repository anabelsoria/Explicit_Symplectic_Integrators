function [xout,pout,xlout,plout] = SI4_CR3BP_5stage_ddInc(t0,tf,dt,S0,...
    w,partial_U,vargin)
% Five-stage SI4 with the state carried as a double-double pair and the
% increments evaluated in ordinary double precision.
%
% This is the pattern of sympl.csi_step_inc. Each increment is derived
% algebraically so that it is of order dt and no large quantities cancel, it
% is evaluated in double at the compensated point xh+xl, and it is then
% accumulated into the pair error-free with TwoSum followed by a Fast2Sum
% renormalization.
%
% It differs from SI4_CR3BP_5stage_dd, which propagates the full map through
% double-double multiplication and division. Here only the accumulation is
% done in extended precision, which is where the error that grows with the
% step count enters.
%
% The low words are returned as well, since the round-off being measured can
% fall below the resolution of the double xout on its own.

if ~exist('vargin','var')
    vargin = [];
end

tspan  = t0:dt:tf;
ntspan = length(tspan);

xout  = nan(3,ntspan);  pout  = nan(3,ntspan);
xlout = nan(3,ntspan);  plout = nan(3,ntspan);

xh = S0(1:3);  xl = zeros(3,1);
ph = S0(4:6);  pl = zeros(3,1);

xout(:,1) = xh;  xlout(:,1) = xl;
pout(:,1) = ph;  plout(:,1) = pl;

s3    = 1/(4-4^(1/3));
phi_l = [s3, s3, 1-4*s3, s3, s3];

for ii = 2:ntspan
    for jj = 1:length(phi_l)
        [xh,xl,ph,pl] = step(dt,xh,xl,ph,pl,partial_U,phi_l(jj),vargin);
    end
    xout(:,ii) = xh;  xlout(:,ii) = xl;
    pout(:,ii) = ph;  plout(:,ii) = pl;
end
end

%% ---------------------------------------------------------------- one step
function [xh,xl,ph,pl] = step(dt,xh,xl,ph,pl,partial_U,phi_l,vargin)
dt  = phi_l*dt;
hdt = dt/2;
dt2_4 = 4 + dt^2;

% Increments are evaluated at the compensated point.
x = xh + xl;
p = ph + pl;

%% Stage 1, x <- x_n2
dx1 = dt*(2*p(1) + 2*x(2) + dt*p(2) - dt*x(1))/dt2_4;
dx2 = dt*(2*p(2) - 2*x(1) - dt*p(1) - dt*x(2))/dt2_4;
dx3 = hdt*p(3);

[xh(1),xl(1)] = acc(xh(1),xl(1),dx1);
[xh(2),xl(2)] = acc(xh(2),xl(2),dx2);
[xh(3),xl(3)] = acc(xh(3),xl(3),dx3);

x = xh + xl;

%% Force at x_n2
dU = partial_U(x,vargin);

%% Stage 2, p <- p_n1
dp1 = dt*(-2*dt*p(1) + 4*p(2) - 4*dU(1) - 2*dt*dU(2))/dt2_4;
dp2 = dt*(-4*p(1) - 2*dt*p(2) + 2*dt*dU(1) - 4*dU(2))/dt2_4;
dp3 = -dt*dU(3);

[ph(1),pl(1)] = acc(ph(1),pl(1),dp1);
[ph(2),pl(2)] = acc(ph(2),pl(2),dp2);
[ph(3),pl(3)] = acc(ph(3),pl(3),dp3);

p = ph + pl;

%% Stage 3, x <- x_n1
dxx1 = hdt*(x(2) + p(1));
dxx2 = hdt*(p(2) - x(1));
dxx3 = hdt*p(3);

[xh(1),xl(1)] = acc(xh(1),xl(1),dxx1);
[xh(2),xl(2)] = acc(xh(2),xl(2),dxx2);
[xh(3),xl(3)] = acc(xh(3),xl(3),dxx3);
end

%% ------------------------------------------------- error-free accumulation
function [zh, zl] = acc(zh, zl, d)
s  = zh + d;
v  = s - zh;
e  = (zh - (s - v)) + (d - v);      % exact rounding error of zh + d
zl = zl + e;
t  = s + zl;
zl = zl - (t - s);
zh = t;
end
