function [q,p,e_q,e_p] = SI2_ER3BP_onetimestep_CompSumScalar(dt,q,p,U,partial_U,phi_l,vargin,e_q,e_p)
% ER3BP counterpart of SI2_CR3BP_onetimestep_Scheme2_CompSumScalar.
%
% The update of SI2_ER3BP_onetimestep is rearranged so that each stage is
% evaluated as an analytic increment and accumulated by compensated
% summation, as described in Section 4.4. Results agree with the closed form
% to round-off; only the accumulated error differs.
%
% The state is eight dimensional, q = (x,y,z,q4) and p = (px,py,pz,nu). The
% spatial part has the same structure as the CR3BP scheme of Scheme 2, with
% the same T and D, so the spatial increments carry over unchanged except
% that the single force term h*dU is replaced by
%
%   G = (h/2) [ (e cos(nu) q_n2 + dU)/(1 + e cos(nu))
%             + (e cos(nu') q_n2 + dU)/(1 + e cos(nu')) ]
%
% evaluated before and after the true anomaly is advanced. The fourth
% components are already written as increments in the closed form and need no
% rearrangement.
%
% The gradient is evaluated once and reused for both terms of G. The closed
% form calls partial_U twice with the same argument, which is redundant.

mu = vargin(2);
ec = vargin(end);

h   = phi_l*dt;
hd  = h/2;
den = 1 + hd*hd;
d4  = h*h + 4;

%% Stage 1, q(1:3) advances to q_n2. Increments formed before any update.
dq1 = h*( 2*p(1) + 2*q(2) + h*p(2) - h*q(1) )/d4;
dq2 = h*( 2*p(2) - 2*q(1) - h*p(1) - h*q(2) )/d4;
dq3 = hd*p(3);

[q(1), e_q(1)] = comp_sum(q(1),e_q(1),dq1);
[q(2), e_q(2)] = comp_sum(q(2),e_q(2),dq2);
[q(3), e_q(3)] = comp_sum(q(3),e_q(3),dq3);

% Values at q_n2, needed again in stage 3 before q(1:3) advances further.
n1 = q(1);  n2 = q(2);  n3 = q(3);
qn2 = [n1;n2;n3];

W  = n1*n1 + n2*n2 + n3*n3 - 2*U(qn2,mu);
dU = partial_U(qn2,vargin);

%% Stage 1b, the coordinate conjugate to the true anomaly
nu  = p(4);
c   = cos(nu);
dq4 = -h*ec*sin(nu)/(4*(1 + ec*c)^2) * W;
[q(4), e_q(4)] = comp_sum(q(4),e_q(4),dq4);

%% Stage 2, the true anomaly advances, then the momenta
[p(4), e_p(4)] = comp_sum(p(4),e_p(4),h);
nu2 = p(4);
c2  = cos(nu2);

G = hd*( (ec*c *qn2 + dU)/(1 + ec*c ) ...
       + (ec*c2*qn2 + dU)/(1 + ec*c2) );

dp1 = ( -h*hd*p(1) + h*p(2) - G(1) - hd*G(2) )/den;
dp2 = ( -h*p(1) - h*hd*p(2) + hd*G(1) - G(2) )/den;
dp3 = -G(3);

[p(1), e_p(1)] = comp_sum(p(1),e_p(1),dp1);
[p(2), e_p(2)] = comp_sum(p(2),e_p(2),dp2);
[p(3), e_p(3)] = comp_sum(p(3),e_p(3),dp3);

%% Stage 3, q(1:3) advances to q_n1 using q_n2 and the updated momenta
dr1 = hd*( n2 + p(1) );
dr2 = hd*( p(2) - n1 );
dr3 = hd*p(3);

[q(1), e_q(1)] = comp_sum(q(1),e_q(1),dr1);
[q(2), e_q(2)] = comp_sum(q(2),e_q(2),dr2);
[q(3), e_q(3)] = comp_sum(q(3),e_q(3),dr3);

%% Stage 3b
dq4b = -hd*(ec/2)*sin(nu2)/(1 + ec*c2)^2 * W;
[q(4), e_q(4)] = comp_sum(q(4),e_q(4),dq4b);
end
