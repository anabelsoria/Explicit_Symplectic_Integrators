function diff_g = find_G_heggie
clc; clear all
syms q [1 3] real
syms rM [1 3] real
syms rE [1 3] real


r12 = sqrt(sum((q - rE).^2));
r13 = sqrt(sum((q - rM).^2));
r23 = sqrt(sum((rM - rE).^2));

g = r12*r13*r23 / (r12 + r13 + r23)^(3/2);
% 
% diff_r12_q1 = simplify(diff(r12,q(1)));
% diff_r13_q1 = diff(r13,q(1));
% diff_r23_q1 = diff(r23,q(1));

diff_g_q1 = simplify(diff(g,q(1)))
diff_g_q2 = simplify(diff(g,q(2)))
diff_g_q3 = simplify(diff(g,q(3)))

diff_g = [diff_g_q1; diff_g_q2; diff_g_q3]';

% G = -1/g * [diff_g_q1; diff_g_q2; diff_g_q3]' * Hp_CR3BP_centerAtP2(mu_EM,p,q)';