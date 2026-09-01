function diff_g = find_G_russell_HR4BP
clc; clear all

syms q [1 3] real
syms rM [1 3] real
syms rE [1 3] real

syms A_E B_E C_E rH_E A_M B_M C_M rH_M alpha_E alpha_M real

dE = norm(q-rE);
dM = norm(q-rM);

rho_1 = rho_fun(A_E,B_E,C_E,dE,rH_E);
rho_2 = rho_fun(A_M,B_M,C_M,dM,rH_M);

g = rho_1^alpha_E * rho_2^alpha_M;

diff_g_q1 = simplify(diff(g,q(1)));
diff_g_q2 = simplify(diff(g,q(2)));
diff_g_q3 = simplify(diff(g,q(3)));

diff_g = [diff_g_q1; diff_g_q2; diff_g_q3]';

end

function rho = rho_fun(A,B,C,r,rH)

rho = (A^2 + A) / ( A + (A*(1-C)/(C+A))^(r/(B*rH)) ) - A;

end