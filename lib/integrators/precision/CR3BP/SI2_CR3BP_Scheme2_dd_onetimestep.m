function [xh,xl,ph,pl] = SI2_CR3BP_Scheme2_dd_onetimestep(dt,xh,xl,ph,pl,...
    phi_l,vargin)

dt = phi_l*dt;
hdt = dt/2;

GM1 = vargin(1);
GM2 = vargin(2);
x1  = vargin(3);
x2  = vargin(4);

% den = 1 + hdt^2  (w=1)
[hdt2h, hdt2l] = TwoProd(hdt, hdt);
[denh, denl] = dd_add(1, 0, hdt2h, hdt2l);

% inv_den = 1/den
[invdh, invdl] = dd_div(1, 0, denh, denl);

% Tinv coefficients: 1/den and hdt/den
[hdtdenh, hdtdenl] = dd_mul(hdt, 0, invdh, invdl);

%% v = x + hdt*p  (intermediate vector)
[t1h, t1l] = dd_mul(hdt, 0, ph(1), pl(1));
[v1h, v1l] = dd_add(xh(1), xl(1), t1h, t1l);

[t2h, t2l] = dd_mul(hdt, 0, ph(2), pl(2));
[v2h, v2l] = dd_add(xh(2), xl(2), t2h, t2l);

[t3h, t3l] = dd_mul(hdt, 0, ph(3), pl(3));
[v3h, v3l] = dd_add(xh(3), xl(3), t3h, t3l);

%% x_n2 = Tinv * v
% x_n2(1) = v1/den + hdt*v2/den
[a1h, a1l] = dd_mul(invdh, invdl, v1h, v1l);
[a2h, a2l] = dd_mul(hdtdenh, hdtdenl, v2h, v2l);
[xn2h1, xn2l1] = dd_add(a1h, a1l, a2h, a2l);

% x_n2(2) = -hdt*v1/den + v2/den
[a3h, a3l] = dd_mul(-hdtdenh, -hdtdenl, v1h, v1l);
[a4h, a4l] = dd_mul(invdh, invdl, v2h, v2l);
[xn2h2, xn2l2] = dd_add(a3h, a3l, a4h, a4l);

% x_n2(3) = v3
xn2h3 = v3h; xn2l3 = v3l;

%% Force evaluation at x_n2 (double precision)
xn2_1 = xn2h1; xn2_2 = xn2h2; xn2_3 = xn2h3;
r_norm = ((x1 - xn2_1)^2 + xn2_2^2 + xn2_3^2)^(1/2);
d_norm = ((x2 - xn2_1)^2 + xn2_2^2 + xn2_3^2)^(1/2);
r3 = r_norm^3;
d3 = d_norm^3;

dU1 = GM1*(xn2_1 - x1)/r3 + GM2*(xn2_1 - x2)/d3;
dU2 = GM1*xn2_2/r3 + GM2*xn2_2/d3;
dU3 = GM1*xn2_3/r3 + GM2*xn2_3/d3;

%% w = D*p - dt*dU
% D*p: w1 = p1 + hdt*p2, w2 = -hdt*p1 + p2, w3 = p3
[d1h, d1l] = dd_mul(hdt, 0, ph(2), pl(2));
[w1h, w1l] = dd_add(ph(1), pl(1), d1h, d1l);

[d2h, d2l] = dd_mul(-hdt, 0, ph(1), pl(1));
[w2h, w2l] = dd_add(d2h, d2l, ph(2), pl(2));

w3h = ph(3); w3l = pl(3);

% w = D*p - dt*dU
[f1h, f1l] = TwoProd(-dt, dU1);
[w1h, w1l] = dd_add(w1h, w1l, f1h, f1l);

[f2h, f2l] = TwoProd(-dt, dU2);
[w2h, w2l] = dd_add(w2h, w2l, f2h, f2l);

[f3h, f3l] = TwoProd(-dt, dU3);
[w3h, w3l] = dd_add(w3h, w3l, f3h, f3l);

%% p_n1 = Tinv * w
[b1h, b1l] = dd_mul(invdh, invdl, w1h, w1l);
[b2h, b2l] = dd_mul(hdtdenh, hdtdenl, w2h, w2l);
[pn1h1, pn1l1] = dd_add(b1h, b1l, b2h, b2l);

[b3h, b3l] = dd_mul(-hdtdenh, -hdtdenl, w1h, w1l);
[b4h, b4l] = dd_mul(invdh, invdl, w2h, w2l);
[pn1h2, pn1l2] = dd_add(b3h, b3l, b4h, b4l);

pn1h3 = w3h; pn1l3 = w3l;

%% x_n1 = D*x_n2 + hdt*p_n1
% D*x_n2: g1 = xn2(1) + hdt*xn2(2), g2 = -hdt*xn2(1) + xn2(2)
[e1h, e1l] = dd_mul(hdt, 0, xn2h2, xn2l2);
[g1h, g1l] = dd_add(xn2h1, xn2l1, e1h, e1l);

[e2h, e2l] = dd_mul(-hdt, 0, xn2h1, xn2l1);
[g2h, g2l] = dd_add(e2h, e2l, xn2h2, xn2l2);

% + hdt*p_n1
[h1h, h1l] = dd_mul(hdt, 0, pn1h1, pn1l1);
[xn1h1, xn1l1] = dd_add(g1h, g1l, h1h, h1l);

[h2h, h2l] = dd_mul(hdt, 0, pn1h2, pn1l2);
[xn1h2, xn1l2] = dd_add(g2h, g2l, h2h, h2l);

[h3h, h3l] = dd_mul(hdt, 0, pn1h3, pn1l3);
[xn1h3, xn1l3] = dd_add(xn2h3, xn2l3, h3h, h3l);

%% Output
xh = [xn1h1; xn1h2; xn1h3];
xl = [xn1l1; xn1l2; xn1l3];
ph = [pn1h1; pn1h2; pn1h3];
pl = [pn1l1; pn1l2; pn1l3];

end
