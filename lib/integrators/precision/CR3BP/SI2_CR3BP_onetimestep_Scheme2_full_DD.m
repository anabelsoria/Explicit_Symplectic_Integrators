function [x_n1_hi, x_n1_lo, p_n1_hi, p_n1_lo, Phi_hi, Phi_lo] = ...
    SI2_CR3BP_onetimestep_Scheme2_full_DD(dt, x_hi, x_lo, p_hi, p_lo, phi_l, vargin)
% One Scheme 2 step with BOTH state AND STM in double-double precision.
%
% State (x, p) carried as DD pairs (hi, lo); true value = hi + lo.
% STM  Phi = d[x_n1; p_n1]/d[x; p] also in DD.
% Requires dd_matmul.m on the path.
%
% Inputs:
%   dt              macro time step
%   x_hi, x_lo     3x1 position DD pair
%   p_hi, p_lo     3x1 canonical momentum DD pair
%   phi_l           Suzuki sub-step fraction
%   vargin          [GM1, GM2, x1, x2]
%
% Outputs:
%   x_n1_hi/lo, p_n1_hi/lo   updated state as DD pairs
%   Phi_hi, Phi_lo            6x6 per-step STM as DD pair

    % --- DD step size h = phi_l * dt ---
    [h_hi, h_lo] = tp(phi_l, dt);
    hw2_hi = h_hi / 2;
    hw2_lo = h_lo / 2;

    % --- DD matrices T, D ---
    [hw2sq_hi, hw2sq_lo] = ddm(hw2_hi, hw2_lo, hw2_hi, hw2_lo);
    [den_hi, den_lo] = dda(1, 0, hw2sq_hi, hw2sq_lo);

    [invd_hi, invd_lo] = ddd(1, 0, den_hi, den_lo);
    [hw2i_hi, hw2i_lo] = ddm(hw2_hi, hw2_lo, invd_hi, invd_lo);

    T_hi = [ invd_hi,  hw2i_hi, 0;
            -hw2i_hi,  invd_hi, 0;
             0,        0,       1];
    T_lo = [ invd_lo,  hw2i_lo, 0;
            -hw2i_lo,  invd_lo, 0;
             0,        0,       0];

    D_hi = [ 1,       hw2_hi, 0;
            -hw2_hi,  1,      0;
             0,       0,      1];
    D_lo = [ 0,       hw2_lo, 0;
            -hw2_lo,  0,      0;
             0,       0,      0];

    % ================================================================
    %  DD STATE UPDATE
    % ================================================================

    % Step 1: x_mid = T * (x + h/2 * p)
    [hp_hi, hp_lo] = dd_smat(hw2_hi, hw2_lo, p_hi, p_lo);
    [xhp_hi, xhp_lo] = dd_madd(x_hi, x_lo, hp_hi, hp_lo);
    [xm_hi, xm_lo] = dd_matmul(T_hi, T_lo, xhp_hi, xhp_lo);

    % Step 2: p_n1 = T * (D*p - h * dU(x_mid))
    [Dp_hi, Dp_lo] = dd_matmul(D_hi, D_lo, p_hi, p_lo);
    [dU_hi, dU_lo] = partial_U_dd(xm_hi, xm_lo, vargin);
    [hdU_hi, hdU_lo] = dd_smat(h_hi, h_lo, dU_hi, dU_lo);
    [inn_hi, inn_lo] = dd_msub(Dp_hi, Dp_lo, hdU_hi, hdU_lo);
    [p_n1_hi, p_n1_lo] = dd_matmul(T_hi, T_lo, inn_hi, inn_lo);

    % Step 3: x_n1 = D*x_mid + h/2 * p_n1
    [Dxm_hi, Dxm_lo] = dd_matmul(D_hi, D_lo, xm_hi, xm_lo);
    [hp1_hi, hp1_lo] = dd_smat(hw2_hi, hw2_lo, p_n1_hi, p_n1_lo);
    [x_n1_hi, x_n1_lo] = dd_madd(Dxm_hi, Dxm_lo, hp1_hi, hp1_lo);

    % ================================================================
    %  DD STM (same as SI2_..._STM_DD but Hessian uses DD midpoint)
    % ================================================================

    [Urr_hi, Urr_lo] = hessian_U_grav_dd(xm_hi, xm_lo, vargin);

    [TU_hi, TU_lo] = dd_matmul(T_hi, T_lo, Urr_hi, Urr_lo);
    [V_hi, V_lo]   = dd_matmul(TU_hi, TU_lo, T_hi, T_lo);

    [DT_hi, DT_lo] = dd_matmul(D_hi, D_lo, T_hi, T_lo);
    [TD_hi, TD_lo] = dd_matmul(T_hi, T_lo, D_hi, D_lo);

    [h2_hi, h2_lo]   = ddm(h_hi, h_lo, h_hi, h_lo);
    [h2h_hi, h2h_lo] = ddm(h2_hi, h2_lo, 0.5, 0);
    [hh_hi, hh_lo]   = ddm(h_hi, h_lo, 0.5, 0);
    [h3_hi, h3_lo]   = ddm(h2_hi, h2_lo, h_hi, h_lo);
    [h3q_hi, h3q_lo] = ddm(h3_hi, h3_lo, 0.25, 0);

    % Phi_21 = -h * V
    [P21_hi, P21_lo] = dd_smat(-h_hi, -h_lo, V_hi, V_lo);

    % (h^2/2)*V (shared by Phi_11 and Phi_22)
    [sV_hi, sV_lo] = dd_smat(h2h_hi, h2h_lo, V_hi, V_lo);

    % Phi_11 = DT - (h^2/2)*V
    [P11_hi, P11_lo] = dd_msub(DT_hi, DT_lo, sV_hi, sV_lo);

    % Phi_22 = TD - (h^2/2)*V
    [P22_hi, P22_lo] = dd_msub(TD_hi, TD_lo, sV_hi, sV_lo);

    % Phi_12 = (h/2)*(DT + TD) - (h^3/4)*V
    [sum_hi, sum_lo] = dd_madd(DT_hi, DT_lo, TD_hi, TD_lo);
    [t1_hi, t1_lo]   = dd_smat(hh_hi, hh_lo, sum_hi, sum_lo);
    [t2_hi, t2_lo]   = dd_smat(h3q_hi, h3q_lo, V_hi, V_lo);
    [P12_hi, P12_lo] = dd_msub(t1_hi, t1_lo, t2_hi, t2_lo);

    Phi_hi = [P11_hi, P12_hi; P21_hi, P22_hi];
    Phi_lo = [P11_lo, P12_lo; P21_lo, P22_lo];
end


% =====================================================================
%  DD gravitational gradient  (takes DD position)
% =====================================================================

function [dU_hi, dU_lo] = partial_U_dd(x_hi, x_lo, vargin)
    GM1 = vargin(1); GM2 = vargin(2);
    x1 = vargin(3);  x2 = vargin(4);

    % DD displacements from each body
    [dx1_hi, dx1_lo] = dda(x_hi(1), x_lo(1), -x1, 0);
    [dx2_hi, dx2_lo] = dda(x_hi(1), x_lo(1), -x2, 0);
    yh = x_hi(2); yl = x_lo(2);
    zh = x_hi(3); zl = x_lo(3);

    % R1^2 = dx1^2 + y^2 + z^2
    [a2h, a2l] = ddm(dx1_hi, dx1_lo, dx1_hi, dx1_lo);
    [b2h, b2l] = ddm(yh, yl, yh, yl);
    [c2h, c2l] = ddm(zh, zl, zh, zl);
    [R1sq_hi, R1sq_lo] = dda(a2h, a2l, b2h, b2l);
    [R1sq_hi, R1sq_lo] = dda(R1sq_hi, R1sq_lo, c2h, c2l);

    % R2^2
    [a2h, a2l] = ddm(dx2_hi, dx2_lo, dx2_hi, dx2_lo);
    [R2sq_hi, R2sq_lo] = dda(a2h, a2l, b2h, b2l);
    [R2sq_hi, R2sq_lo] = dda(R2sq_hi, R2sq_lo, c2h, c2l);

    % R1, R2
    [R1_hi, R1_lo] = ddsqrt(R1sq_hi, R1sq_lo);
    [R2_hi, R2_lo] = ddsqrt(R2sq_hi, R2sq_lo);

    % R^3 = R * R^2
    [R1c_hi, R1c_lo] = ddm(R1_hi, R1_lo, R1sq_hi, R1sq_lo);
    [R2c_hi, R2c_lo] = ddm(R2_hi, R2_lo, R2sq_hi, R2sq_lo);

    % GM / R^3
    [a1_hi, a1_lo] = ddd(GM1, 0, R1c_hi, R1c_lo);
    [a2_hi, a2_lo] = ddd(GM2, 0, R2c_hi, R2c_lo);

    % dU(i) = a1 * r1(i) + a2 * r2(i)
    r1_hi = [dx1_hi; yh; zh];  r1_lo = [dx1_lo; yl; zl];
    r2_hi = [dx2_hi; yh; zh];  r2_lo = [dx2_lo; yl; zl];

    dU_hi = zeros(3,1); dU_lo = zeros(3,1);
    for i = 1:3
        [t1h, t1l] = ddm(a1_hi, a1_lo, r1_hi(i), r1_lo(i));
        [t2h, t2l] = ddm(a2_hi, a2_lo, r2_hi(i), r2_lo(i));
        [dU_hi(i), dU_lo(i)] = dda(t1h, t1l, t2h, t2l);
    end
end


% =====================================================================
%  DD Hessian of gravitational potential  (takes DD position)
% =====================================================================

function [Urr_hi, Urr_lo] = hessian_U_grav_dd(x_hi, x_lo, vargin)
    GM1 = vargin(1); GM2 = vargin(2);
    x1 = vargin(3);  x2 = vargin(4);

    [dx1_hi, dx1_lo] = dda(x_hi(1), x_lo(1), -x1, 0);
    [dx2_hi, dx2_lo] = dda(x_hi(1), x_lo(1), -x2, 0);
    yh = x_hi(2); yl = x_lo(2);
    zh = x_hi(3); zl = x_lo(3);

    r1_hi = [dx1_hi; yh; zh];  r1_lo = [dx1_lo; yl; zl];
    r2_hi = [dx2_hi; yh; zh];  r2_lo = [dx2_lo; yl; zl];

    [a2h, a2l] = ddm(dx1_hi, dx1_lo, dx1_hi, dx1_lo);
    [b2h, b2l] = ddm(yh, yl, yh, yl);
    [c2h, c2l] = ddm(zh, zl, zh, zl);
    [R1sq_hi, R1sq_lo] = dda(a2h, a2l, b2h, b2l);
    [R1sq_hi, R1sq_lo] = dda(R1sq_hi, R1sq_lo, c2h, c2l);

    [a2h, a2l] = ddm(dx2_hi, dx2_lo, dx2_hi, dx2_lo);
    [R2sq_hi, R2sq_lo] = dda(a2h, a2l, b2h, b2l);
    [R2sq_hi, R2sq_lo] = dda(R2sq_hi, R2sq_lo, c2h, c2l);

    [R1_hi, R1_lo] = ddsqrt(R1sq_hi, R1sq_lo);
    [R2_hi, R2_lo] = ddsqrt(R2sq_hi, R2sq_lo);

    [R1c_hi, R1c_lo] = ddm(R1_hi, R1_lo, R1sq_hi, R1sq_lo);   % R^3
    [R2c_hi, R2c_lo] = ddm(R2_hi, R2_lo, R2sq_hi, R2sq_lo);
    [R1f_hi, R1f_lo] = ddm(R1c_hi, R1c_lo, R1sq_hi, R1sq_lo); % R^5
    [R2f_hi, R2f_lo] = ddm(R2c_hi, R2c_lo, R2sq_hi, R2sq_lo);

    [a1_hi, a1_lo] = ddd(GM1, 0, R1c_hi, R1c_lo);   % GM/R^3
    [a2_hi, a2_lo] = ddd(GM2, 0, R2c_hi, R2c_lo);

    [g1_hi, g1_lo] = tp(3, GM1);
    [b1_hi, b1_lo] = ddd(g1_hi, g1_lo, R1f_hi, R1f_lo);  % 3*GM/R^5
    [g2_hi, g2_lo] = tp(3, GM2);
    [b2_hi, b2_lo] = ddd(g2_hi, g2_lo, R2f_hi, R2f_lo);

    Urr_hi = zeros(3); Urr_lo = zeros(3);

    for i = 1:3
        for j = i:3
            if i == j
                [d_hi, d_lo] = dda(a1_hi, a1_lo, a2_hi, a2_lo);
            else
                d_hi = 0; d_lo = 0;
            end

            [p1_hi, p1_lo] = ddm(r1_hi(i), r1_lo(i), r1_hi(j), r1_lo(j));
            [bp1_hi, bp1_lo] = ddm(b1_hi, b1_lo, p1_hi, p1_lo);

            [p2_hi, p2_lo] = ddm(r2_hi(i), r2_lo(i), r2_hi(j), r2_lo(j));
            [bp2_hi, bp2_lo] = ddm(b2_hi, b2_lo, p2_hi, p2_lo);

            [e_hi, e_lo] = dda(d_hi, d_lo, -bp1_hi, -bp1_lo);
            [e_hi, e_lo] = dda(e_hi, e_lo, -bp2_hi, -bp2_lo);

            Urr_hi(i,j) = e_hi; Urr_lo(i,j) = e_lo;
            if i ~= j
                Urr_hi(j,i) = e_hi; Urr_lo(j,i) = e_lo;
            end
        end
    end
end


% =====================================================================
%  DD scalar primitives
% =====================================================================

function [s, e] = ts(a, b)
    s = a + b; v = s - a; e = (a - (s - v)) + (b - v);
end

function [p, e] = tp(a, b)
    p = a * b;
    c = 134217729 * a; ahi = c - (c - a); alo = a - ahi;
    c = 134217729 * b; bhi = c - (c - b); blo = b - bhi;
    e = ((ahi*bhi - p) + ahi*blo + alo*bhi) + alo*blo;
end

function [sh, sl] = dda(ah, al, bh, bl)
    [sh, eh] = ts(ah, bh);
    sl = eh + al + bl;
    [sh, sl] = ts(sh, sl);
end

function [sh, sl] = ddm(ah, al, bh, bl)
    [ph, pl] = tp(ah, bh);
    pl = pl + ah*bl + al*bh;
    [sh, sl] = ts(ph, pl);
end

function [qh, ql] = ddd(ah, al, bh, bl)
    qh = ah / bh;
    [ph, pl] = tp(qh, bh);
    ql = ((ah - ph) - pl + al - qh*bl) / bh;
    [qh, ql] = ts(qh, ql);
end

function [sh, sl] = ddsqrt(ah, al)
    sh = sqrt(ah);
    [ph, pl] = tp(sh, sh);
    sl = ((ah - ph) - pl + al) / (2*sh);
    [sh, sl] = ts(sh, sl);
end


% =====================================================================
%  DD matrix helpers (element-wise vectorized)
% =====================================================================

function [Ch, Cl] = dd_smat(sh, sl, Mh, Ml)
    fac = 134217729;
    P = sh * Mh;
    c = fac * sh; ahi = c - (c - sh); alo = sh - ahi;
    c = fac * Mh; bhi = c - (c - Mh); blo = Mh - bhi;
    E = ((ahi*bhi - P) + ahi*blo + alo*bhi) + alo*blo;
    pl = E + sh*Ml + sl*Mh;
    S = P + pl; V = S - P;
    Cl = (P - (S - V)) + (pl - V);
    Ch = S;
end

function [Ch, Cl] = dd_madd(Ah, Al, Bh, Bl)
    S = Ah + Bh; V = S - Ah;
    E = (Ah - (S - V)) + (Bh - V);
    sl = E + Al + Bl;
    Ch = S + sl; V = Ch - S;
    Cl = (S - (Ch - V)) + (sl - V);
end

function [Ch, Cl] = dd_msub(Ah, Al, Bh, Bl)
    [Ch, Cl] = dd_madd(Ah, Al, -Bh, -Bl);
end
