% ========================================================================
% --- DOUBLE-DOUBLE ARITHMETIC UTILITIES ---
% ========================================================================

function dd = dd_from_double(x)
    dd = [x, 0.0];
end

function d = dd2double(dd)
    d = dd(1) + dd(2);
end

function dd = dd_add(a, b)
    [s, e1] = twoSum(a(1), b(1));
    e = a(2) + b(2) + e1;
    [hi, lo] = quickTwoSum(s, e);
    dd = [hi, lo];
end

function dd = dd_sub(a, b)
    [s, e1] = twoSum(a(1), -b(1));
    e = a(2) - b(2) + e1;
    [hi, lo] = quickTwoSum(s, e);
    dd = [hi, lo];
end

function dd = dd_mul(a, b)
    [p1, e1] = twoProd(a(1), b(1));
    t = a(1)*b(2) + a(2)*b(1);
    s = p1 + t;
    e = e1 + (t - (s - p1)) + a(2)*b(2);
    [hi, lo] = quickTwoSum(s, e);
    dd = [hi, lo];
end

function dd = dd_neg(a)
    dd = [-a(1), -a(2)];
end

function dd = dd_recip(d)
    r0 = 1.0 / d(1);
    r = dd_from_double(r0);
    one = dd_from_double(1.0);
    dr = dd_mul(d, r);
    tmp = dd_sub(one, dr);
    corr = dd_mul(r, tmp);
    r = dd_add(r, corr);
    [hi, lo] = quickTwoSum(r(1), r(2));
    dd = [hi, lo];
end

% ========================
% Error-free transforms
% ========================

function [s, e] = twoSum(a, b)
    s = a + b;
    bp = s - a;
    e = (a - (s - bp)) + (b - bp);
end

function [s, e] = quickTwoSum(a, b)
    s = a + b;
    e = b - (s - a);
end

function [p, e] = twoProd(a, b)
    c = 134217729 * a; % 2^27+1
    a_hi = c - (c - a);
    a_lo = a - a_hi;
    c = 134217729 * b;
    b_hi = c - (c - b);
    b_lo = b - b_hi;
    p = a * b;
    e = ((a_hi * b_hi - p) + a_hi*b_lo + a_lo*b_hi) + a_lo*b_lo;
end