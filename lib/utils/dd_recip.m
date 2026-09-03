function dd = dd_recip(d)
% Double-double reciprocal: 1/d, via one Newton-Raphson correction step.
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
