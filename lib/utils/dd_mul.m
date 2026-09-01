function dd = dd_mul(a, b)
% Double-double multiplication: a * b, both [hi, lo] pairs.
[p1, e1] = twoProd(a(1), b(1));
t = a(1)*b(2) + a(2)*b(1);
s = p1 + t;
e = e1 + (t - (s - p1)) + a(2)*b(2);
[hi, lo] = quickTwoSum(s, e);
dd = [hi, lo];
end
