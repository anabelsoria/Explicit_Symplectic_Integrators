function dd = dd_add(a, b)
% Double-double addition: a + b, both [hi, lo] pairs.
[s, e1] = twoSum(a(1), b(1));
e = a(2) + b(2) + e1;
[hi, lo] = quickTwoSum(s, e);
dd = [hi, lo];
end
