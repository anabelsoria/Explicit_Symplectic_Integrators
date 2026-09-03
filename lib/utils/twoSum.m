function [s, e] = twoSum(a, b)
% Error-free transformation: s = fl(a+b), e = exact rounding error,
% so that a+b == s+e in infinite precision.
s = a + b;
bp = s - a;
e = (a - (s - bp)) + (b - bp);
end
