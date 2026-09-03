function [s, e] = quickTwoSum(a, b)
% Error-free transformation, cheaper than twoSum but requires |a| >= |b|.
s = a + b;
e = b - (s - a);
end
