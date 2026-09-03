function [zh, zl] = dd_accum(zh, zl, d)
% Error-free accumulation of a plain-double d into dd pair (zh,zl).
% dd's counterpart to comp_sum.
[s, e] = twoSum(zh, d);
zl = zl + e;
[zh, zl] = quickTwoSum(s, zl);
end
