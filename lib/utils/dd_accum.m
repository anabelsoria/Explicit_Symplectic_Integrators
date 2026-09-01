function [zh, zl] = dd_accum(zh, zl, d)
% Error-free accumulation of a plain-double increment d into a
% double-double pair (zh,zl): TwoSum captures the rounding
% error of zh+d, then a quickTwoSum renormalization folds it
% back in. Used by SI_EOM_ddInc in place of Kahan compensation
% (comp_sum's dd counterpart). Ported from the local `acc`
% helper in CHANCE's SI4_CR3BP_5stage_ddInc.m.
[s, e] = twoSum(zh, d);
zl = zl + e;
[zh, zl] = quickTwoSum(s, zl);
end
