function [y_n1, e] = comp_sum(a,e,delta)

e = e + delta;

y_n1 = a + e;

e = e + (a-y_n1);