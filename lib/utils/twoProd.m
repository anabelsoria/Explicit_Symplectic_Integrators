function [p, e] = twoProd(a, b)
% Error-free transformation: p = fl(a*b), e = exact rounding error,
% so that a*b == p+e in infinite precision. Uses Dekker's splitting.
c = 134217729 * a; % 2^27+1
a_hi = c - (c - a);
a_lo = a - a_hi;
c = 134217729 * b;
b_hi = c - (c - b);
b_lo = b - b_hi;
p = a * b;
e = ((a_hi * b_hi - p) + a_hi*b_lo + a_lo*b_hi) + a_lo*b_lo;
end
