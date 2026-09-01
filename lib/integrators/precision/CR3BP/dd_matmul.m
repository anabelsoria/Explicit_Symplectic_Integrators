function [Ch, Cl] = dd_matmul(Ah, Al, Bh, Bl)
% DD_MATMUL  Double-double precision matrix multiply: C = A * B
%
%   A = (Ah + Al),  B = (Bh + Bl),  C = (Ch + Cl)
%
% Each matrix is stored as a pair (hi, lo) of double-precision matrices.
% The true value is hi + lo.  This gives ~32 digits of effective precision
% for the accumulated product.
%
% Fully vectorized over matrix dimensions — only loops over the shared
% inner dimension K.

    [n, K] = size(Ah);
    m      = size(Bh, 2);
    Ch     = zeros(n, m);
    Cl     = zeros(n, m);

    fac = 134217729;   % 2^27 + 1 (for Split)

    for k = 1:K
        ah = Ah(:, k);   al = Al(:, k);
        bh = Bh(k, :);   bl = Bl(k, :);

        % TwoProd(ah, bh): error-free outer product
        P = ah * bh;

        c   = fac * ah;
        ahi = c - (c - ah);
        alo = ah - ahi;

        c   = fac * bh;
        bhi = c - (c - bh);
        blo = bh - bhi;

        E = ((ahi * bhi - P) + ahi * blo + alo * bhi) + alo * blo;

        pl = E + ah * bl + al * bh;

        % TwoSum(P, pl)
        S = P + pl;
        V = S - P;
        pl = (P - (S - V)) + (pl - V);
        ph = S;

        % dd_add: accumulate (Ch, Cl) += (ph, pl)
        S  = Ch + ph;
        V  = S - Ch;
        sl = (Ch - (S - V)) + (ph - V) + (Cl + pl);

        Ch = S + sl;
        V  = Ch - S;
        Cl = (S - (Ch - V)) + (sl - V);
    end
end
