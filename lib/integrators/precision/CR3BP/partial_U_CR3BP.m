function dU = partial_U_CR3BP(s,vargin)

GM1 = vargin(1);
GM2 = vargin(2);
x1  = vargin(3);
x2  = vargin(4);

if length(s) > 2
    x = s(1);
    y = s(2);
    z = s(3);
else
    x = s(1);
    y = s(2);
    z = 0;
end

R1 = sqrt((x-x1)^2 + y^2 + z^2);      % R1: distance to m1, LARGER MASS
R2 = sqrt((x-x2)^2 + y^2 + z^2);      % R2: distance to m2, smaller mass

Ux =   GM1*(x-x1)/R1^(3) + GM2*(x-x2)/R2^(3) ;
Uy =   GM1* y    /R1^(3) + GM2* y    /R2^(3) ;
Uz =   GM1* z    /R1^3   + GM2* z    /R2^3;

if length(s) > 2
    dU = [Ux; Uy; Uz];
else
    dU = [Ux; Uy];
end
