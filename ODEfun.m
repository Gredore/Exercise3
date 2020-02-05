function f = ODEfun(V, Y, Other_in)

Fa = Y(1);
Fb = Y(2);
Fc = Y(3);
Fd = Y(4);
Fe = Y(5);
Ff = Y(6);
Fg = Y(7);

const1 = Other_in(1);
const2 = Other_in(2);
const3 = Other_in(3);
K = Other_in(4);
Fh = Other_in(5);
P = Other_in(6);
two_and_three = Other_in(7);

Pa = P .* Fa ./ (Fa + Fb + Fc + Fd + Fe + Ff + Fg + Fh) ;
Pb = P .* Fb ./ (Fa + Fb + Fc + Fd + Fe + Ff + Fg + Fh) ;
Pc = P .* Fc ./ (Fa + Fb + Fc + Fd + Fe + Ff + Fg + Fh) ;




r1 = const1.*(Pa - (Pb.*Pc)/K);
r2 = const2.*Pa.*(two_and_three);
r3 = const3.*Pa.*Pc.*(two_and_three);


dFadV = -r1 - r2 - r3;
dFbdV = r1;
dFcdV = r1 - r3;
dFddV = r2;
dFedV = r2;
dFfdV = r3;
dFgdV = r3;

f = [dFadV; dFbdV; dFcdV; dFddV; dFedV; dFfdV; dFgdV];
end
