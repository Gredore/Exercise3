function f = ODEfun_adiabatic(V, Y, Other_in, Thermo_data)

Fa = Y(1);
Fb = Y(2);
Fc = Y(3);
Fd = Y(4);
Fe = Y(5);
Ff = Y(6);
Fg = Y(7);
T  = Y(8);

K_con = [-13.1 -13050 4.27 -3.0*10^(-3) 3.0*10^(-7)];

const1 = exp(13.4 - (10930./T));
const2 = exp(23.9 - (25000./T));
const3 = exp(14.3 - (12000./T));
K = exp(K_con(1) + (K_con(2)./T) + (K_con(3).*log(T)) +(K_con(4) .*T) + (K_con(5) .* (T.^2)));

Fh = Other_in(1);
P = Other_in(2);

Pa = P .* Fa ./ (Fa + Fb + Fc + Fd + Fe + Ff + Fg + Fh) ;
Pb = P .* Fb ./ (Fa + Fb + Fc + Fd + Fe + Ff + Fg + Fh) ;
Pc = P .* Fc ./ (Fa + Fb + Fc + Fd + Fe + Ff + Fg + Fh) ;




r1 = const1.*(Pa - (Pb.*Pc)/K);
r2 = const2.*Pa;
r3 = const3.*Pa.*Pc;



Cp = Thermo_data(:,2) + (Thermo_data(:,3).*T) + (Thermo_data(:,4).*(T^2));

deltaHa = Thermo_data(1,1) + integral(@(T) Thermo_data(1,2) + (Thermo_data(1,3).*T) + (Thermo_data(1,4).*(T.^2)),298,T,'RelTol', 1e-4, 'AbsTol', 1e-6);
deltaHb = Thermo_data(2,1) + integral(@(T) Thermo_data(2,2) + (Thermo_data(2,3).*T) + (Thermo_data(2,4).*(T.^2)),298,T,'RelTol', 1e-4, 'AbsTol', 1e-6);
deltaHc = Thermo_data(3,1) + integral(@(T) Thermo_data(3,2) + (Thermo_data(3,3).*T) + (Thermo_data(3,4).*(T.^2)),298,T,'RelTol', 1e-4, 'AbsTol', 1e-6);
deltaHd = Thermo_data(4,1) + integral(@(T) Thermo_data(4,2) + (Thermo_data(4,3).*T) + (Thermo_data(4,4).*(T.^2)),298,T,'RelTol', 1e-4, 'AbsTol', 1e-6);
deltaHe = Thermo_data(5,1) + integral(@(T) Thermo_data(5,2) + (Thermo_data(5,3).*T) + (Thermo_data(5,4).*(T.^2)),298,T,'RelTol', 1e-4, 'AbsTol', 1e-6);
deltaHf = Thermo_data(6,1) + integral(@(T) Thermo_data(6,2) + (Thermo_data(6,3).*T) + (Thermo_data(6,4).*(T.^2)),298,T,'RelTol', 1e-4, 'AbsTol', 1e-6);
deltaHg = Thermo_data(7,1) + integral(@(T) Thermo_data(7,2) + (Thermo_data(7,3).*T) + (Thermo_data(7,4).*(T.^2)),298,T,'RelTol', 1e-4, 'AbsTol', 1e-6);


dFadV = -r1 - r2 - r3;
dFbdV = r1;
dFcdV = r1 - r3;
dFddV = r2;
dFedV = r2;
dFfdV = r3;
dFgdV = r3;

dFadV_H = dFadV .* deltaHa;
dFbdV_H = dFbdV .* deltaHb;
dFcdV_H = dFcdV .* deltaHc;
dFddV_H = dFddV .* deltaHd;
dFedV_H = dFedV .* deltaHe;
dFfdV_H = dFfdV .* deltaHf;
dFgdV_H = dFgdV .* deltaHg;


dTdV = - (dFadV_H + dFbdV_H + dFcdV_H + dFddV_H + dFedV_H + dFfdV_H + dFgdV_H ) ./ ...
         ((Fa.*Cp(1)) + (Fb.*Cp(2)) + (Fc.*Cp(3)) + (Fd.*Cp(4)) + (Fe.*Cp(5)) + (Ff.*Cp(6)) + (Fg.*Cp(7)) + (Fh.*Cp(8)));

f = [dFadV; dFbdV; dFcdV; dFddV; dFedV; dFfdV; dFgdV; dTdV];
end
