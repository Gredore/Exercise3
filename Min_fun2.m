function [h, ceq]= Min_fun2(in)

V = in(1);
T = in(2);
Steam_frac = in(3);


V = in(1);
T = in(2);
Steam_frac = in(3);


P = 2.4;
Fa0 = 10.4;

two_and_three = 1;
Vspan = [0 V];

Thermo_data(:,1) = [30 147 0 83 52 50 -75 -242];
Thermo_data(:,1) = Thermo_data(:,1) .*1000;
Thermo_data(:,2) = [6.3 12.8 28.8 -1.8 11.1 1.6 13.9 30.7];
Thermo_data(:,3) = [0.469 0.420 0.000 0.330 0.119 0.396 0.075 0.089];
Thermo_data(:,4) = [-15.7 -14.4 0.138 -11.4 -3.61 -13.3 -1.57 0.146];
Thermo_data(:,4) = Thermo_data(:,4) .*(10^(-5));

y0 = [Fa0 0 0 0 0 0 0 T]; %Initial conditions
    
Other_in = [
    Steam_frac*Fa0 %Steam,
    P,
];


opts = odeset('RelTol',1e-10,'AbsTol',1e-8);
[v, y] = ode45(@(V,Y) ODEfun_adiabatic(V, Y, Other_in, Thermo_data), Vspan, y0, opts);

Selec = y(end,2) ./ (y(end,2) + y(end,4) + y(end,6));
Conv = (Fa0  - y(end,1))./Fa0;

h = 0.95 - Selec;
ceq = 0;
end