clear all
clf
close all

P = 2.4;
Fa0 = 10.4;
Steam_frac = 12;
two_and_three = 1;
Vspan = [0 10];
y0 = [Fa0 0 0 0 0 0 0];

K_con = [-13.1 -13050 4.27 -3.0*10^(-3) 3.0*10^(-7)];

counter = 1;

Ts = linspace(800, 1100, 100);

for T = Ts

Other_in = [
    exp(13.4 - (10930./T)),
    exp(23.9 - (25000./T)),
    exp(14.3 - (12000./T)),
    exp(K_con(1) + (K_con(2)./T) + (K_con(3).*log(T)) +(K_con(4) .*T) + (K_con(5) .* (T.^2))),
    Steam_frac*Fa0 %Steam,
    P,
    two_and_three % enable or disable reactions 2 and 3
];

opts = odeset('RelTol',1e-10,'AbsTol',1e-8);
[v, y] = ode45(@(V,Y) ODEfun(V, Y, Other_in), Vspan, y0, opts);

Selec1(counter) = y(end,2) ./ (y(end,2) + y(end,4) + y(end,6));
Conv1(counter) = (Fa0  - y(end,1))./Fa0;

counter = counter + 1;

%plot(v, (Fa0  - y(:,1))./Fa0);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars -except Conv1

P = 2.4;
Fa0 = 10.4;
Steam_frac = 12;
two_and_three = 0;
Vspan = [0 10];
y0 = [Fa0 0 0 0 0 0 0];

K_con = [-13.1 -13050 4.27 -3.0*10^(-3) 3.0*10^(-7)];

counter = 1;

Ts = linspace(800, 1100, 100);

for T = Ts

Other_in = [
    exp(13.4 - (10930./T)),
    exp(23.9 - (25000./T)),
    exp(14.3 - (12000./T)),
    exp(K_con(1) + (K_con(2)./T) + (K_con(3).*log(T)) +(K_con(4) .*T) + (K_con(5) .* (T.^2))),
    Steam_frac*Fa0 %Steam,
    P,
    two_and_three % enable or disable reactions 2 and 3
];

opts = odeset('RelTol',1e-10,'AbsTol',1e-8);
[v, y] = ode45(@(V,Y) ODEfun(V, Y, Other_in), Vspan, y0, opts);

Selec2(counter) = y(end,2) ./ (y(end,2) + y(end,4) + y(end,6));
Conv2(counter) = (Fa0  - y(end,1))./Fa0;

counter = counter + 1;

%plot(v, (Fa0  - y(:,1))./Fa0);

end



set(gcf,'Position',[200 200 650 400])


plot(Ts, Conv1, '-k')
ylabel('Conversion in 10 m^{3} reactor')
hold on
plot(Ts, Conv2, '--r')

xlabel('Temperature / K')

title({'Conversion achieved in a 10 m^{3} isothermal reactor', 'for all reactions and for solely reaction 1'})
legend('All reactions', 'Only reaction 1', 'location', 'southeast')


set(gca,'fontname','times')

%print('Part3aexten.png', '-dpng', '-r1000')

