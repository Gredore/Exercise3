clear all
clf
close all


P = 2.4;
Fa0 = 10.4;
Steam_frac = 0;
two_and_three = 0;
Vspan = [0 1000];
Conv_to_achieve_of_max = 0.95;
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

conv = (Fa0  - y(:,1))./Fa0;



[conv, index] = unique(conv); 

Con_req = ((Fa0  - y(end,1))./Fa0).*Conv_to_achieve_of_max;

v = v(index);

V_95conv1(counter) = interp1(conv,v,Con_req,'spline');

counter = counter + 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

counter = 1;

Steam_frac = 12;

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

conv = (Fa0  - y(:,1))./Fa0;


[conv, index] = unique(conv); 

Con_req = ((Fa0  - y(end,1))./Fa0).*Conv_to_achieve_of_max;

v = v(index);

%plot(v, conv);
%xlim([0 50])

V_95conv2(counter) = interp1(conv,v,Con_req,'spline');

counter = counter + 1;
end


left_color = [0 0 0];
right_color = [0.7 0 0];
set(figure,'defaultAxesColorOrder',[left_color; right_color]);
set(gcf,'Position',[200 200 650 400])

yyaxis right
plot(Ts, V_95conv1, '--r')
ylabel({sprintf('Volume to achieve %.f %% of max conversion',Conv_to_achieve_of_max*100), 'No steam'})
hold on
yyaxis left
plot(Ts, V_95conv2, '-k')
ylabel({sprintf('Volume to achieve %.f %% of max conversion',Conv_to_achieve_of_max*100), 'Steam-ethylbenzene inlet ratio of 12'})
hold off

xlabel('Temperature / K')
title('Volume to achieve 95% of maximum Equilibrium Conversion for Reaction 1')
legend('Steam/ethylbenzene initial ratio 12', 'Absence of steam', 'location', 'east')

set(gca,'fontname','times')

print('Part2.png', '-dpng', '-r1000')

