clear all
clf
close all

P = 2.4;
Fa0 = 10.4;
Steam_frac = 12;
two_and_three = 1;
Vspan = [0 10];


Thermo_data(:,1) = [30 147 0 83 52 50 -75 -242];
Thermo_data(:,1) = Thermo_data(:,1) .*1000;
Thermo_data(:,2) = [6.3 12.8 28.8 -1.8 11.1 1.6 13.9 30.7];
Thermo_data(:,3) = [0.469 0.420 0.000 0.330 0.119 0.396 0.075 0.089];
Thermo_data(:,4) = [-15.7 -14.4 0.138 -11.4 -3.61 -13.3 -1.57 0.146];
Thermo_data(:,4) = Thermo_data(:,4) .*(10^(-5));


counter = 1;

T = 1000;

Steam_fracs2 = linspace(12, 32, 40);

for Steam_frac = Steam_fracs2

y0 = [Fa0 0 0 0 0 0 0 T]; %Initial conditions
    
Other_in = [
    Steam_frac*Fa0 %Steam,
    P,
];

opts = odeset('RelTol',1e-10,'AbsTol',1e-8);
[v, y] = ode45(@(V,Y) ODEfun_adiabatic(V, Y, Other_in, Thermo_data), Vspan, y0, opts);

Selec2(counter) = y(end,2) ./ (y(end,2) + y(end,4) + y(end,6));
Conv2(counter) = (Fa0  - y(end,1))./Fa0;

%plot(v, (Fa0  - y(:,1))./Fa0);
%plot(v, y(:, 8))
T_variation_i = imresize(y(:, 8), [100 1], 'bilinear');
v_variation_i = imresize(v, [100 1], 'bilinear');

T_variation(:,counter) = T_variation_i;
v_variation(:,counter) = v_variation_i;

counter = counter + 1;


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('Part3bFilesOut.mat')

% left_color = [0 0 0];
% right_color = [0.7 0 0];
% set(figure,'defaultAxesColorOrder',[left_color; right_color]);
% set(gcf,'Position',[200 200 650 400])
% 
% yyaxis left
% plot(Steam_fracs2, Conv2, '-k')
% ylabel('Conversion in 10 m^{3} reactor')
% hold on
% plot(Steam_fracs, Conv1, '--k')
% yyaxis right
% ylabel('Selectivity in 10 m^{3} reactor')
% plot(Steam_fracs2, Selec2, '-r')
% plot(Steam_fracs, Selec1, '--r')
% ylim([0.8 1])
% hold off
% 
% xlabel('Steam/ethylbenzene ratio')
% set(gca,'fontname','times')
% 
% title({'Conversion and Selectivity achieved in a 10 m^{3} adiabatic and isothermal reactor', 'with inlet temperature of 1000 K'})
% legend('Conversion (adiabatic)', 'Conversion (isothermal)', 'Selectivity (adiabatic)', 'Selectivity (isothermal)', 'location', 'southwest')
%%print('Part4b.png', '-dpng', '-r1000')

plot(v_variation(:,1), T_variation(:,1) ./ T_variation(1,1), '-', 'color', '[0 0 1]')
hold on
plot(v_variation(:,5), T_variation(:,5) ./ T_variation(1,5), '-', 'color', '[0.25 0 0.75]')
plot(v_variation(:,10), T_variation(:,10) ./ T_variation(1,10), '-','color', '[0.5 0 0.5]')
plot(v_variation(:,20), T_variation(:,20) ./ T_variation(1,20), '-','color', '[0.75 0 0.25]')
plot(v_variation(:,40), T_variation(:,40) ./ T_variation(1,40), '-','color', '[1 0 0]')
hold off

set(gcf,'Position',[200 200 650 400])
set(gca,'fontname','times')

ylabel('Temperature as a fraction of inlet temperature')
xlabel('Volume of reactor reactants already passed through / m^3')
title({'Temperature variation along the reactor with inlet T = 1000 K', 'for different steam-ethylbenzene ratios'})
legend('Fraction 12', 'Fraction 14', 'Fraction 17', 'Fraction 22', 'Fraction 32', 'location', 'southwest')

