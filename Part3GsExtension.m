clear all
clf
close all

T_conmax = [];
Steam_conmax = [];

VMAXs = round(linspace(1,8,20).^2);
for VMAX = VMAXs
    
clearvars -except VMAX VMAXs Steam_conmax T_conmax
clf
close all

P = 2.4;
Fa0 = 10.4;
Steam_frac = 12;
two_and_three = 1;
Vspan = [0 VMAX];
y0 = [Fa0 0 0 0 0 0 0];

K_con = [-13.1 -13050 4.27 -3.0*10^(-3) 3.0*10^(-7)];

inner_counter = 1;
outer_counter = 1;

T = 1000;

Steam_fracs = linspace(0, 32, 50);
T_s = linspace(800, 1100, 50);

for Steam_frac = Steam_fracs
    for T = T_s

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

        Selec1(outer_counter,inner_counter) = y(end,2) ./ (y(end,2) + y(end,4) + y(end,6));
        
        Conv1(outer_counter,inner_counter) = (Fa0  - y(end,1))./Fa0;
        
%         if Selec1(outer_counter,inner_counter) > 0.95
%             Conv1(outer_counter,inner_counter) = (Fa0  - y(end,1))./Fa0;
%         else
%             Conv1(outer_counter,inner_counter) = 0;
%         end
        
        inner_counter = inner_counter + 1;

        %plot(v, (Fa0  - y(:,1))./Fa0);
    end
    inner_counter = 1;
    outer_counter = outer_counter + 1;
    if ~rem(outer_counter,10)
        outer_counter 
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Conv1 = imresize(Conv1, [1000 1000], 'bilinear');
Selec1 = imresize(Selec1, [1000 1000], 'bilinear');

for x = [1:1:length(Conv1)]
    for y = [1:1:length(Conv1)]
        if Selec1(x,y) < 0.95
            Conv1(x,y) = 0;
        end
    end
end

T_s = imresize(T_s', [1000 1]);
Steam_fracs = imresize(Steam_fracs', [1000 1]);

% [row, col] = find(ismember(Conv1, max(Conv1(:))));
% T_conmax = [T_conmax; {T_s(col)}];
% Steam_conmax = [Steam_conmax; {Steam_fracs(row)}];

s = pcolor(T_s, Steam_fracs, Conv1);

%s = surf(T_s, Steam_fracs, Conv1, 'FaceAlpha',1)
xlabel('Temperature / K')
ylabel('Steam-ethylbenzene ratio')
colorbar
ax = gca;
ax.YDir = 'reverse';
s.EdgeColor = 'none';

caxis([0 1])

vec = [      0;   5;    35; 50;      80;             100];
hex = ['#ffffff';'#ffffff';'#FFEE00';'#FF8300';'#D30000'; '#000000'];
raw = sscanf(hex','#%2x%2x%2x',[3,size(hex,1)]).' / 255;
N = 128;
%N = size(get(gcf,'colormap'),1) % size of the current colormap
map = interp1(vec,raw,linspace(100,0,N),'pchip');

set(gcf,'Position',[200 200 650 400])
set(gca,'fontname','times')

colormap(map)
title({'Conversion variation with temperature and steam-ethylbenzene ratio', sprintf('Reactor volume %.f m^{3}',Vspan(2))})
%print(sprintf('Part3Extension(ver2)V%.f.png', Vspan(2)), '-dpng', '-r1000')

% left_color = [0 0 0];
% right_color = [0.7 0 0];
% set(figure,'defaultAxesColorOrder',[left_color; right_color]);
% set(gcf,'Position',[200 200 650 400])
% 
% 
% yyaxis left
% plot(Steam_fracs, Conv1, '-k')
% ylabel('Conversion in 10 m^{3} reactor')
% hold on
% yyaxis right
% ylabel('Selectivity in 10 m^{3} reactor')
% plot(Steam_fracs, Selec1, '-r')
% ylim([0.8 1])
% hold off
% 
% xlabel('Steam/ethylbezene ratio')
% 
% title({'Conversion and Selectivity achieved in a 10 m^{3} isothermal reactor', 'with a temperature of 1000 K'})
% legend('Conversion', 'Selectivity', 'location', 'east')
% 
% set(gca,'fontname','times')

%print('Part3b.png', '-dpng', '-r1000')

end
