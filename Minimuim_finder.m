clear all
close all

num_points = 50;
Volumes = linspace(1, 3, 10);
Volumes(:,end+1:end+20) = linspace(3, 3.2, 20);
Volumes(:,end+1:end+num_points-30) = linspace(3.2,64,num_points-30);
counter = 1;
for V = Volumes
    Aeq = [];
    beq = [];
    lb = [];
    ub = [];
    nonlcon = @(in) Min_fun2(in);
    options = optimoptions('fmincon', 'Algorithm','interior-point', 'TolFun', 10^(-9), 'Display', 'off');

    in0 = [9 1100 22];
    A = [1 0 0 
         -1 0 0 ;
         0 1 0 ;
         0 0 1 ;
         0 0 -1] ;
    B = [V ; -0.1; 1100; 32; -12];
    x = fmincon(@(in) Min_fun(in)*-1, in0, A, B, Aeq,beq,lb,ub,nonlcon,options);
    x(:,end+1) = Min_fun(x);
    output(counter, :) = x; 
    counter = counter + 1
end

output = imresize(output, [4000 4], 'bilinear');

vec = [      0;   20;    50; 70;      90;             100];
hex = ['#ffffff';'#ffffff';'#FFEE00';'#FF8300';'#D30000'; '#000000'];
raw = sscanf(hex','#%2x%2x%2x',[3,size(hex,1)]).' / 255;
N = 1000;
%N = size(get(gcf,'colormap'),1) % size of the current colormap
map = interp1(vec,raw,linspace(100,0,N),'pchip');
c = colormap(map);

%h=colormapline(output(:,2), output(:,3), [], c)
h = scatter(output(:,2),output(:,3),[],output(:,1), '.');
ax = gca;
ax.YDir = 'reverse'; 

col = colorbar;
ylabel(col, 'Volume / m^3')
caxis([1,64]);

set(gcf,'Position',[200 200 650 400])
xlabel('Temperature / K')
ylabel('Steam-ethylbenzene ratio')
title({'Optimal conditions as volume is changed', '(ensuring selectivity > 0.95)'})
set(gca,'fontname','times')

% drawnow
% 
% vec = [      0;   20;    50; 70;      90;             100];
% hex = ['#ffffff';'#ffffff';'#FFEE00';'#FF8300';'#D30000'; '#000000'];
% raw = sscanf(hex','#%2x%2x%2x',[3,size(hex,1)]).' / 255;
% N = num_points;
% %N = size(get(gcf,'colormap'),1) % size of the current colormap
% map = interp1(vec,raw,linspace(100,0,N),'pchip');
% 
% map = uint8(map'*255); % need a 4xN uint8 array
% map(4,:) = 255;
% 
% set(h.Edge,'ColorBinding','interpolated','ColorData',map)
% 
% map = interp1(vec,raw,linspace(100,0,N),'pchip');
% 
% vec = [      0;   5;    35; 50;      80;             100];
% hex = ['#ffffff';'#ffffff';'#FFEE00';'#FF8300';'#D30000'; '#000000'];
% raw = sscanf(hex','#%2x%2x%2x',[3,size(hex,1)]).' / 255;
% N = 100;
% %N = size(get(gcf,'colormap'),1) % size of the current colormap
% map = interp1(vec,raw,linspace(100,0,N),'pchip');
% map = map;
% colormap(map)
% caxis([1,64]);
% colorbar