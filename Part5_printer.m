

output2 = output(40:3961,:);
output2 = [output2(1:1800:2410,:) ; output2(2410:50:end,:) ];
%output2 = output2(1:1:end, :);

f = fit(output2(:,1), output2(:,4), 'smoothingspline');

plot(f, output2(:,1), output2(:,4), '+w')

ylabel('Optimal conversion')
xlabel('Volume')
title({'The optimal conversion achieved for selectivity > 0.95', 'by varying T and S-E ratio'})
xlabel('Volume / m^3')
set(gca,'fontname','times')
set(gcf,'Position',[200 200 650 400])
xlim([0 65])
set(gcf,'Position',[200 200 650 400])
%print('Part5_Optimcon_vol.png', '-dpng', '-r1000')
