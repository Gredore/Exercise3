clear all

T = linspace(800, 1100, 100);
k1 = exp(13.4 - (10930./T));
k2 = exp(23.9 - (25000./T));
k3 = exp(14.3 - (12000./T));
    
plot(T, k1, '-k')
hold on
plot(T,k2, '--k')
plot(T,k3, '.r')
hold off
legend