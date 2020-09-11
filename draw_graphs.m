y = 1:100;
x = randsample(5, 100, true);
k = x;
w = randsample(5, 100, true);
q = randsample(5, 100, true);
e = randsample(5, 100, true);

figure;
plot(y,x,'ko', 'MarkerSize', 4, 'MarkerFaceColor', 'r');
hold on;
grid on;
grid minor;
plot(y,w,'ko', 'MarkerSize', 8);
hold on;

ylabel('Selected precoder index');
xlabel('Number of trial');

xlim([0 100])
ylim([0 6])

xbounds = ylim;
set(gca,'YTick',xbounds(1):xbounds(2));

hold off;

% legend('Optimal','Proposed');

legend('Optimal','LR-based');


figure;

plot(y,x,'ko', 'MarkerSize', 4, 'MarkerFaceColor', 'r');
hold on;
grid on;
grid minor;
plot(y,k,'ko', 'MarkerSize', 8);
hold on;

ylabel('Selected precoder index');
xlabel('Number of trial');

xlim([0 100])
ylim([0 6])

xbounds = ylim;
set(gca,'YTick',xbounds(1):xbounds(2));

hold off;

legend('Optimal','LR-based');