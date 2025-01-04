function plot_heart_work_space_value(e)

%e1 = coarsegraining(e,1);
e1 = coarsegraining(e,679); % max 679 factor ,Time scale 3395
figure
plot(e1)
ax = gca;
ax.FontSize = 40;
grid on;
%ylim([0.18 0.85]);
%title('Heart Rate Multiscale Fuzzy Entropy');
xlabel('Scale [sample]');
ylabel('Heart rate (Z score)');
