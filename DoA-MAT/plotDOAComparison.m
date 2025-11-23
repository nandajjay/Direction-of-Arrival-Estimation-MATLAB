function plotDOAComparison(frameT, doaGCC, doaMUSIC, doaKF, trueDOA, savePath)
% plotDOAComparison
% Plots GCC, MUSIC, Kalman Filter, and True DOA vs Time.

figure('Visible','off');

plot(frameT, rad2deg(trueDOA), 'k', 'LineWidth', 2); hold on;
plot(frameT, rad2deg(doaGCC), 'r--', 'LineWidth', 1.2);
plot(frameT, rad2deg(doaMUSIC), 'b-.', 'LineWidth', 1.2);
plot(frameT, rad2deg(doaKF), 'g', 'LineWidth', 2);

xlabel('Time (s)');
ylabel('DOA (degrees)');
title('Direction of Arrival: GCC vs MUSIC vs Kalman');
legend('True DOA','GCC','MUSIC','Kalman','Location','best');
grid on;

saveas(gcf, fullfile(savePath,'doa_vs_time.png'));
close(gcf);

end
