function plotDOARose(doaGCC, doaMUSIC, doaKF, savePath)

edges = deg2rad(0:10:360);

figure('Visible','off');
pax = polaraxes;
hold on;

% rose plots for each method
[countG,~] = histcounts(doaGCC, edges);
[countM,~] = histcounts(doaMUSIC, edges);
[countK,~] = histcounts(doaKF, edges);

theta = deg2rad(0:10:350);

polarplot(theta, countG, 'r', 'LineWidth',1.5);
polarplot(theta, countM, 'b', 'LineWidth',1.5);
polarplot(theta, countK, 'g', 'LineWidth',2);

legend('GCC-PHAT','MUSIC','Kalman');
title('Polar Histogram of DOA Estimates (Rose Plot)');

set(gca,'ThetaZeroLocation','top','ThetaDir','clockwise');
saveas(gcf, fullfile(savePath,'doa_rose_all.png'));
close(gcf);

end
