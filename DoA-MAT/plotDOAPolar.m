function plotDOAPolar(doaGCC, doaMUSIC, doaKF, savePath)

figure('Visible','off');
pax = polaraxes;

hold on;
polarplot(pax, doaGCC, ones(size(doaGCC))*1, 'r.', 'MarkerSize', 12);
polarplot(pax, doaMUSIC, ones(size(doaMUSIC))*1, 'b.', 'MarkerSize', 12);
polarplot(pax, doaKF, ones(size(doaKF))*1.1, 'g.', 'MarkerSize', 15);

legend('GCC-PHAT','MUSIC','Kalman Filter');

title('Polar DOA Plot (All Methods)');
set(gca,'ThetaZeroLocation','top','ThetaDir','clockwise');

saveas(gcf, fullfile(savePath,'doa_polar_all.png'));
close(gcf);

end
