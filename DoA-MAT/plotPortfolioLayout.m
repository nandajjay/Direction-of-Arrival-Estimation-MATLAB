function plotPortfolioLayout(frameT, doaGCC, doaMUSIC, doaKF, trueDOA, mic1, beamDS, beamMVDR, fs, outDir)
% plotPortfolioLayout Create a polished multi-panel figure for portfolio.

if ~exist(outDir,'dir'), mkdir(outDir); end

fig = figure('Visible','off','Units','normalized',...
    'Position',[0.05 0.05 0.9 0.8],'Color','w');

%% ----------------------------
% Panel 1 — DOA vs Time
%% ----------------------------
ax1 = subplot(2,3,1);
plot(frameT, rad2deg(trueDOA), 'k-', 'LineWidth', 1.8); hold on;
plot(frameT, rad2deg(doaGCC), 'r--', 'LineWidth', 1);
plot(frameT, rad2deg(doaMUSIC), 'b-.', 'LineWidth', 1);
plot(frameT, rad2deg(doaKF), 'g-', 'LineWidth', 1.6);
xlabel('Time (s)'); ylabel('Angle (°)');
title('DOA vs Time');
legend('True','GCC','MUSIC','Kalman','Location','best');
grid on; box on;

%% ----------------------------
% Panel 2 — Polar DOA (GCC/MUSIC/KF)
%% ----------------------------
temp = subplot(2,3,2);   % placeholder
pos = get(temp,'Position');
delete(temp);

ax2 = polaraxes('Position',pos);
hold(ax2,'on');
polarplot(ax2, doaGCC, ones(size(doaGCC))*1.05, 'r.', 'MarkerSize', 10);
polarplot(ax2, doaMUSIC, ones(size(doaMUSIC))*1.0, 'b.', 'MarkerSize', 12);
polarplot(ax2, doaKF, ones(size(doaKF))*0.95, 'g.', 'MarkerSize', 14);
title(ax2,'Polar DOA (GCC / MUSIC / KF)');
set(ax2,'ThetaZeroLocation','top','ThetaDir','clockwise');

%% ----------------------------
% Panel 3 — Rose Plot
%% ----------------------------
temp = subplot(2,3,3);
pos = get(temp,'Position');
delete(temp);

ax3 = polaraxes('Position',pos);
edges = deg2rad(0:15:360);
[countK,~] = histcounts(doaKF, edges);
theta = deg2rad(0:15:345);
polarplot(ax3, theta, countK, 'g', 'LineWidth', 2);
title(ax3,'Rose Plot (KF)');
set(ax3,'ThetaZeroLocation','top','ThetaDir','clockwise');

%% ----------------------------
% Panel 4 — Reference Mic waveform
%% ----------------------------
ax4 = subplot(2,3,4);
tSig = (0:length(mic1)-1)/fs;
plot(tSig, mic1, 'Color',[0.6 0.6 0.6]);
xlabel('Time (s)'); ylabel('Amplitude');
title('Mic 1 (Reference)');
grid on;

%% ----------------------------
% Panel 5 — Beamformed vs Mic1
%% ----------------------------
ax5 = subplot(2,3,5);
plot(tSig, mic1, 'Color',[0.7 0.7 0.7]); hold on;
plot(tSig, beamDS, 'b', 'LineWidth', 1);
plot(tSig, beamMVDR, 'r', 'LineWidth', 1);
xlabel('Time (s)');
title('Beamformed: DS (blue) vs MVDR (red)');
legend('Mic1','DS','MVDR');
grid on;

%% ----------------------------
% Panel 6 — CDF + Histogram (Embedded Images)
%% ----------------------------
ax6 = subplot(2,3,6);

cdfPNG  = '/mnt/data/cdf_error.png';
histPNG = '/mnt/data/hist_error.png';

if exist(cdfPNG,'file') && exist(histPNG,'file')
    imCDF = imread(cdfPNG);
    imshow(imCDF, 'Parent', ax6);
    title(ax6,'CDF of Error');
else
    text(0.5,0.5,'CDF/Hist images not found','HorizontalAlignment','center');
end

%% ----------------------------
% Save
%% ----------------------------
outFile = fullfile(outDir, 'portfolio_layout.png');
saveas(fig, outFile);
close(fig);

fprintf('Saved portfolio layout: %s\n', outFile);
end