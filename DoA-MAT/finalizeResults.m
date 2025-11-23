function metrics = finalizeResults(angleKF_rad, trueAngles, frameTimes, outDir)
% FINALIZERESULTS  Compute metrics and save plots for DOA project
%
% INPUTS:
%   angleKF_rad : nFrames x 1 tracked angles (rad)
%   trueAngles  : Nsamples x 1 true angle per sample (rad) OR true angle at frame centers
%   frameTimes  : nFrames x 1 times of frame centers (s)
%   outDir      : folder to save outputs (string)
%
% OUTPUT:
%   metrics : struct with fields meanErrDeg, medianErrDeg, pctWithin5, pctWithin10, errorsDeg

if nargin < 4 || isempty(outDir)
    outDir = 'outputs';
end
if ~exist(outDir,'dir'), mkdir(outDir); end

% Interpolate true angles to frameTimes if needed
if length(trueAngles) > length(frameTimes)
    trueAtFrames = interp1( linspace(0,1,length(trueAngles)), trueAngles, linspace(0,1,length(frameTimes)), 'linear');
else
    trueAtFrames = trueAngles;
end

% Compute circular angular error
diff = angleKF_rad - trueAtFrames(:);
diffWrapped = angle(exp(1j*(diff)));  % equivalent to wrapToPi
errorsDeg = abs(rad2deg(diffWrapped));

metrics.meanErrDeg = mean(errorsDeg(~isnan(errorsDeg)));
metrics.medianErrDeg = median(errorsDeg(~isnan(errorsDeg)));
metrics.pctWithin5 = 100 * sum(errorsDeg <= 5)/numel(errorsDeg);
metrics.pctWithin10 = 100 * sum(errorsDeg <= 10)/numel(errorsDeg);
metrics.errorsDeg = errorsDeg;

% Save CDF plot
figure('Visible','off');
% Manual CDF (no toolbox required)
err = errorsDeg(~isnan(errorsDeg));
xvals = sort(err);
fvals = (1:length(xvals)) / length(xvals);

plot(xvals, fvals, 'LineWidth', 1.6);
xlabel('Absolute Angular Error (deg)');
ylabel('CDF');
grid on;
title('CDF of DOA Absolute Error');
saveas(gcf, fullfile(outDir, 'cdf_error.png'));
close(gcf);

% Histogram
figure('Visible','off');
histogram(errorsDeg, 0:2:60);
xlabel('Absolute Error (deg)');
ylabel('Count');
title('Histogram of DOA Error (deg)');
grid on;
saveas(gcf, fullfile(outDir, 'hist_error.png'));
close(gcf);

% Save metrics to text
fid = fopen(fullfile(outDir,'metrics.txt'),'w');
fprintf(fid,'Mean ABS error (deg): %.3f\n', metrics.meanErrDeg);
fprintf(fid,'Median ABS error (deg): %.3f\n', metrics.medianErrDeg);
fprintf(fid,'Percent within 5 deg: %.2f%%\n', metrics.pctWithin5);
fprintf(fid,'Percent within 10 deg: %.2f%%\n', metrics.pctWithin10);
fclose(fid);

end
