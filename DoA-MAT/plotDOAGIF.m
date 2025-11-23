function plotDOAGIF(frameT, doaKF_rad, outfile, fps)
% plotDOAGIF Create a GIF animation showing a rotating DOA arrow over time.
%   frameT       - time vector (seconds) for each frame (Nx1)
%   doaKF_rad    - tracked DOA (radians) for each frame (Nx1)
%   outfile      - path to output GIF file, e.g. 'outputs/doa_animation.gif'
%   fps          - frames per second for GIF (optional, default 10)

if nargin < 4 || isempty(fps), fps = 10; end
delay = 1 / fps;

% Create outputs folder if necessary
[outdir, ~, ~] = fileparts(outfile);
if ~isempty(outdir) && ~exist(outdir,'dir')
    mkdir(outdir);
end

% figure settings
h = figure('Visible','off','Color','w','Position',[200 200 500 500]);
ax = polaraxes(h);
thetaticks(-180:45:180);
rlim([0 1.5]);
set(ax, 'ThetaZeroLocation','top','ThetaDir','clockwise');
hold on;

% Normalize radius for arrow
r = 1.0;

% Convert to degrees for annotation
degLabels = rad2deg(doaKF_rad);

% Prepare GIF
frameCount = length(frameT);
for k = 1:frameCount
    cla(ax);
    th = doaKF_rad(k);
    % plot arrow
    hLine = polarplot(ax, [th th], [0 r], 'g-', 'LineWidth', 3);
    hold on;
    % plot dot at tip
    polarplot(ax, th, r, 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 10);

    % show time and angle as text
    txt = sprintf('t = %.2fs\nangle = %.1f°', frameT(k), rad2deg(th));
    % place annotation using axes coordinates
    text(0.02, 0.02, txt, 'Units','normalized', 'FontSize', 12, 'Color','k', 'Parent', ax);

    % draw and capture
    drawnow;
    frame = getframe(h);
    im = frame2im(frame);
    [A,map] = rgb2ind(im,256);

    if k == 1
        imwrite(A,map,outfile,'gif','LoopCount',Inf,'DelayTime',delay);
    else
        imwrite(A,map,outfile,'gif','WriteMode','append','DelayTime',delay);
    end
end

close(h);
fprintf('GIF saved: %s\n', outfile);
end
