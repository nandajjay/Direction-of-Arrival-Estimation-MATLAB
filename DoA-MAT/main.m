%% MAIN.M — Full DOA Pipeline Runner
% Runs the complete DOA estimation + tracking + beamforming pipeline.

clear; clc; close all;

disp('==============================');
disp('  DOA PROJECT — FULL PIPELINE  ');
disp('==============================');

%% ----------------------------
% USER CONFIGURATION
%% ----------------------------

inputAudio = '';          % leave '' to use synthetic chirp
fs = 16000;               % sample rate
duration = 4;             % seconds (used only if no inputAudio)
numMics = 4;              % ULA with 4 microphones
spacing = 0.05;           % 5 cm spacing
angleSweep = [-50 50];    % true angle from -50° to +50°

% Frame parameters
win = 1024;
hop = 512;

% MUSIC parameters
freqRange = [300 3000];
nSource = 1;

%% ----------------------------
% STEP 1 — SIMULATION
%% ----------------------------

disp('Step 1: Simulating microphone array...');
[multich, fs, anglesTrue, micPos] = simulateArray( ...
        inputAudio, fs, duration, numMics, spacing, angleSweep);

if ~exist('demos','dir'), mkdir demos; end
audiowrite('demos/mic1_ref.wav', multich(1,:), fs);

disp('   ✔ Simulation complete.');

%% ----------------------------
% STEP 2 — GCC-PHAT
%% ----------------------------

disp('Step 2: Running GCC-PHAT...');
[doaGCC_rad, frameCenters_s, gccParams] = runGCC(multich, fs, micPos, ...
                            'win', win, 'hop', hop, 'interp', 16);
disp('   ✔ GCC-PHAT complete.');

%% ----------------------------
% STEP 3 — MUSIC
%% ----------------------------

disp('Step 3: Running MUSIC estimator...');
[doaMUSIC_rad, frameCenters_s_music, musicParams] = runMUSIC(multich, fs, micPos, ...
                            'win', win, 'hop', hop, 'freqRange', freqRange, 'nSource', nSource);
disp('   ✔ MUSIC complete.');

%% ----------------------------
% Prepare unified measurements
%% ----------------------------

anglesTrue_music = interp1( ...
    (0:length(multich)-1)/fs, anglesTrue, frameCenters_s_music, 'linear');

measAngles = doaMUSIC_rad;
nanIdx = isnan(measAngles);
measAngles(nanIdx) = doaGCC_rad(nanIdx);

%% ----------------------------
% STEP 4 — KALMAN FILTER
%% ----------------------------

disp('Step 4: Running Kalman filter...');
angleKF_rad = angularKalman(measAngles, frameCenters_s_music, ...
                            'sigmaProcess',5, 'sigmaMeas',8);
disp('   ✔ Kalman tracking complete.');

%% ----------------------------
% GIF ANIMATION OF DOA
%% ----------------------------

if ~exist('outputs','dir'), mkdir outputs; end
gifPath = fullfile('outputs','doa_animation.gif');
plotDOAGIF(frameCenters_s_music, angleKF_rad, gifPath, 12);  % 12 fps

%% ----------------------------
% STEP 5 — DELAY-AND-SUM BEAMFORMER
%% ----------------------------

disp('Step 5: Running Delay-and-Sum beamformer...');
beamformedDS = beamformDS(multich, fs, micPos, angleKF_rad, win, hop);
audiowrite('demos/beamformed_ds.wav', beamformedDS, fs);
disp('   ✔ Delay-and-Sum beamforming complete.');

%% ----------------------------
% STEP 6 — MVDR BEAMFORMER
%% ----------------------------

disp('Step 6: Running MVDR beamformer...');
beamformedMVDR = beamformMVDR(multich, fs, micPos, angleKF_rad, win, hop, freqRange);
audiowrite('demos/beamformed_mvdr.wav', beamformedMVDR, fs);
disp('   ✔ MVDR beamforming complete.');

%% ----------------------------
% STEP 7 — FINAL METRICS + PLOTS
%% ----------------------------

disp('Step 7: Computing metrics and saving results...');

metrics = finalizeResults(angleKF_rad, anglesTrue_music, frameCenters_s_music, 'outputs');

% --- Resample GCC DOA to MUSIC timeline ---
doaGCC_resampled = interp1(frameCenters_s, doaGCC_rad, frameCenters_s_music, 'nearest', 'extrap');

% === DOA Comparison Plot ===
plotDOAComparison(frameCenters_s_music, doaGCC_resampled, doaMUSIC_rad, angleKF_rad, anglesTrue_music, 'outputs');

% === Polar DOA Plot ===
plotDOAPolar(doaGCC_resampled, doaMUSIC_rad, angleKF_rad, 'outputs');

% === Rose Plot ===
plotDOARose(doaGCC_resampled, doaMUSIC_rad, angleKF_rad, 'outputs');

% === Portfolio Layout ===
mic1 = multich(1,:);
plotPortfolioLayout(frameCenters_s_music, doaGCC_resampled, doaMUSIC_rad, ...
    angleKF_rad, anglesTrue_music, mic1, beamformedDS, beamformedMVDR, fs, 'outputs');

disp(' ');
disp('========== FINAL METRICS ==========');
disp(['Mean absolute error       : ', num2str(metrics.meanErrDeg, '%.2f'), '°']);
disp(['Median absolute error     : ', num2str(metrics.medianErrDeg, '%.2f'), '°']);
disp(['% frames within 5°        : ', num2str(metrics.pctWithin5, '%.1f'), '%']);
disp(['% frames within 10°       : ', num2str(metrics.pctWithin10, '%.1f'), '%']);
disp('===================================');

disp('All outputs saved under:');
disp('   demos/   (audio)');
disp('   outputs/ (plots, GIF, metrics)');
disp(' ');

disp('Pipeline complete ✔');
