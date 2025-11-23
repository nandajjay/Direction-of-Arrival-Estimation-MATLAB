function [multich, fs, anglesTrue, micPos] = simulateArray(inputAudioPath, fs, duration, numMics, spacing, angleSweep)
% SIMULATEARRAY  Generate multichannel microphone array signals
%   from a single mono input using fractional delays.
%
%   INPUTS:
%       inputAudioPath : path to mono audio (.wav). If empty, uses MATLAB chirp.
%       fs             : sampling rate (Hz)
%       duration       : duration (s)
%       numMics        : number of microphones in ULA
%       spacing        : spacing in meters (e.g., 0.05m)
%       angleSweep     : [startAngle endAngle] in degrees (e.g., [-50 50])
%
%   OUTPUTS:
%       multich   : multi-channel audio (numMics × N)
%       fs        : sample rate
%       anglesTrue: instantaneous angle per sample (radians)
%       micPos    : array of mic positions (numMics × 2)

% Speed of sound
c = 343;

% --- Load or generate source ---
if ~isempty(inputAudioPath)
    [src, fsOrig] = audioread(inputAudioPath);
    
    if fsOrig ~= fs
        src = resample(src, fs, fsOrig);
    end
    
    if size(src,2) > 1
        src = mean(src,2); % convert stereo to mono
    end
    
    src = src(:)';
else
    % Generate synthetic chirp if no input is given
    t = 0:1/fs:duration;
    src = chirp(t, 300, duration, 2000);
end

N = length(src);

% --- Microphone positions (ULA) ---
micPos = zeros(numMics,2);
micPos(:,1) = (0:numMics-1) * spacing;  % x positions
micPos(:,2) = 0;                        % y = 0 for ULA

% --- Create angle sweep ---
anglesTrue = linspace(angleSweep(1), angleSweep(2), N);
anglesTrue = deg2rad(anglesTrue);

% --- Output matrix ---
multich = zeros(numMics, N);

% --- Fractional delay for each mic ---
for n = 1:N
    theta = anglesTrue(n);
    u = [cos(theta), sin(theta)];

    taus = - (micPos * u') / c;   % delays for each mic (seconds)

    for m = 1:numMics
        multich(m,n) = fractionalDelay(src, n, fs, taus(m));
    end
end

% Normalize output
maxVal = max(abs(multich), [], 'all');
multich = 0.98 * multich / maxVal;

end

% ------------------------------------------------------------
% Helper: apply fractional delay via short-window interpolation
% ------------------------------------------------------------
function out = fractionalDelay(src, n, fs, tau)
% Apply fractional delay using linear interpolation (simple, fast)
    tDelayed = n - tau*fs;
    n1 = floor(tDelayed);
    n2 = n1 + 1;

    if n1 < 1 || n2 > length(src)
        out = 0;
        return;
    end
    
    frac = tDelayed - n1;
    out = (1-frac)*src(n1) + frac*src(n2);
end
