function [doaGCC_rad, frameCenters_s, params] = runGCC(multich, fs, micPos, varargin)
% RUNGCC  Estimate per-frame DoA using GCC-PHAT across mic pairs
%
% USAGE:
%   [doaGCC_rad, frameCenters_s, params] = runGCC(multich, fs, micPos, 'win',1024,'hop',512,'interp',8);
%
% INPUTS:
%   multich  : M x N multichannel audio (rows = mics)
%   fs       : sampling rate (Hz)
%   micPos   : M x 2 matrix of mic (x,y) coordinates (meters)
%
% OPTIONAL NAME-VALUE:
%   'win'    : window length in samples (default 1024)
%   'hop'    : hop length in samples (default 512)
%   'interp' : integer interpolation factor for FFT zero-padding (default 8)
%   'maxTau' : maximum expected TDOA in seconds (default computed from array aperture)
%
% OUTPUTS:
%   doaGCC_rad    : nFrames x 1 array of DOA estimates (radians, NaN when no estimate)
%   frameCenters_s: nFrames x 1 time of each frame center (seconds)
%   params        : struct of used parameters (win, hop, interp, nFrames, etc.)
%
% NOTES:
%   - Uses PHAT weighting and parabolic peak interpolation to refine lag estimates.
%   - Converts pairwise TDOAs to angles using far-field plane-wave assumption.

% Parse args
p = inputParser;
addRequired(p,'multich');
addRequired(p,'fs');
addRequired(p,'micPos');
addParameter(p,'win',1024);
addParameter(p,'hop',512);
addParameter(p,'interp',8);
addParameter(p,'maxTau',[]);
parse(p, multich, fs, micPos, varargin{:});
win = p.Results.win;
hop = p.Results.hop;
interp = p.Results.interp;
maxTau = p.Results.maxTau;

[M, N] = size(multich);
if isempty(maxTau)
    % conservative max TDOA based on array diameter / c
    M = size(micPos,1);
    aperture = 0;
    for ii = 1:M
        for jj = ii+1:M
            d = norm(micPos(ii,:) - micPos(jj,:));
            if d > aperture
                aperture = d;
            end
        end
    end
    maxTau = aperture / 343; % seconds
end

% Prepare frames per channel
window = hann(win,'periodic')';
nFrames = ceil((N - win) / hop) + 1;
frames = zeros(M, nFrames, win);
for m = 1:M
    % extract frames (zero-pad last)
    for f = 1:nFrames
        startIdx = (f-1)*hop + 1;
        endIdx = startIdx + win - 1;
        seg = zeros(1,win);
        if startIdx <= N
            seg(1:min(win, N-startIdx+1)) = multich(m, startIdx:min(endIdx, N));
        end
        frames(m,f,:) = seg .* window;
    end
end

% list mic pairs
pairs = [];
for i = 1:M
    for j = i+1:M
        pairs = [pairs; i, j]; %#ok<AGROW>
    end
end

% pre-alloc
doaEst = nan(nFrames,1);
frameCenters_s = ((0:(nFrames-1))*hop + win/2) / fs;

% helper: PHAT GCC + parabolic interpolation
    function [tau] = gcc_phat_subsample(x, y, fs, interpFactor, maxTau_s)
        % x,y: 1 x L
        L = length(x);
        Nfft = 2^nextpow2(L);
        Nfft = Nfft * interpFactor; % zero-pad factor
        X = fft(x, Nfft);
        Y = fft(y, Nfft);
        R = X .* conj(Y);
        denom = abs(R);
        denom(denom < 1e-12) = 1e-12;
        R = R ./ denom;
        cc = real(ifft(R, Nfft));
        % center
        cc = [cc(end - Nfft/2 + 1:end), cc(1:Nfft/2)];
        % search within max shift
        maxShift = min(round(maxTau_s * fs * interpFactor), Nfft/2 - 1);
        center = Nfft/2 + 1;
        [~, idx] = max(abs(cc(center - maxShift:center + maxShift)));
        peakIdx = idx + (center - maxShift) - 1;
        % parabolic interpolation around peak
        if peakIdx > 1 && peakIdx < length(cc)
            y0 = cc(peakIdx-1);
            y1 = cc(peakIdx);
            y2 = cc(peakIdx+1);
            denomPar = (y0 - 2*y1 + y2);
            if abs(denomPar) > 1e-12
                delta = 0.5*(y0 - y2)/denomPar; % correction in samples
            else
                delta = 0;
            end
        else
            delta = 0;
        end
        shiftSamples = (peakIdx - center) + delta;
        tau = shiftSamples / (fs * interpFactor);
    end

% helper: convert pairwise tau -> angle (far-field)
    function theta = tdoa_to_angle(tau, pos_i, pos_j)
        % Solve (p_j - p_i)·u = c * tau  with u = [cos theta; sin theta]
        d = pos_j - pos_i;
        dist = norm(d);
        if dist < 1e-8
            theta = NaN;
            return;
        end
        phi = atan2(d(2), d(1));
        val = (343 * tau) / dist;
        if abs(val) > 1
            theta = NaN;
            return;
        end
        % two possible solutions: phi +/- acos(val). Choose one in [-pi, pi].
        th1 = phi + acos(val);
        th2 = phi - acos(val);
        % map to [-pi,pi]
        th1 = mod(th1 + pi, 2*pi) - pi;
        th2 = mod(th2 + pi, 2*pi) - pi;
        % pick one closer to broadside (heuristic) - choose th2
        theta = th2;
    end

% Process frames
for f = 1:nFrames
    thetas = [];
    for pIdx = 1:size(pairs,1)
        i = pairs(pIdx,1); j = pairs(pIdx,2);
        x = squeeze(frames(i,f,:))';
        y = squeeze(frames(j,f,:))';
        try
            tau = gcc_phat_subsample(x, y, fs, interp, maxTau);
            th = tdoa_to_angle(tau, micPos(i,:), micPos(j,:));
            if ~isnan(th)
                thetas(end+1) = th; %#ok<AGROW>
            end
        catch
            % on numerical issues skip
        end
    end
    if ~isempty(thetas)
        % circular mean
        s = mean(sin(thetas));
        cval = mean(cos(thetas));
        doaEst(f) = atan2(s, cval);
    else
        doaEst(f) = NaN;
    end
end

% return
doaGCC_rad = doaEst;
params.win = win;
params.hop = hop;
params.interp = interp;
params.nFrames = nFrames;
params.pairs = pairs;
params.frameCenters = frameCenters_s;

end


