function beamformed = beamformMVDR(multich, fs, micPos, anglesTracked, win, hop, freqRange)
% BEAMFORMMVDR  Narrowband MVDR beamformer applied framewise then overlap-add
%
% USAGE:
%   beamformed = beamformMVDR(multich, fs, micPos, anglesTracked, win, hop, freqRange)
%
% INPUTS:
%   multich       : M x N multichannel audio
%   fs            : sampling rate
%   micPos        : M x 2 microphone positions (meters)
%   anglesTracked : nFrames x 1 tracked angles (radians)
%   win, hop      : frame parameters (samples)
%   freqRange     : [fmin fmax] (Hz) used to select frequency bins for covariance averaging (optional)
%
% OUTPUT:
%   beamformed    : 1 x N beamformed signal
%
% NOTES:
%   - This implements a simple per-frame, narrowband MVDR using covariance estimated
%     from the STFT bins in freqRange at that frame.
%   - More robust implementations estimate continuous frequency covariances across time.
%
% Author: Generated for DOA project

if nargin < 7 || isempty(freqRange)
    freqRange = [300 fs/2*0.9];
end

c = 343;
[M, N] = size(multich);

% STFT parameters (used internally)
nfft = max(1024, 2^nextpow2(win));  % resolution
window = hann(win,'periodic')';
overlap = win - hop;

% Compute STFT for each mic: result dims (nfft/2+1) x nFrames x M
[S, f, t] = stft(multich.', fs, 'Window', window, 'OverlapLength', overlap, 'FFTLength', nfft);
% S: time-frequency across columns: (nFreq x nTime x M) after repositioning
nFreq = size(S,1);
nFrames = size(S,2);

% Permute to (M x nFreq x nFrames)
Sperm = permute(S, [3 1 2]); % M x nFreq x nFrames

% frequency mask
freqMask = (f >= freqRange(1)) & (f <= freqRange(2));
freqIdxs = find(freqMask);
if isempty(freqIdxs)
    freqIdxs = 1:nFreq;
end

% prepare output STFT for beamformed signal (nFreq x nFrames)
Y = zeros(nFreq, nFrames);

% Precompute steering vectors for all angles and freqs used
angles_deg_grid = rad2deg(anglesTracked(:))'; % per-frame angles in degrees
% but steering computed per frame & per freq

for frameIdx = 1:nFrames
    % choose tracked angle for this frame (rad)
    theta = anglesTracked(min(frameIdx, length(anglesTracked)));
    if isnan(theta)
        theta = 0;
    end
    u = [cos(theta); sin(theta)];
    
    % For each selected frequency bin, compute covariance and MVDR weight
    for fi = freqIdxs'
        % Snapshot across mics at freq fi, this frame
        x_f = Sperm(:,fi,frameIdx);  % M x 1 vector (complex)
        % Use neighboring freqs (freqRange) around fi to estimate covariance robustly
        % Collect snapshots for covariance: use same fi across small freq neighborhood
        % We'll take indices in freqRange for covariance estimation
        snapshots = squeeze(Sperm(:,freqIdxs,frameIdx)); % M x K
        if size(snapshots,2) < 1
            R = (x_f * x_f') + 1e-6*eye(M);
        else
            R = (snapshots * snapshots') / size(snapshots,2);
            % regularize
            R = R + 1e-6 * trace(R)/M * eye(M);
        end
        % steering vector
        freqHz = f(fi);
        k = 2*pi*freqHz / c;
        a = exp(-1j * k * (micPos * u));
        a = a(:); % ensure column
        
        % MVDR weight: w = R^{-1} a / (a' R^{-1} a)
        % use backslash for stability
        Rinva = R \ a;
        denom = (a' * Rinva);
        if abs(denom) < 1e-12
            w = (Rinva);
        else
            w = Rinva / denom;
        end
        
        % apply weights to current frequency snapshot across mics
        Y(fi, frameIdx) = (w' * x_f);
    end
end

% For frequencies not computed (outside freqIdxs), fallback to sum across mics (delay-and-sum style)
otherIdxs = setdiff(1:nFreq, freqIdxs);
if ~isempty(otherIdxs)
    for fi = otherIdxs
        % simple average across mics
        x_f = squeeze(Sperm(:,fi,:)); % M x nFrames
        % average across mics for each frame
        Y(fi,:) = mean(x_f,1);
    end
end

% Inverse STFT back to time-domain
% Build complex STFT matrix needed by istft: size (nFreq x nFrames)
% S (original) was nFreq x nFrames x channels; Y is nFreq x nFrames
y = istft(Y, fs, 'Window', window, 'OverlapLength', overlap, 'FFTLength', nfft);

% Ensure length matches original
y = real(y(:))';
len = size(multich,2);
if length(y) < len
    y = [y, zeros(1, len-length(y))];
elseif length(y) > len
    y = y(1:len);
end

% Normalize
if max(abs(y)) > 0
    beamformed = 0.98 * y / max(abs(y));
else
    beamformed = y;
end

end

function p = nextpow2(x)
    p = ceil(log2(x));
    % used by calling code as 2^p sometimes; here we return exponent for convenience
end
