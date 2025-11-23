function beamformed = beamformDS(multich, fs, micPos, anglesTracked, win, hop)
% BEAMFORMDS  Delay-and-sum beamformer (framewise fractional-delay)
%
%   beamformed = beamformDS(multich, fs, micPos, anglesTracked, win, hop)
%
% INPUTS:
%   multich       : M x N multichannel audio
%   fs            : sample rate
%   micPos        : M x 2 microphone positions (meters)
%   anglesTracked : nFrames x 1 tracked angles (radians)
%   win           : window length (samples)
%   hop           : hop length (samples)
%
% OUTPUT:
%   beamformed    : 1 x N beamformed (enhanced) signal
%

c = 343;
[M, N] = size(multich);

% compute number of frames
nFrames = length(anglesTracked);

beamformed = zeros(1, N);
winVec = hann(win, 'periodic')';

% FIXED Nfft (MATLAB version)
Nfft = 2^(nextpow2(win*4));   % <-- FIXED HERE

for fIdx = 1:nFrames
    startIdx = (fIdx-1)*hop + 1;
    endIdx   = startIdx + win - 1;
    if startIdx > N, break; end

    % extract windowed mic frames
    segs = zeros(M, win);
    for m = 1:M
        idx2 = min(endIdx, N);
        seg = zeros(1, win);
        seg(1:(idx2-startIdx+1)) = multich(m, startIdx:idx2);
        segs(m,:) = seg .* winVec;
    end

    theta = anglesTracked(fIdx);

    if isnan(theta)
        alignedSum = sum(segs,1);
    else
        u = [cos(theta); sin(theta)];
        taus = -(micPos * u) / c;

        aligned = zeros(M, win);
        for m = 1:M
            X = fft(segs(m,:), Nfft);
            freqs = (0:Nfft-1)*(fs/Nfft);
            phase = exp(-1j*2*pi*freqs*(-taus(m)));
            y = real(ifft(X .* phase, Nfft));
            aligned(m,:) = y(1:win);
        end
        alignedSum = sum(aligned,1);
    end

    % overlap-add
    outEnd = min(endIdx, N);
    beamformed(startIdx:outEnd) = beamformed(startIdx:outEnd) + alignedSum(1:(outEnd-startIdx+1));
end

% normalize
if max(abs(beamformed)) > 0
    beamformed = 0.98 * beamformed / max(abs(beamformed));
end
end

% helper
function p = nextpow2(x)
    p = ceil(log2(x));
end