function [doaMUSIC_rad, frameCenters_s, params] = runMUSIC(multich, fs, micPos, varargin)
% RUNMUSIC  Manual MUSIC DoA estimator (no toolbox required)
%
%   Broadband MUSIC using manual covariance eigendecomposition and
%   steering vector scanning.
%
% INPUTS:
%   multich  : M x N multichannel audio
%   fs       : sampling rate
%   micPos   : M x 2 microphone positions
%
% OPTIONS:
%   'win'      : frame size (default 1024)
%   'hop'      : hop size (default 512)
%   'freqRange': [fmin fmax] (default [300 3000])
%   'nSource'  : number of sources (default 1)
%
% OUTPUTS:
%   doaMUSIC_rad : DoA estimates per frame
%   frameCenters_s : time stamps
%

p = inputParser;
addParameter(p,'win',1024);
addParameter(p,'hop',512);
addParameter(p,'freqRange',[300 3000]);
addParameter(p,'nSource',1);
parse(p,varargin{:});

win       = p.Results.win;
hop       = p.Results.hop;
freqRange = p.Results.freqRange;
nSource   = p.Results.nSource;

[M, N] = size(multich);

%% -----------------------------
% STFT
%% -----------------------------
window = hann(win,'periodic')';
overlap = win - hop;
nfft = 2^nextpow2(win);

[S,freqs,frameCenters_s] = stft(multich.', fs, ...
    'Window', window, 'OverlapLength', overlap, 'FFTLength', nfft);

% S original shape: nFreq x nFrames x M
% After permute: M x nFreq x nFrames
S = permute(S, [3 1 2]);

% correct number of frames
nFrames = size(S, 3);

%% -----------------------------
% MUSIC SETUP
%% -----------------------------
% Frequency bins
freqMask = freqs >= freqRange(1) & freqs <= freqRange(2);
fBins = find(freqMask);

% Angle scan grid
scanAngles = -90:1:90;
doaMUSIC_rad = nan(nFrames,1);

%% -----------------------------
% MAIN LOOP
%% -----------------------------
for t = 1:nFrames

    % Spectral snapshot
    Xf = S(:,:,t);   % M x nFreq

    % Broadband covariance
    Xb = Xf(:,fBins);
    if size(Xb,2) < 2
        continue;
    end

    R = (Xb * Xb.') / size(Xb,2);

    % Eigen-decomposition
    [V,D] = eig(R);
    [~,idx] = sort(diag(D),'ascend');

    % Noise subspace
    En = V(:, idx(1:end-nSource));

    % MUSIC spectrum
    Pmusic = zeros(1,length(scanAngles));
    for ai = 1:length(scanAngles)
        th = deg2rad(scanAngles(ai));
        u = [cos(th); sin(th)];
        a = exp(-1j * 2*pi*(freqRange(1)/343) * (micPos*u));
        Pmusic(ai) = 1 / (norm(En' * a)^2);
    end

    % Peak angle
    [~,idxMax] = max(Pmusic);
    doaMUSIC_rad(t) = deg2rad(scanAngles(idxMax));

end

% return parameters
params.win = win;
params.hop = hop;
params.freqRange = freqRange;

end