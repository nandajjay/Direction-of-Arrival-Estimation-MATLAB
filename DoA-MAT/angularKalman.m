
function angleKF_rad = angularKalman(measAngles, frameTimes, varargin)
% ANGULARKALMAN  Simple 2-state Kalman Filter for angle tracking
%
% State: [theta; omega]
% Measurement: theta_measured
%
% Handles angle wrapping and missing values (NaN).
%
% USAGE:
%   angleKF = angularKalman(measAngles, frameTimes);
%
% OPTIONALS:
%   'sigmaProcess' : process noise (default 5°/s)
%   'sigmaMeas'    : measurement noise in degrees (default 10°)

p = inputParser;
addRequired(p,'measAngles');
addRequired(p,'frameTimes');
addParameter(p,'sigmaProcess',5);
addParameter(p,'sigmaMeas',10);
parse(p, measAngles, frameTimes, varargin{:});

sigmaProcess = deg2rad(p.Results.sigmaProcess);
sigmaMeas = deg2rad(p.Results.sigmaMeas);

N = length(measAngles);
angleKF_rad = nan(N,1);

% initial state
theta = measAngles(1);
if isnan(theta)
    theta = 0;
end
omega = 0;

x = [theta; omega];
P = eye(2);

for k = 2:N
    dt = frameTimes(k) - frameTimes(k-1);
    if dt <= 0, dt = 1e-3; end

    % State transition
    A = [1 dt; 0 1];
    Q = [sigmaProcess^2 0; 0 (sigmaProcess/2)^2];

    % Predict
    x = A*x;
    P = A*P*A' + Q;

    % Measurement
    z = measAngles(k);

    if ~isnan(z)
        % Angle wrapping
        innov = wrapToPi(z - x(1));

        % Measurement model
        H = [1 0];
        R = sigmaMeas^2;

        % Kalman gain
        S = H*P*H' + R;
        K = P*H'/S;

        % Update
        x = x + K*innov;
        P = (eye(2) - K*H) * P;
    end

    % Save
    angleKF_rad(k) = wrapToPi(x(1));
end

end
