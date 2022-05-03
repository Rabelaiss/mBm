function [mbm, ts, hs] = mBm(n, H, interval, fig)
% mBm    Riemann-Liouville Multrifractional Brownian motion
%    mbm = mBm(n,H,interval) produces a mBm path of length n with Hurst
%    function H evaluated at the interval. If interval = [] then it is set
%    to [0 1].
%
%    [mbm, ts] = mBm(n,H,interval) also produces the vector of the time
%    steps.
%
%    [mbm, ts, hs] = mBm(n,H,interval) also produces the vector of the
%    Hurst steps, i.e. the Hurst function evaluated at the interval.
%
%    [...] = mBm(n,H,interval,fig) plots the path and H if fig = true.
%
%           n = integer bigger than 1
%           H = function or real number between 0 and 1
%    interval = vector with two increasing components
%         fig = boolean
%
% Examples
%    mBm(500, 0.8, [], true);
%    mBm(500, @(t) 0.6*t + 0.3, [], true);
%    mBm(500, @(t) 0.7 - 0.4 * exp(-64*(t-0.75).^2), [], true);
%    mBm(500, @(t) atan(t) / 3 + 0.5, [-pi pi], true);
%    mBm(500, @(t) sin(t) / 3 + 1/2, [0 4*pi], true);
%
% Reference
%    S. V. Muniandy and S. C. Lim (2001)
%    Modeling of locally self-similar processes using multifractional Brownian
%    motion of Riemann-Liouville type.
%    Physical Review E 63(4 Pt 2):046104
%    DOI: 10.1103/PhysRevE.63.046104
%

if nargin < 3
    error("At least 3 inputs required: mBm(n,H,interval) with 'n' = length of the path, 'H' = Hurst function, 'interval' = vector with two components")
elseif nargin == 3
    fig = 0;
end

if ~isscalar(n)
    error("'n' must be a number")
elseif n < 2
    error("'n' must be bigger than 1")
elseif n - floor(n) > 0
    n = round(n);
    warning("'n' must be an integer, it was rounded to the nearest integer (%d)", n)
end

if ~isa(H,'function_handle')
    if isscalar(H)
        if (H<=0) || (H>=1)
            error("'H' must be 0 < H < 1")
        else
            H = @(t) 0*t + H; % constant vector
        end
    else
        error("'H' must be a function or a number")
    end
end

if ~isnumeric(interval)
    error("'interval' must be a numeric vector")
elseif isempty(interval)
    interval = [0 1];
elseif (max(size(interval)) ~= 2) || (min(size(interval)) ~= 1)
    error("'interval' must be a vector with two components")
elseif diff(interval) <= 0
    error("'interval' must be a vector with two increasing components")
end

if ~isa(fig,'logical')
    if isscalar(fig)
        if (fig ~= 0) && (fig ~= 1)
            error("'fig' must be 0 or 1 (false or true)")
        end
    else
        error("'fig' must be 0, 1, false or true")
    end
end

ts = linspace(interval(1), interval(2), n)'; % time steps for the Hurst function
hs = H(ts); % Hurst steps

if min(hs) <= 0
    warning("The Hurst function goes below 0 (min = %.2f) while it should be in the interval (0,1)", min(hs))
elseif max(hs) >= 1
    warning("The Hurst function goes above 1 (max = %.2f) while it should be in the interval (0,1)", max(hs))
end

t0 = 0; % starting time for the mBm (0 in the paper)
tn = 1; % final time for the mBm (1 in the paper)
dt = (tn - t0) / (n - 1); % step size
ts = linspace(t0, tn, n)'; % time steps for the mBm

mbm = zeros(n,1);
xi  = normrnd(0,1 , n-1,1); % gaussian white noise = vector of random numbers from the normal distribution with mean 0 and standard deviation 1
w   = @(t, H) sqrt((t .^ (2 * H) - (t - dt) .^ (2 * H)) ./ (2 * H * dt)) ./ gamma(H + 1/2); % eq 19 of the paper https://sci-hub.se/10.1103/PhysRevE.63.046104

for k = 2:n % skip k=1 since by definition the Brownian motion starts at 0
    weights  = w(ts(2:k), hs(k));
    mbm(k) = sqrt(dt) * sum( xi(1:k-1) .* flip(weights) ); % eq 17 of the paper
end

if fig
    figure
    plot(ts,mbm,ts,hs)
    ylim([min(mbm)-0.1 max(1,max(mbm)+0.1)])
end
