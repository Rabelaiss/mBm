
% this code was used to create the video youtu.be/2xaTte3XxYU

paths = 300;
pause_time = 0.02;

for h = linspace(0.1,0.9,70)
    mBm(paths, @(t) 0*t+h, [], true);
    pause(pause_time)
    clf
end
for h = linspace(0.9,0.1,70)
    mBm(paths, @(t) 0*t+h, [], true);
    pause(pause_time)
    clf
end

for h = linspace(0,0.8,70)
    mBm(paths, @(t) 0.1 + h*t, [], true);
    pause(pause_time)
    clf
end
for h = linspace(0.8,0.4,45)
    mBm(paths, @(t) (0.9-h) + (2*h-0.8)*t, [], true);
    pause(pause_time)
    clf
end

for h = linspace(0.5,0.757,45)
    mBm(paths, @(t) h + 0*t, [], true);
    pause(pause_time)
    clf
end
for h = [linspace(20,5,15) linspace(5,-5,70)]
    mBm(paths, @(t) atan(2*t) / 6 + 0.5, [-pi pi]+h, true);
    pause(pause_time)
    clf
end

for h = linspace(0,1,45)
    mBm(paths, @(t) sin(h*t) / 3 + (h+1)/4, [0 4*pi], true);
    pause(pause_time)
    clf
end
for h = linspace(1,4,90)
    mBm(paths, @(t) sin(h*t) / 3 + 1/2, [0 4*pi], true);
    pause(pause_time)
    if h < 4
        clf
    end
end


function [mbm, ts, hs] = mBm(n, H, interval, fig)
% this function is a slightly modified version of the function in the file mBm.m
% in particular, the only differences are the following
% NEW    @ 126  rng('default')
% NEW    @ 134  mbm = mbm / 4;
% REMOVE @ 137  figure
% EDIT   @ 138  ylim([-0.7 1])

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

t0 = 0; % starting time for the mBm (0 in the paper, see Reference)
tn = 1; % final time for the mBm (1 in the paper)
dt = (tn - t0) / (n - 1); % step size
ts = linspace(t0, tn, n)'; % time steps for the mBm

mbm = zeros(n,1);
% rng('default') fix the random number generator used by RAND, RANDI, and RANDN
% so that xi (and so the generated mBm trajectory too) is always the same
rng('default')
xi  = normrnd(0,1 , n-1,1); % gaussian white noise = vector of random numbers from the normal distribution with mean 0 and standard deviation 1
w   = @(t, H) sqrt((t .^ (2 * H) - (t - dt) .^ (2 * H)) ./ (2 * H * dt)) ./ gamma(H + 1/2); % eq 19 of the paper

for k = 2:n % skip k=1 since by definition the Brownian motion starts at 0
    weights  = w(ts(2:k), hs(k));
    mbm(k) = sqrt(dt) * sum( xi(1:k-1) .* flip(weights) ); % eq 17 of the paper
end
mbm = mbm / 4;

if fig
    plot(ts,mbm,ts,hs)
    ylim([-0.7 1])
end

end