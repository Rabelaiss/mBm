# mBm
Generates Riemann-Liouville multifractional Brownian motion paths with a given Hurst function.

## Usage
* `mbm = mBm(n,H,interval)` produces a mBm path of length `n` with Hurst function `H` evaluated at the `interval`. If `interval = []` then it is set to `[0 1]`.
* `[mbm, ts] = mBm(n,H,interval)` also produces the vector of the time steps.
* `[mbm, ts, hs] = mBm(n,H,interval)` also produces the vector of the Hurst steps, i.e. the Hurst function evaluated at the `interval`.
* `[...] = mBm(n,H,interval,fig)` plots the path and `H` if `fig = true`.

`n` = integer bigger than 1<br>
`H` = function or real number between 0 and 1<br>
`interval` = vector with two increasing components<br>
`fig` = boolean

## Examples
```
mBm(500, 0.8, [], true);
mBm(500, @(t) 0.6*t + 0.3, [], true);
mBm(500, @(t) 0.7 - 0.4 * exp(-64*(t-0.75).^2), [], true);
mBm(500, @(t) atan(t) / 3 + 0.5, [-pi pi], true);
mBm(500, @(t) sin(t) / 3 + 1/2, [0 4*pi], true);
```

## Reference
S. V. Muniandy and S. C. Lim (2001)<br>
Modeling of locally self-similar processes using multifractional Brownian motion of Riemann-Liouville type.<br>
Physical Review E 63(4 Pt 2):046104<br>
DOI: 10.1103/PhysRevE.63.046104
