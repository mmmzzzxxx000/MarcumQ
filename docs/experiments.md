# Experiment-to-manuscript map

| Manuscript content | Code | Recorded data |
|---|---|---|
| Original Algorithm 939 partition, orders 1 and 20 | plots/section2/regions.py | Generated from the original partition rules |
| Four-method timing maps and summary | experiments/section51/run.py; plots/section51/run.py | data/section51 |
| Updated hybrid timing maps | experiments/section52/grid/run.py; plots/section52/grid.py | data/section52/grid |
| Bound B heatmaps | plots/section4/bound.py | Generated from analytic bound |
| Original coefficient decay figure | plots/section3/coefficient_decay.py | data/section3/coefficient_decay.csv (200 original MATLAB figure values) |
| Coefficient error theorem | experiments/section3/recurrence.py | Recompute with exact rational arithmetic |
| Cancellation examples | experiments/section4/cancellation.py | data/section4/cancellation.json |
| Full-domain random hybrid timing | experiments/section52/random/run.py; plots/section52/speedup.py | data/section52/random |
| Large random test | experiments/section53/validate.f90; run.py | data/section53 |
| CDF comparison | experiments/section61/run.py | data/section61 |
| Outage comparison | experiments/section62/run.py | data/section62 |

The implemented early branch switch in B is M<0.998, as used in the accepted figures and tests. The bound threshold is 1000. This value is a sufficient routing criterion, not the empirical success definition. Success requires at least one tested N to attain the actual relative-error target; monotonic error decrease is not required.

The random test samples a,b,nu independently on logarithmic scales over [0.01,100], [0.001,7], [0.1,1000], with the recorded Fortran seed 20261109. The reference is released Algorithm 939, unlike the Bessel-integral performance references. Its documented order domain starts at 1; results for nu<1 and reported underflow outcomes remain separately classified. Do not describe every reference value as arbitrary-precision or every nu<1 outcome as certified by the 939 domain. Fortran intrinsic random streams can differ between compiler versions.

Full-plane gray cells retain their original meaning: failure to attain accuracy within the tested work range, or timing above the display scale. The color limits are 100 us for Python and 2 us for C++/Fortran. No measurements are interpolated or smoothed.

CDF: six fixed (M,a) configurations, 1000 thresholds including both endpoints of [0,5]. Both GammaCHI and the hybrid use scalar and batch modes and externally selected work counts at relative error 1e-10. Outage: 31 integer K values plus 29.5 and 30.5, each with three thresholds (99 cases), with direct lower-tail output. Recorded times are hardware-dependent.

Only paths, imports, source organization and build orchestration have been changed during collection. Numerical kernels and experiment calibration/timing bodies are carried over from their selected current sources. Build-time changes exposing fixed work and caches are retained and recorded in generated build manifests. Historical absolute paths inside recorded result metadata identify the original experiment; they are not runtime dependencies.

`python -B run.py paper-figures` recreates all ten manuscript image files and exports the numerical tables. `python -B run.py paper-audit` checks all recorded table/text statistics. The complete image mapping is in `manuscript_coverage.json`. The public hybrid interfaces are in `src/python/marcum/hybrid.py`, `src/cpp/hybrid.cpp`, and `src/fortran/hybrid.f90`. See `interfaces.md` for language sharing and fixed-work semantics. The random timing sample and the random accuracy sample have separate drivers and records. See the README for their domains and the treatment of tiny reference values.

Chapter layout: random timing belongs to `section52/random`, random accuracy belongs to `section53`, and coefficient/cancellation checks belong to `section3`/`section4`. Matching data and plot directories follow the same chapter names. The public language wrappers were extracted around the retained component kernels; no manuscript timing record was recomputed or replaced by the reorganization.
