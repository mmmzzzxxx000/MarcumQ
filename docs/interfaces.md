# Hybrid entry points

Run `python -B run.py build` from the checkout to build all native libraries. The code currently targets Windows x86-64 and requires GNU C++ and Fortran. Generated libraries and patched source copies go under `build/`; original third-party files are not edited.

## Python

```python
from marcum.hybrid import select_method, hybrid_q
method = select_method(8.0, 2.0, 1.0)
value = hybrid_q(8.0, 2.0, 1.0, 80, language="python")
```

`language` can be `python`, `cpp` or `fortran`. It defaults to `fortran`, preserving the earlier interface. `proposed=False` uses the original Algorithm 939 partition, including the approved diagonal treatment. `serr.py` implements Algorithm 1, `algorithm939.py` provides SE/AF/NI component ports, and `rr.py` provides RR.

## C++

Include `src/cpp/hybrid.hpp` and link against the generated `build/hybrid/hybrid_cpp.dll` using the GNU toolchain (or load its C ABI). The C++ functions are:

```cpp
int marcum::select_method(double a, double b, double nu, bool proposed=true);
double marcum::hybrid_q(double a, double b, double nu, std::int64_t setting,
                       bool proposed=true);
```

The C ABI exports `marcum_cpp_select_method` and `marcum_cpp_hybrid_q`. Keep `measure.dll` alongside `hybrid_cpp.dll` because LAE is bridged to the official Fortran implementation. The build command in the manifest records the exact link inputs.

## Fortran

The `marcum_hybrid` module in `src/fortran/hybrid.f90` exposes:

```fortran
method = select_method(a, b, nu, proposed)
q = hybrid_q(a, b, nu, setting, proposed)
```

Use `real(c_double)` parameters, `integer(c_int64_t)` work counts, and an `integer(c_int)` proposed flag (1 revised, 0 original). The C ABI exports `marcum_fortran_select_method` and `marcum_fortran_hybrid_q`. The component adapter is `src/fortran/algorithm939.f90`. Internal module names in `serr.f90` are retained to preserve the existing callers and numerical kernels.

## Shared conventions

Method IDs are SE=0, AF=1, NI=2, SE-RR=3, RR=4 and LAE=5. The parameter-dependent method choice happens before evaluation. `setting` is a caller-selected fixed truncation/refinement count, not a requested tolerance. The experiment drivers choose this count using the independent reference before timing. SE-RR requires N>=1; all interfaces require finite positive a,b,nu. Scalar wrappers do not share coefficient caches between calls. Fortran components use module state and are called sequentially within a process; experiment concurrency uses separate processes.

SE, AF, NI, SE-RR and RR retain their existing per-language calculations. Python/C++ LAE uses the official Fortran routine explicitly. C++ bound evaluation uses the existing C++ incomplete-gamma primitive. The standalone language wrappers and their source paths are public; language sharing is documented rather than hidden. Paper benchmark kernels, calibrated settings and archived timing records retain their original behavior.
