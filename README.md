# Marcum SE-RR — manuscript revision code

Python, C++ and Fortran implementations and reproducibility programs for the following paper:

**Zhong-Xuan Mao and Jing-Feng Tian.**
*Efficient numerical computation of the generalized Marcum Q-function with applications.*
arXiv:2611.29209 [math.NA], submitted 22 November 2026.
[Read the preprint](https://arxiv.org/abs/2611.29209).

Project authors: **Zhong-Xuan Mao and Jing-Feng Tian**. The source headers, citation files,
and preprint identify the authors and the research associated with this code.
Algorithm 939 and GammaCHI retain their original authorship and licensing notices.

## Cite this work

If you use SE-RR, the proposed hybrid algorithm, or these experimental results,
please cite the paper. Machine-readable references are provided in
[`CITATION.cff`](CITATION.cff) and [`CITATION.bib`](CITATION.bib).
For uses of Algorithm 939 or GammaCHI, please also cite their original papers listed below.
This citation request does not replace or add to the applicable software license terms.

```bibtex
@misc{MaoTian2026MarcumQ,
  author        = {Zhong-Xuan Mao and Jing-Feng Tian},
  title         = {Efficient numerical computation of the generalized {Marcum} {$Q$}-function with applications},
  year          = {2026},
  eprint        = {2611.29209},
  archivePrefix = {arXiv},
  primaryClass  = {math.NA},
  url           = {https://arxiv.org/abs/2611.29209}
}
```


## Requirements and build

Tested platform: Windows x86-64, Python with NumPy/SciPy/mpmath/Matplotlib, GNU C++ and GNU Fortran 13.2.0. Install Python 3.9 or newer for a new environment (the recorded measurement environment used Python 3.7). Both `g++` and `gfortran` must be on PATH. The current native loaders and build scripts target Windows DLLs; other operating systems have not been validated.

Run from this source directory:

```text
python -m pip install -e .
python -B run.py build
python -B run.py test
python -B run.py test-hybrid
python examples/single_point.py
python examples/batch_cdf.py
python examples/hybrid.py
```

The editable installation is intentional: the experiment programs also use source files and third-party sources in this checkout. A standalone Python wheel is not provided.

## Reproduce manuscript figures from recorded data

```text
python -B run.py plot-methods
python -B run.py plot-hybrid
python -B run.py plot-bound
python -B run.py plot-regions
python -B run.py plot-coefficients
python -B run.py plot-random-speedup
```

To generate every manuscript image under its manuscript filename and export the numerical tables, run `python -B run.py paper-figures`. To check the recorded numerical claims without rerunning timings, run `python -B run.py paper-audit`. Outputs are written to `results/figures/`, `results/manuscript_figures/`, and `results/tables/`. The exact ten image-file mappings are listed in `docs/manuscript_coverage.json`. Recorded experimental inputs and results are in `data/`; numerical records are never silently replaced by figure generation.

## Recompute experiments

| Command | Experiment |
|---|---|
| `python -B run.py section51` | Four fixed-setting methods, three languages, orders 1/10/20, 40 by 40 grid |
| `python -B run.py section52` | Original and revised hybrid partitions on that grid |
| `python -B run.py cdf` | 1000 CDF thresholds on [0,5], pointwise and batch, Fortran |
| `python -B run.py outage` | 99 outage cases: K=5,...,35 plus 29.5 and 30.5, with three thresholds |
| `python -B run.py random-timing --count 10000000 --workers 8` | Section 5.2 Fortran random timing on the full manuscript parameter range |
| `python -B run.py random-timing --count 100 --workers 2` | Small end-to-end random-timing check |
| `python -B run.py outage-reference` | Regenerate the 99 independent Bessel-integral references |
| `python -B run.py recurrence` | Exact-rational coefficient-error experiment |
| `python -B run.py cancellation` | Three fixed-N cancellation examples |
| `python -B run.py random --count 1000000 --threads 4` | Smaller run of the Fortran random-point experiment |
| `python -B run.py section53 --count 100000000 --threads 4` | Full random-point experiment |

Section 5.1 reuses the archived independent high-precision Bessel integral reference grid. Section 5.2 reuses its archived external calibration; it recalibrates RR seed work. Application reference programs evaluate the defining Bessel integral at different working precisions. Recomputed outputs go under `results/`. The original published measurements are retained under `data/`.

The chapter aliases `section53`, `section61`, and `section62` run the same drivers as `random`, `cdf`, and `outage`, respectively. The established command names remain available. Configuration files are in `configs/section51` and `configs/section61`.

## Numerical and timing conventions

SE-RR implements Algorithm 1 with the highest index N supplied by the caller (N+1 terms). It contains no reference-based stopping test. The approved equivalent log-prefactor expression is used for orders above 100. Coefficients are prepared before the ascending summation. Direct CDF output and shared batch preparation follow the approved application implementation.

SE/AF/NI follow the corresponding Algorithm 939 components. Python uses SciPy special functions as library primitives. Original default stopping and externally fixed work counts are distinct entry modes. The build records every generated source substitution and preserves the original 939 source files unchanged.

Timing uses the smallest tested setting attaining relative error at most 1e-10, followed by a separate numerical run. References, routing, bound evaluation, setting searches and error checks are outside timers. CDF scalar/batch timings include the corresponding coefficient preparation. The revised hybrid uses SE-RR for b<=7 and B<=1000; both strategies send AF diagonal points to SE. See `docs/experiments.md` for details.

The low-level source kernels do not add clipping, fallback or parameter-validation policies. Callers must supply valid positive parameters and N>=1 for SE-RR. A fixed N alone does not guarantee target accuracy everywhere. `marcum.hybrid.select_method` and `marcum.hybrid.hybrid_q` expose the revised hybrid interfaces, including RR and LAE; the work count is supplied by the caller. They do not introduce a new stopping rule.

## Source layout

```text
MarcumR1_code/
|-- src/
|   |-- python/marcum/{serr.py,algorithm939.py,hybrid.py}
|   |-- cpp/{serr.cpp,algorithm939.cpp,hybrid.cpp,hybrid.hpp}
|   `-- fortran/{serr.f90,algorithm939.f90,hybrid.f90}
|-- experiments/
|   |-- section3/                 Coefficient recurrence/error checks
|   |-- section4/                 Cancellation checks
|   |-- section51/                Individual-method comparisons
|   |-- section52/
|   |   |-- grid/                 Hybrid comparisons on grids
|   |   `-- random/               Random timing comparisons
|   |-- section53/                Random accuracy tests
|   |-- section61/                Noncentral chi CDF
|   |-- section62/                Outage probability
|   `-- common/                   Build, reference and timing support
|-- plots/                        Figure programs by manuscript section
|-- data/                         Recorded results by manuscript section
|-- configs/                      Configuration files by section
|-- tests/                        Numerical and result checks
|-- examples/                     Usage examples
|-- docs/                         Experiment and interface documentation
|-- third_party/                  Original sources and notices
|-- run.py                        Common execution entry point
`-- README.md
```

The brace groups above list the principal entry files, not additional directories. Supporting files remain beside their language implementations. Section 5.2 has matching `grid` and `random` directories in the experiments, data and figure programs; Section 5.3 contains the separate random accuracy test. Shared Bessel reference grids are in `data/section51/references`. Historical metadata inside archived data retains its original filenames and paths.

## Hybrid interfaces

The three public entry files are `src/python/marcum/hybrid.py`, `src/cpp/hybrid.cpp` (declarations in `hybrid.hpp`), and `src/fortran/hybrid.f90` (module `marcum_hybrid`). All use the manuscript's region rules, B<=1000 criterion, M<0.998 branch switch, and approved AF-diagonal handling. `serr` contains Algorithm 1; `algorithm939` contains or adapts the baseline component routines; RR support is kept with the language source. The original upstream sources remain unchanged in `third_party`.

Python exposes `select_method(a,b,nu,proposed=True)` and `hybrid_q(a,b,nu,setting,proposed=True,language='fortran')`. Select `language='python'`, `'cpp'`, or `'fortran'` explicitly, as shown in `examples/hybrid.py`. C++ exposes `marcum::select_method` and `marcum::hybrid_q`; Fortran exposes `select_method` and `hybrid_q`. The `proposed` flag selects the revised strategy or the Algorithm 939 partition. Inputs are finite positive a,b,nu; work counts are externally supplied, with N>=1 for SE-RR. These interfaces add no adaptive stopping rule. GammaCHI CDF batching uses its separate application entry points.

Python and C++ evaluate SE/AF/NI/SE-RR/RR with their existing language implementations. Their LAE branch calls the official Fortran component through an explicit bridge. This shared LAE implementation is not presented as an independent Python or C++ port. The Fortran public wrapper and the random-timing driver share the same fixed-work component adapter. The C++ wrapper also uses the existing incomplete-gamma primitive to evaluate B. See `docs/interfaces.md` for build and call details.

No precompiled executables, runtime DLLs, old launchers or temporary diagnostic scripts are distributed.

## Random timing and archived results

The Section 5.2 timing sample uses seed 2026111702 and independent log-uniform sampling on a=[0.01,100], b=[0.001,100], nu=[1,1000]. It is distinct from the Section 5.3 accuracy sample. The archived timing summaries, batch totals, candidate timings and audits are in `data/section52/random`. All random points and reference values can be regenerated by the supplied driver; the full 10-million-point intermediate arrays are generated under `results/` rather than bundled with the source checkout.

The archived summaries retain all 10,000,000 draws and identify 600,935 points below Algorithm 939's 1e-290 cutoff. Paired timings use 9,399,065 points, including all 7,402,955 points selected for SE-RR. Reference calculation, work-count selection, routing and checks remain outside the numerical timers. The reference uses x87 80-bit Bessel quadrature; when Q is sufficiently close to one, it uses one together with an explicit probability-bound error allowance. These details and classifications are preserved rather than silently relabelled.

All figures can be regenerated from archived numerical records without rerunning either large experiment. Newly measured times depend on the hardware and are never substituted for archived manuscript measurements. Response-only component-timing diagrams are outside this release.

## Attribution and publication

Algorithm 939 is by Amparo Gil, Javier Segura and Nico M. Temme, ACM TOMS 40(3), Article 20 (2014), DOI 10.1145/2591004. Its source is https://calgo.acm.org/939.zip. GammaCHI is by the same authors; its original README is included. Upstream notices remain applicable to upstream code and derived ports. See `third_party/NOTICE.md`.

The authors' paper is publicly available at https://arxiv.org/abs/2611.29209.
The authors' original SE-RR code and independently authored contributions use
**BSD-3-Clause**; see [`LICENSE.md`](LICENSE.md). Algorithm 939 and its ports/adapters
retain the ACM terms, and the bundled 2021 GammaCHI source retains GPLv2.
Those licenses are included under `third_party/`; see
[`docs/licensing.md`](docs/licensing.md) for scope and component details.
Please cite the paper when using our code in academic work.
