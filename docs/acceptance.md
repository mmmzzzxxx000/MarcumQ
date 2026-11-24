# Collection acceptance

- Numeric kernels retain their selected source computations; paths and the C++ local header include have been reorganized.
- Python Algorithm 1 function ASTs match the selected implementation.
- Raw Python SE/AF/NI source matches the selected source bytes.
- Extracted private helper function ASTs match their source functions.
- The Fortran random-test numerical source matches the completed experiment source.
- GNU builds completed for C++, Fortran, official 939 and scalar/batch application kernels.
- The original 939 driver passed 18 reference cases (maximum relative error 9.302e-13).
- 65 numerical, reference, batch, region, RR-port and CDF component checks passed.
- All manuscript plotting entry points completed using archived data.
- The 99-case outage driver and the 10,000-triple random driver completed.
- The three cancellation cases reproduced their high-precision finite-sum checks.
- The recurrence-error experiment completed 107 parameter cases and 27,285 exact-rational checks with zero violations.

The full 100-million experiment and all full-plane timings were not rerun during collection. Archived experimental records are preserved separately from acceptance runs. No new acceptance timing is substituted into manuscript data.

## Manuscript synchronization checks

- Three-language builds and all 65 numerical/route/component checks passed in the relocated checkout.
- The official Algorithm 939 driver passed all 18 reference cases.
- The 99-case outage driver completed with both methods attaining relative error 1e-10.
- The random accuracy driver completed a 10,000-triple check with 9,453 admitted and passed triples.
- The new random-timing driver completed a 100-triple end-to-end check, including independent reference construction, calibration, timed batches and candidate remeasurement.
- All six plotting entry points completed from archived data, including the random speedup figure.
- The archive-to-manuscript audit checked the individual-method common sets, hybrid grid means, CDF table, 99-point outage means, random timing speedups, 100-million accuracy totals, and cancellation table.
- The full 10-million timing and 100-million accuracy experiments were not repeated during packaging. Their archived measurements remain unchanged.
- Response-only component diagrams and their temporary scripts are not part of this package.

- All 24 complete CDF curves (GammaCHI/hybrid, pointwise/batch, six configurations), totaling 24,000 values, met relative error 1e-10 at their archived settings.
- The public hybrid wrapper was checked against archived independent references for SE, AF, NI, SE-RR, RR and LAE. All six branches met the target.

## Chapter-based source layout checks

- The three public hybrid source files are present and build in the reorganized checkout.
- 10,588 original/revised routing cases per native selector agreed with the existing region rules; 18 language/method evaluations passed independent Bessel-integral checks.
- All 65 existing numerical/route/component checks passed after relocation.
- All ten manuscript image files and the table/statistic exports were regenerated from the reorganized archive.
- The relocated 99-case outage driver, 100-triple random-timing driver and 10,000-triple random-accuracy driver completed.
- The Section 3 recurrence check and Section 4 cancellation examples completed through their new paths.
- All 288 archived data and third-party files retain their original SHA256 hashes. Existing SE-RR numerical kernels and baseline component calculations are unchanged.
- Response-only component figures remain excluded. Full 10-million and 100-million experiments were not rerun for this directory-only reorganization.
