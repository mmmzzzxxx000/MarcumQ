# Release readiness review (24 November 2026)

## Paper and authorship

Paper: Zhong-Xuan Mao and Jing-Feng Tian, *Efficient numerical computation of the
generalized Marcum Q-function with applications*, arXiv:2611.29209 [math.NA].
Submitted 22 November 2026: https://arxiv.org/abs/2611.29209 .

The README, CITATION.cff, CITATION.bib, package metadata and 83 project source
files identify the paper and project authors. CITATION.cff was validated against
the official CFF 1.2.0 JSON schema. Upstream author comments remain intact.

## Technical checks passed

Checks ran in an isolated source copy on the documented Windows environment.

- Python/C++/Fortran build and the upstream Algorithm 939 driver: passed.
- Numerical acceptance tests: 65 checks passed.
- Hybrid interfaces: 10,588 route triples per native selector and 18 six-method
  language evaluations passed.
- Archived manuscript-result audit: passed for Sections 4, 5.1, 5.2, 5.3, 6.1
  and 6.2, including the reported speedups and sample counts.
- All ten manuscript figures and the numerical tables were regenerated.
- No executable statements were changed by the paper-attribution update.
  Python ASTs match; C++/Fortran changes are comment prefixes only.
- Archived experimental data and bundled upstream source files are byte-identical.
- The release file list contains no compiled binaries or bytecode. A credential
  pattern scan of source/documentation found no embedded private keys or tokens.

The large random experiments were not rerun. Their archived records were audited.
Linux/macOS were not tested; the README describes the Windows native loaders.

## License preparation (24 November 2026)

The authors selected BSD-3-Clause for original SE-RR code and independently
authored contributions. The root LICENSE.md states its scope; the three original
SE-RR kernels carry SPDX notices. ACM and GPLv2 license texts are included under
third_party, and adapted Algorithm 939 components are identified as modified.
No executable statements were changed by the license update.

The ACM agreement requires notifying ACM and the original authors of modified
distributions. A notification draft is kept with the local release preparation;
it has not been sent by the assistant. See docs/licensing.md for component scope.
The actual public repository URL should replace the response-letter placeholder
after publication. A public upload must be verified before claiming release.
