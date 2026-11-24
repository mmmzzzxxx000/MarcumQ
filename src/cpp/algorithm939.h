// Project authors: Zhong-Xuan Mao and Jing-Feng Tian.
// Paper: Efficient numerical computation of the generalized Marcum
// Q-function with applications (2026).
// Preprint: https://arxiv.org/abs/2611.29209
// Citation: CITATION.cff / CITATION.bib at the repository root.
// Licensing and upstream attribution: LICENSE.md / third_party/NOTICE.md.

/* Implemented thin adapter to the unmodified official Fortran Algorithm 939.
 * Uses ORIGINAL variables mu=nu, x=a*a/2, y=b*b/2, NOT public physical a,b.
 * ierr=0: native completed; 1: upstream underflow; 2: unsupported domain.
 * Always validate returned p,q; upstream can return NaN even with ierr=0.
 */
#ifndef MARCUM_ALGORITHM939_H
#define MARCUM_ALGORITHM939_H
#ifdef __cplusplus
extern "C" {
#endif
int algorithm939_abi_version_c(void);
void algorithm939_eval_c(double mu, double x, double y,
                         double *p, double *q, int *ierr);
#ifdef __cplusplus
}
#endif
#endif
