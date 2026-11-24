// Project authors: Zhong-Xuan Mao and Jing-Feng Tian.
// Paper: Efficient numerical computation of the generalized Marcum
// Q-function with applications (2026).
// Preprint: https://arxiv.org/abs/2611.29209
// Citation: CITATION.cff / CITATION.bib at the repository root.
// Licensing and upstream attribution: LICENSE.md / third_party/NOTICE.md.

/* Scalar same-language component benchmark ABI. Methods: SE=0,AF=1,NI=2,SERR=3,
 * original-route baseline939=4, SERR-first hybrid939=5. Target: P=0,Q=1.
 * Each library implements all four methods in its own language.
 * Status: 0 numerical criterion met,1 limit,2 precision,3 underflow,
 * 4 unsupported,5 numerical range,6 exact boundary,7 bounded approximation.
 */
#pragma once
#include <stdint.h>
#ifdef __cplusplus
extern "C" {
#endif
typedef struct {
    double value,p,q,log_value,truncation_bound,roundoff_proxy,condition,sum_abs;
    int64_t iterations,evaluations;
    int32_t status,reserved;
} components_result;
int components_abi_version_c(void);
int components_result_size_c(void);
void components_eval_c(int method,double a,double b,double nu,int target,
                      double atol,double rtol,int64_t max_terms,
                      int source_defaults,int max_refinements,components_result *out);
void components_time_c(int method,double a,double b,double nu,int target,
                      double atol,double rtol,int64_t max_terms,
                      int source_defaults,int max_refinements,int64_t repeats,double *seconds);
/* CDF batch plans. 0..5 keep scalar methods; 6=SE cached, 7=SERR cached,
 * 8=cached hybrid. Values returned in component_result are direct P.
 * Caller owns plan, input and output buffers. No allocator crosses the ABI.
 * Timing mode: 0 apply an existing plan, 1 prepare+apply, 2 preparation only. */
typedef struct {
    double a,nu,atol,rtol,bmin,bmax;
    int64_t max_terms,method,se_count,serr_count,status;
    double weights[513],coefficients[161];
} components_batch_plan;
int components_batch_plan_size_c(void);
void components_batch_prepare_c(int method,double a,double nu,double atol,double rtol,
                               int64_t max_terms,double bmin,double bmax,components_batch_plan *plan);
void components_batch_apply_c(const components_batch_plan *plan,int64_t count,
                             const double *b,components_result *out);
void components_batch_time_c(const components_batch_plan *plan,int64_t count,
                            const double *b,int64_t repeats,int mode,double *seconds);
/* Streaming validation of varying triples; calls the SAME scalar kernel. */
void components_triples_c(int method,int64_t count,const double *a,const double *b,
                         const double *nu,int target,double atol,double rtol,
                         int64_t max_terms,components_result *out);
#ifdef __cplusplus
}
#endif
