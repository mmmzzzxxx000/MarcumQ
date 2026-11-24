// Project authors: Zhong-Xuan Mao and Jing-Feng Tian.
// Paper: Efficient numerical computation of the generalized Marcum
// Q-function with applications (2026).
// Preprint: https://arxiv.org/abs/2611.29209
// Citation: CITATION.cff / CITATION.bib at the repository root.
// Licensing and upstream attribution: LICENSE.md / third_party/NOTICE.md.

/* SE-RR native-port interface reservation ONLY; no SE-RR native kernels yet.
 * The implemented original-939 adapter has a separate algorithm939.h ABI.
 * Names/layout form a draft v1 contract for later C++ / Fortran BIND(C) ports.
 * No exceptions, Python objects, allocation ownership, or string objects cross ABI.
 */
#ifndef MARCUM_SE_RR_H
#define MARCUM_SE_RR_H
#include <stdint.h>
#ifdef __cplusplus
extern "C" {
#endif

enum { MARCUM_TARGET_P = 0, MARCUM_TARGET_Q = 1 };
enum {
    MARCUM_EXACT = 0, MARCUM_BOUNDED_APPROXIMATION = 1,
    MARCUM_TRUNCATION_CONVERGED = 2, MARCUM_SPECIAL_CASE = 3,
    MARCUM_NATIVE_COMPLETED = 4,
    MARCUM_MAX_TERMS = 10, MARCUM_PRECISION_LIMITED = 11,
    MARCUM_UNDERFLOW = 12, MARCUM_NUMERICAL_RANGE = 13, MARCUM_NATIVE_ERROR = 14,
    MARCUM_INVALID_ARGUMENT = 20, MARCUM_UNAVAILABLE = 21,
    MARCUM_UNSUPPORTED_DOMAIN = 22
};

typedef struct {
    int32_t abi_version; /* must be 1; signed type also maps to Fortran c_int32_t */
    int32_t struct_size;
    double atol;
    double rtol;
    int64_t max_terms;    /* includes n=0 */
} marcum_options_v1;

typedef struct {
    int32_t abi_version;
    int32_t struct_size;
    int32_t status;
    int32_t target;
    int64_t terms;
    int64_t last_index;
    double value;
    double log_value;
    double truncation_bound;
    double log_truncation_bound;
    double summation_roundoff_estimate;
    double condition_estimate;
    double kappa_q_bound;
    double log_envelope_bound;
    int32_t truncation_satisfied;
    int32_t scaling_steps;
} marcum_result_v1;

/* All inputs and returned values use the same a,b,nu convention as Python.
 * Missing diagnostics are NaN. Status, never clipping, signals failure.
 * Implementations return status and also write it in result->status.
 */
int32_t marcum_eval_v1(double a, double b, double nu, int32_t target,
                     const marcum_options_v1 *options, marcum_result_v1 *result);

/* Future batch entry: length-n contiguous arrays, one result per triple.
 * This does not yet specify the separate coefficient-reuse batch plan.
 */
int32_t marcum_eval_batch_v1(int64_t n, const double *a, const double *b,
                           const double *nu, int32_t target,
                           const marcum_options_v1 *options,
                           marcum_result_v1 *results);
#ifdef __cplusplus
}
#endif
#endif
