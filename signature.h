#ifndef _SPENT_H
#define _SPENT_H

#include <stdint.h>
#include "poly_q.h"
#include "param.h"

typedef struct
{
    POLY_Q pk[N_BAR];
} RPK;

typedef struct
{
    POLY_R sk[N_BAR + K];
} SCK;

typedef struct
{
    POLY_R x;
    POLY_R f1[N_SPENT];
    POLY_R za[N_HAT + K_HAT];
    POLY_R z[N_BAR + K];
    POLY_Q E[N_HAT];
    POLY_Q F[N_HAT];
    POLY_QHAT H[N_HAT];
} SIGNATURE_OUT;

void ringsign(SIGNATURE_OUT *out, unsigned int *mu, const RPK rpk_in[N_SPENT], uint64_t l, SCK *ask, const POLY_Q mat1[][N], const POLY_Q mat2[][3], const POLY_QHAT mat3[][N_HAT], unsigned char *pp);

#endif
