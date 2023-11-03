#ifndef _VERIFY_H
#define _VERIFY_H

#include <stdint.h>
#include "poly_q.h"
#include "param.h"
#include "signature.h"

uint64_t verify(const RPK a_in[],  const SIGNATURE_OUT *sn_out, unsigned int *mu, const POLY_Q mat1[][N], const POLY_Q mat2[][3], const POLY_QHAT mat3[][N_HAT], unsigned char *pp);

#endif
