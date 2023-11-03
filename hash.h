#ifndef _HASH_H
#define _HASH_H

#include "poly_q.h"
#include "signature.h"
#include "param.h"

void hash_x(POLY_R *out, unsigned int *mu, const POLY_QHAT *E, const POLY_QHAT *F, const POLY_Q *H);

#endif
