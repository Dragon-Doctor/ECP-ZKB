#include <string.h>

#include "verify.h"
#include "poly_q.h"
#include "hash.h"
#include "poly_mult.h"
#include "poly_red.h"
#include "signature.h"
#include "param.h"


#define Z_NORM 3161803581581766LL


uint64_t verify(const RPK a_in[], const SIGNATURE_OUT *in, unsigned int *mu, const POLY_Q mat1[][N],
                const POLY_Q mat2[][3], const POLY_QHAT mat3[][N_HAT], unsigned char *pp) {
    uint64_t i, j;
    uint64_t tmp;

    static POLY_R f0;
    static POLY_R g[N_SPENT];

    static POLY_Q x_ntt;
    static POLY_QHAT x_ntt_qhat;

    static POLY_Q f_ntt[N_SPENT];
    static POLY_QHAT f_ntt_qhat[N_SPENT];

    static POLY_Q tmp_aut;

    static POLY_Q tmp_q;
    static POLY_QHAT a[N_HAT];
    static POLY_QHAT zb_ntt[N_HAT + K_HAT];
    static POLY_QHAT g_ntt[N_SPENT];

    static POLY_Q e[N_BAR];
    static POLY_Q z_ntt[N_BAR + K];
    static POLY_R x;
    static POLY_R tmp_aut_r;

    /* check z_a norm */
    for (i = 0; i < N_HAT + K_HAT; i++) {
        for (j = 0; j < D; j++) {
            if (((int64_t) ((in->za)[i].poly[j]) > B_6) || ((int64_t) ((in->za)[i].poly[j]) < -B_6)) {
                return 0;
            }
        }
    }

    /* check f_1 norm */
    for (i = 0; i < N_SPENT ; i++) {
        for (j = 0; j < D; j++) {
            if (((int64_t) ((in->f1)[i].poly[j]) > BA_10_6) || ((int64_t) ((in->f1)[i].poly[j]) < -BA_10_6)) {
                return 0;
            }
        }
    }

    /* check z norm */
    tmp = 0;
    for (i = 0; i < N_BAR + K; i++) {
        for (j = 0; j < D; j++) {
            tmp += (int64_t) ((in->z)[i].poly[j]) * (int64_t) ((in->z)[i].poly[j]);
            if (((int64_t) ((in->z)[i].poly[j]) > BRS_13_6) || ((int64_t) ((in->z)[i].poly[j]) < -BRS_13_6)) {
                return 0;
            }
        }
    }
    if (tmp > Z_NORM) {
        return 0;
    }

    /* f_0 = x - \sum_{j = 0}^{N_SPENT - 1} f_j */
    memcpy(&f0, &(in->x), sizeof(POLY_Q));
    for (i = 0; i < N_SPENT - 1; i++) {
        for (j = 0; j < D; j++) {
            f0.poly[j] -= (in->f1)[i].poly[j];
        }
    }

    /* g_0 = f_0(x - f_0) */
    mult_fxf(g, &f0, &(in->x));

    /* check g_0 norm */
    for (j = 0; j < D; j++) {
        if ((((int64_t) g[0].poly[j]) > BG_0) || (((int64_t) g[0].poly[j]) < -BG_0)) {
            return 0;
        }
    }

    /* g_1 = (f_1(x - f_1),...,f_{N_SPENT - 1}(X - f_{N_SPENT - 1})) */
    for (i = 1; i < N_SPENT ; i++) {
        mult_fxf(g + i, (in->f1) + i - 1, &(in->x));
    }

    /* check g_1 norm */
    for (i = 1; i < N_SPENT; i++) {
        for (j = 0; j < D; j++) {
            if ((((int64_t) g[i].poly[j]) > BG_1) || (((int64_t) g[i].poly[j]) < -BG_1)) {
                return 0;
            }
        }
    }

    /* \sigma^{-1}(x) */
    iaut_r(&tmp_aut_r, &(in->x));

    /* ntt(x) */
    for (i = 0; i < D; i++) {
        x_ntt.poly[i] = con_add((in->x).poly[i], Q);
        x_ntt_qhat.poly[i] = con_add((in->x).poly[i], QHAT);
        tmp_aut.poly[i] = con_add(tmp_aut_r.poly[i], Q);
    }
    ntt_q(&x_ntt);
    ntt_qhat(&x_ntt_qhat);
    ntt_q(&tmp_aut);



    /* ntt(f) */
    for (i = 0; i < D; i++) {
        f_ntt[0].poly[i] = con_add(f0.poly[i], Q);
        f_ntt_qhat[0].poly[i] = con_add(f0.poly[i], QHAT);
    }
    ntt_q(f_ntt);
    ntt_qhat(f_ntt_qhat);
    for (i = 1; i < N_SPENT; i++) {
        for (j = 0; j < D; j++) {
            f_ntt[i].poly[j] = con_add((in->f1)[i - 1].poly[j], Q);
            f_ntt_qhat[i].poly[j] = con_add((in->f1)[i - 1].poly[j], QHAT);
        }

        ntt_q(f_ntt + i);
        ntt_qhat(f_ntt_qhat + i);
    }


    /* ntt(z_b) */
    for (i = 0; i < N_HAT + K_HAT; i++) {
        for (j = 0; j < D; j++) {
            zb_ntt[i].poly[j] = con_add((in->za)[i].poly[j], QHAT);
        }

        ntt_qhat(zb_ntt + i);
    }

    /* ntt(g) */
    for (i = 0; i < N_SPENT ; i++) {
        for (j = 0; j < D; j++) {
            tmp_q.poly[j] = con_add(x_ntt_qhat.poly[j] - f_ntt_qhat[i].poly[j], QHAT);
        }

        mult_rqhat(g_ntt + i, f_ntt_qhat + i, &tmp_q);
    }

    /* F = G' * (z_b, f, g) - x * E */
    for (i = 0; i < N_HAT; i++) {
        memcpy(a + i, zb_ntt + i, sizeof(POLY_QHAT));
    }
    for (i = 0; i < K_HAT; i++) {
        for (j = 0; j < N_HAT; j++) {
            mult_plus_rqhat(a + j, mat3[i] + j, zb_ntt + N_HAT + i);
        }
    }
    for (i = 0; i < N_SPENT ; i++) {
        for (j = 0; j < N_HAT; j++) {
            mult_plus_rqhat(a + j, mat3[K + i] + j, f_ntt_qhat + i);
        }
    }
    for (i = 0; i < N_SPENT; i++) {
        for (j = 0; j < N_HAT; j++) {
            mult_plus_rqhat(a + j, mat3[K + N_SPENT + 2 + i] + j, g_ntt + i);
        }
    }
    for (i = 0; i < N_HAT; i++) {
        mult_minus_rqhat(a + i, &x_ntt_qhat, (in->E) + i);
    }


    /* ntt(z) */
    for (i = 0; i < N_BAR + K; i++) {
        for (j = 0; j < D; j++) {
            z_ntt[i].poly[j] = con_add((in->z)[i].poly[j], Q);
        }

        ntt_q(z_ntt + i);
    }

    /* E = G * z - \sum_{j = 0}^{N_SPENT - 1} (f_j * P_j) */
    for (i = 0; i < N_BAR; i++) {
        memcpy(e + i, z_ntt + i, sizeof(POLY_Q));
    }
    for (i = 0; i < K; i++) {
        for (j = 0; j < N_BAR; j++) {
            mult_plus_rq(e + j, mat1[i] + j + (N - N_BAR), z_ntt + N_BAR + i);
        }
    }
    for (i = 0; i < N_SPENT; i++) {
        for (j = 0; j < N_BAR; j++) {
            mult_minus_rq(e + j, f_ntt + i, a_in[i].pk + j);
        }
    }

    /* x <-- H(\alpha, v_0, v_1, A, B, E) */
    hash_x(&x, mu, in->E, in->F, in->H);
    for (i = 0; i < D; i++) {
        if (x.poly[i] != (in->x).poly[i]) {
            return 0;
        }
    }

    return 1;
}
