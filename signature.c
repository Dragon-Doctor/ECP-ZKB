#include <string.h>

#include "signature.h"
#include "poly_q.h"
#include "param.h"
#include "poly_red.h"
#include "sampleb.h"
#include "hash.h"
#include "poly_mult.h"
#include "comp.h"

#include "gaussian_avx.h"


void ringsign(SIGNATURE_OUT *out, unsigned int *mu, const RPK rpk_in[], uint64_t l, SCK *ask,
              const POLY_Q mat1[][N], const POLY_Q mat2[][3], const POLY_QHAT mat3[][N_HAT], unsigned char *pp) {
    uint64_t i, j, k, rej;

    static POLY_R a[N_SPENT];
    static POLY_QHAT a_ntt_qhat[N_SPENT];
    static POLY_Q a_ntt[N_SPENT];
    static POLY_QHAT a2_ntt[N_SPENT];
    static uint64_t delta[N_SPENT];
    static POLY_R rb[N_HAT + K_HAT], ra[N_HAT + K_HAT];
    static POLY_QHAT rb_ntt[N_HAT + K_HAT], ra_ntt[N_HAT + K_HAT];
    static POLY_R y[N_BAR + K];
    static POLY_Q y_ntt[N_BAR + K];


    static POLY_R g[N_SPENT];
    static POLY_R f0;
    static POLY_R x_delta[N_SPENT + 1];
    static POLY_R x_s[N_BAR + K];
    static POLY_R x_rb[N_HAT + K_HAT];


    /* b = (\delta_{l, 0},...,\delta_{l, N_SPENT - 1}) */
    for (i = 0; i < N_SPENT; i++) {
        delta[i] = ct_eq(i, l);
    }

    do {
        rej = 0;

        /* a_1,...,a_{N_SPENT - 1} <-- D_{10 * B_a} */
        sample_ba_n1(a - 1);
        for (i = 1; i < N_SPENT; i++) {
            for (j = 0; j < D; j++) {
                a_ntt[i].poly[j] = con_add(a[i].poly[j], Q);
                a_ntt_qhat[i].poly[j] = con_add(a[i].poly[j], QHAT);
            }
            ntt_q(a_ntt + i);
            ntt_qhat(a_ntt_qhat + i);
        }

        /* a_0 = -\sum_{j = 1}^{N_SPENT - 1} a_j */
        memset(a, 0, sizeof(POLY_R));
        for (i = 1; i < N_SPENT; i++)
        {
            for (j = 0; j < D; j++)
            {
                a[0].poly[j] -= a[i].poly[j];
            }
        }
        for (j = 0; j < D; j++)
        {
            a_ntt[0].poly[j] = con_add(a[0].poly[j], Q);
            a_ntt_qhat[0].poly[j] = con_add(a[0].poly[j], QHAT);
        }
        ntt_q(a_ntt);
        ntt_qhat(a_ntt_qhat);

        /* (a_0^2,...,a_{N_SPENT + 1}^2) */
        for (i = 0; i < N_SPENT; i++) {
            mult_rqhat(a2_ntt + i, a_ntt_qhat + i, a_ntt_qhat + i);
        }

        /* rb <-- S_1^{\hat{n} + \hat{k}} */
        /* ra <-- D_B^{\hat{n} + \hat{k}} */
        sample_1(rb, N_HAT + K_HAT);
        sample_b_nhatkhat(ra);
        for (i = 0; i < N_HAT + K_HAT; i++) {
            for (j = 0; j < D; j++) {
                rb_ntt[i].poly[j] = con_add(rb[i].poly[j], QHAT);
                ra_ntt[i].poly[j] = con_add(ra[i].poly[j], QHAT);
            }

            ntt_qhat(rb_ntt + i);
            ntt_qhat(ra_ntt + i);
        }

        /* E = G' * (r_b, b, c) */
        /* F = G' * (r_a, a, d) */

        for (i = 0; i < N_HAT; i++) {
            memcpy((out->E) + i, rb_ntt + i, sizeof(POLY_QHAT));
            memcpy((out->F) + i, ra_ntt + i, sizeof(POLY_QHAT));
        }
        for (i = 0; i < K_HAT; i++) {
            for (j = 0; j < N_HAT; j++) {
                mult_plus_rqhat((out->E) + j, mat3[i] + j, rb_ntt + N_HAT + i);
                mult_plus_rqhat((out->F) + j, mat3[i] + j, ra_ntt + N_HAT + i);
            }
        }
        for (i = 0; i < N_SPENT; i++) {
            for (j = 0; j < N_HAT; j++) {
                for (k = 0; k < D; k++) {
                    (out->E)[j].poly[k] += (-delta[i]) & mat3[K + i][j].poly[k];
                }

                mult_plus_rqhat((out->F) + j, mat3[K + i] + j, a_ntt_qhat + i);
            }
        }
        for (i = 0; i < N_SPENT; i++) {
            for (j = 0; j < N_HAT; j++) {
                mult_plus_rqhat_pm((out->E) + j, mat3[K + N_SPENT  + i] + j, a_ntt_qhat + i, delta[i]);
                mult_minus_rqhat((out->F) + j, mat3[K + N_SPENT  + i] + j, a2_ntt + i);
            }
        }

        /* y <-- D_{13 * B_{rs}}^{\bar{n} + k} */
        sample_brs_13_nbark(y);
        for (i = 0; i < N_BAR + K; i++) {
            for (j = 0; j < D; j++) {
                y_ntt[i].poly[j] = con_add(y[i].poly[j], Q);
            }
            ntt_q(y_ntt + i);
        }

        /* H = G * y - \sum_{j = 0}^{N_SPENT - 1} a_j * P_j */
        for (i = 0; i < N_BAR; i++) {
            memcpy(out->H + i, y_ntt + i, sizeof(POLY_Q));
        }
        for (i = 0; i < K; i++) {
            for (j = 0; j < N_BAR; j++) {
                mult_plus_rq(out->H  + j, mat1[i] + j + (N - N_BAR), y_ntt + N_BAR + i);
            }
        }
        for (i = 0; i < N_SPENT; i++) {
            for (j = 0; j < N_BAR; j++) {
                mult_minus_rq(out->H  + j, a_ntt + i, rpk_in[i].pk + j);
            }
        }

        /* x <-- H(, mu, E, F, H) */
        hash_x(&(out->x), mu, out->E, out->F, out->H);

        /* f_j = x * b_j + a_j */
        for (i = 1; i < N_SPENT; i++) {
            for (j = 0; j < D; j++) {
                x_delta[i - 1].poly[j] = (-delta[i]) & (out->x).poly[j];
                (out->f1)[i - 1].poly[j] = x_delta[i - 1].poly[j] + a[i].poly[j];
            }
        }
        /* Rej(f_1, x * delta, \phi_a, B_a) */
        rej = rej_f(out->f1, x_delta);
        if (rej) {
           continue;
        }

        /* f_0 */
        for (i = 0; i < D; i++) {
            f0.poly[i] = ((-delta[0]) & (out->x).poly[i]) + a[0].poly[i];
        }

        /* g_0 = f_0(x - f_0) */
        mult_fxf(g, &f0, &(out->x));

        /* check g_0 norm here */
        for (j = 0; j < D; j++) {
            rej |= ct_lt(BG_0, g[0].poly[j]) | ct_lt(g[0].poly[j], -BG_0);
        }
        if (rej) {
            continue;
        }

        /* g_1 = (f_1(x - f_1),...,f_{N_SPENT - 1}(X - f_{N_SPENT - 1})) */
        for (i = 1; i < N_SPENT; i++) {
            mult_fxf(g + i, (out->f1) + i - 1, &(out->x));
        }

        /* check g_1 norm here */
        for (i = 1; i < N_SPENT; i++) {
            for (j = 0; j < D; j++) {
                rej |= ct_lt(BG_1, g[i].poly[j]) | ct_lt(g[i].poly[j], -BG_1);
            }
        }
        if (rej) {
            continue;
        }

        /* z = y + x * w (sk) */
        for (i = 0; i < N_BAR + K; i++) {
            mult_r(x_s + i, &(out->x),  ask->sk + i);

            for (j = 0; j < D; j++) {
                (out->z)[i].poly[j] = x_s[i].poly[j] + y[i].poly[j];
            }
        }

        /* Rej(z, x * s, \tau_rs, T_rs) */
        rej = rej_z(out->z, x_s);
        if (rej) {
            continue;
        }

        /* z_a = r_a + x * r_b */
        for (i = 0; i < N_HAT + K_HAT; i++) {
            mult_r(x_rb + i, &(out->x), rb + i);
            for (j = 0; j < D; j++) {
                (out->za)[i].poly[j] = x_rb[i].poly[j] + ra[i].poly[j];
            }
        }
    } while (rej);
}