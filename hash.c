#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "hash.h"
#include "poly_q.h"
#include "param.h"
#include "signature.h"
#include "littleendian.h"

#include "reference/libXKCP.a.headers/SimpleFIPS202.h"

#define HASH_ALPHA_LEN (CRYPTO_BYTES + N_SPENT * M * (N_BAR * D * Q_BYTE + (N + 1) * D * Q_BYTE) + M * CRYPTO_BYTES + S * N_BAR * D * Q_BYTE + S * (N + 1) * D * Q_BYTE + S * N * D * Q_BYTE + N * D * Q_BYTE + N * D * Q_BYTE + D * Q_BYTE + D * Q_BYTE + D * Q_BYTE + N * D * Q_BYTE + N * D * Q_BYTE)

#define HASH_X_INPUT_LEN (  D * Q_BYTE  + 3 * N_HAT * D * QHAT_BYTE )
#define HASH_X_OUTPUT_LEN 137LL
#define HASH_X_COEFF_LEN 6LL
#define HASH_X_GROUP 4LL
#define HASH_X_W (W / HASH_X_GROUP)
#define HASH_X_MASK (D / HASH_X_GROUP - 1)


void hash_x(POLY_R *out, unsigned int *mu, const POLY_QHAT *E, const POLY_QHAT *F, const POLY_Q *H)
{


    static unsigned char hash_input[HASH_X_INPUT_LEN];
    unsigned char hash_output[HASH_X_OUTPUT_LEN];

    unsigned char *hash_input_head;

    uint64_t i, j, boo, x_pos, k;
    uint64_t coeff, count = 0;

    unsigned char *hash_output_head = hash_output + HASH_X_COEFF_LEN;
    uint64_t nonzero_pos[HASH_X_GROUP][HASH_X_W];


    /* mu */
    memcpy(hash_input,mu,(D*QHAT_BYTE));
    hash_input_head = hash_input + D * Q_BYTE;


    /* E */
    for (i = 0; i < N_HAT; i++)
    {
        for (j = 0; j < D; j++)
        {
            STORE_QHAT(hash_input_head, E[i].poly[j]);
            hash_input_head += QHAT_BYTE;
        }
    }

    /* F */
    for (i = 0; i < N_HAT; i++)
    {
        for (j = 0; j < D; j++)
        {
            STORE_QHAT(hash_input_head, F[i].poly[j]);
            hash_input_head += QHAT_BYTE;
        }
    }

    /* H */
    for (i = 0; i < N_BAR; i++)
    {
        for (j = 0; j < D; j++)
        {
            STORE_Q(hash_input_head, H[i].poly[j]);
            hash_input_head += Q_BYTE;
        }
    }

    SHAKE256(hash_output, HASH_X_OUTPUT_LEN, hash_input, HASH_X_INPUT_LEN);

    /* sample and fill the nonzero positions from the hash output
     * since it may generates duplicate positions, we need to do a rejection sampling here */
    memset(out, 0, sizeof(POLY_R));
    coeff = load_48(hash_output);
    for (i = 0; i < HASH_X_GROUP; i++)
    {
        nonzero_pos[i][0] = (*(hash_output_head++)) & HASH_X_MASK;
        out->poly[nonzero_pos[i][0] * HASH_X_GROUP + i] = 1 - 2 * ((coeff >> (count++)) & 0x1);
    }
    for (k = 0; k < HASH_X_GROUP; k++)
    {
        for (i = 1; i < HASH_X_W; i++)
        {
            do
            {
                boo = 0;
                x_pos = (*(hash_output_head++)) & HASH_X_MASK;

                for (j = 0; j < i; j++)
                {
                    if (x_pos == nonzero_pos[k][j])
                    {
                        boo = 1;
                        break;
                    }
                }
            } while (boo);

            nonzero_pos[k][i] = x_pos;
            out->poly[x_pos * HASH_X_GROUP + k] = 1 - 2 * ((coeff >> (count++)) & 0x1);
        }
    }
}
