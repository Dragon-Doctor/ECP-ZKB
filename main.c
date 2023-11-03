#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "param.h"
#include "poly_q.h"
#include "keygen.h"

#include "sammat.h"
#include "randombytes.h"
#include "fastrandombytes.h"
#include "setup.h"
#include "signature.h"
#include "verify.h"
#include "cpucycles.h"


static unsigned char seed[CRYPTO_BYTES];
static POLY_Q mat1[K + 3][N], mat2[K][3];
static POLY_QHAT mat3[K_HAT + 2 * N_SPENT + 2][N_HAT];


#define NUM_ROUNDS 200

int main() {
    static unsigned int mu[D * QHAT_BYTE];
    uint64_t i, p, r;
    uint64_t l;
    long long cycle1, cycle2, cycle3, cycle4, cycle5;
    static RPK a_in[N_SPENT];
    static SCK ask[N_SPENT];
    static POLY_R sk[N_SPENT][N_BAR + K];
    static unsigned char pp[CRYPTO_BYTES];
    static SIGNATURE_OUT sig_out;

    srand(time(NULL));

    for (i = 0; i < D * 4; i++) {
        mu[i] = rand() % 2;
    }
    FILE *f;
    f = fopen("result.txt", "a");
    if (f == NULL) {
        perror("Error opening file");
        return 1;
    }
    for ( p = 0; p < NUM_ROUNDS; p++) {

        randombytes(seed, CRYPTO_BYTES);
        fastrandombytes_setseed_prv(seed);
        randombytes(seed, CRYPTO_BYTES);
        fastrandombytes_setseed_pub(seed);

        setup(pp, seed);


        cycle1 = cpucycles();
        sample_mat1(mat1);
        sample_mat2(mat2);
        sample_mat3(mat3);
        cycle2 = cpucycles();

        l = rand() % N_SPENT;

        for (i = 0; i < N_SPENT; i++) {
            keygen(a_in[i].pk, sk[i], mat1);
        }
        cycle3 = cpucycles();

        for (i = 0; i < N_BAR + K; i++)
        {
            memcpy((ask[i].sk) + i, sk[i] + i, sizeof(POLY_R));
        }

        ringsign(&sig_out, mu, a_in, l, ask[l].sk, mat1, mat2, mat3, pp);

        cycle4 = cpucycles();

        r = verify(a_in, &sig_out, mu, mat1, mat2, mat3, pp);
        cycle5 = cpucycles();

        //printf("Setup:%lld,  Keygen:%lld, RSign %lld,  Verify: %lld", cycle2 - cycle1, cycle3 - cycle2, cycle4 - cycle3, cycle5 - cycle4);

        printf("    scale:us   Setup:%lld,  Keygen:%lld, RSign %lld,  Verify: %lld,  ACCEPT: %lu \n", (cycle2 - cycle1) / 3700, (cycle3 - cycle2) / 3700, (cycle4 - cycle3) / 3700, (cycle5 - cycle4) / 3700, r);
        fprintf(f,"LEVEL:%lld, Setup:%lld,  Keygen:%lld, RSign %lld,  Verify: %lld,  ACCEPT: %lu \n",N_SPENT, (cycle2 - cycle1) / 3700, (cycle3 - cycle2) / 3700, (cycle4 - cycle3) / 3700, (cycle5 - cycle4) / 3700, r);
    }
    fclose(f);
    return 0;
}
