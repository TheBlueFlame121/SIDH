#include "mont_curve.h"
#include <gmp.h>

void MC_init(MC *Res) { fp2_init(&(Res->A)); }

void MP_init(MP *Res) {
    fp2_init(&(Res->X));
    fp2_init(&(Res->Z));
}

void MP_set(MP *Res, const fp2 x, const fp2 z) {
    fp2_copy(&(Res->X), x);
    fp2_copy(&(Res->Z), z);
}

void MP_set_mpz(MP *Res, const mpz_t x1, const mpz_t x0, const mpz_t z1,
                const mpz_t z0) {
    fp2_set_mpz(&(Res->X), x1, x0);
    fp2_set_mpz(&(Res->Z), z1, z0);
}

void MP_set_si(MP *Res, const long int x1, const long int x0, const long int z1,
               const long int z0) {
    fp2_set_si(&(Res->X), x1, x0);
    fp2_set_si(&(Res->Z), z1, z0);
}

void MP_set_x(MP *Res, const fp2 x) {
    fp2_copy(&(Res->X), x);
    fp2_set_si(&(Res->Z), 0, 1);
}

void MP_copy(MP *Res, const MP P) {
    fp2_copy(&(Res->X), P.X);
    fp2_copy(&(Res->Z), P.Z);
}

void MP_clear(MP *Res) {
    fp2_clear(&(Res->X));
    fp2_clear(&(Res->Z));
}

void MP_norm(MP *Res, const MP P) {
    fp2_div(&(Res->X), P.X, P.Z);
    fp2_set_si(&(Res->Z), 0, 1);
}

void xDBL(MP *Res, const MP P, const MP Ap24_C24) {
    fp2_sub(&(global_params.t0), P.X, P.Z);
    fp2_add(&(global_params.t1), P.X, P.Z);
    fp2_sqr(&(global_params.t3), global_params.t0);
    fp2_copy(&(global_params.t0), global_params.t3);

    fp2_sqr(&(global_params.t3), global_params.t1);
    fp2_copy(&(global_params.t1), global_params.t3);
    fp2_mul(&(Res->Z), Ap24_C24.Z, global_params.t0);
    fp2_mul(&(Res->X), Res->Z, global_params.t1);

    fp2_sub(&(global_params.t3), global_params.t1, global_params.t0);
    fp2_copy(&(global_params.t1), global_params.t3);
    fp2_mul(&(global_params.t0), Ap24_C24.X, global_params.t1);
    fp2_add(&(global_params.t3), Res->Z, global_params.t0);
    fp2_copy(&(Res->Z), global_params.t3);

    fp2_mul(&(global_params.t3), Res->Z, global_params.t1);
    fp2_copy(&(Res->Z), global_params.t3);
    return;
}

void xDBLe(MP *Res, const MP P, const MP Ap24_C24, int e) {
    MP_copy(Res, P);
    MP MP_temp;
    MP_init(&MP_temp);
    for (int i = 1; i <= e; i++) {
        xDBL(&MP_temp, *Res, Ap24_C24);
        MP_copy(Res, MP_temp);
    }
    MP_clear(&MP_temp);
    return;
}

// Be careful with the parmams
// In Public Key:= R = P-Q
// In function definition:= R = Q-P
// Might not make a difference as we are omitting y.
void xDBLADD(MP *Dub, MP *Add, const MP P, const MP Q, const MP R,
             const MP Ap24_1) {
    fp2_add(&(global_params.t0), P.X, P.Z);
    fp2_sub(&(global_params.t1), P.X, P.Z);
    fp2_sqr(&(Dub->X), global_params.t0);
    fp2_sub(&(global_params.t2), Q.X, Q.Z);
    fp2_add(&(Add->X), Q.X, Q.Z);
    fp2_mul(&(global_params.t3), global_params.t0, global_params.t2);
    fp2_copy(&(global_params.t0), global_params.t3);
    fp2_sqr(&(Dub->Z), global_params.t1);

    fp2_mul(&(global_params.t3), global_params.t1, Add->X);
    fp2_copy(&(global_params.t1), global_params.t3);
    fp2_sub(&(global_params.t2), Dub->X, Dub->Z);
    fp2_mul(&(global_params.t3), Dub->X, Dub->Z);
    fp2_copy(&(Dub->X), global_params.t3);
    fp2_mul(&(Add->X), Ap24_1.X, global_params.t2);
    fp2_sub(&(Add->Z), global_params.t0, global_params.t1);
    fp2_add(&(global_params.t3), Add->X, Dub->Z);
    fp2_copy(&(Dub->Z), global_params.t3);
    fp2_add(&(Add->X), global_params.t0, global_params.t1);

    fp2_mul(&(global_params.t3), Dub->Z, global_params.t2);
    fp2_copy(&(Dub->Z), global_params.t3);
    fp2_sqr(&(global_params.t3), Add->Z);
    fp2_copy(&(Add->Z), global_params.t3);
    fp2_sqr(&(global_params.t3), Add->X);
    fp2_copy(&(Add->X), global_params.t3);
    fp2_mul(&(global_params.t3), R.X, Add->Z);
    fp2_copy(&(Add->Z), global_params.t3);
    fp2_mul(&(global_params.t3), R.Z, Add->X);
    fp2_copy(&(Add->X), global_params.t3);
    return;
}

void xTPL(MP *Res, const MP P, const MP Ap24_Am24) {
    fp2_sub(&(global_params.t0), P.X, P.Z);
    fp2_sqr(&(global_params.t2), global_params.t0);
    fp2_add(&(global_params.t1), P.X, P.Z);
    fp2_sqr(&(global_params.t3), global_params.t1);
    fp2_add(&(global_params.t4), global_params.t1, global_params.t0);
    fp2_sub(&(global_params.fp2res), global_params.t1, global_params.t0);
    fp2_copy(&(global_params.t0), global_params.fp2res);

    fp2_sqr(&(global_params.t1), global_params.t4);
    fp2_sub(&(global_params.fp2res), global_params.t1, global_params.t3);
    fp2_copy(&(global_params.t1), global_params.fp2res);
    fp2_sub(&(global_params.fp2res), global_params.t1, global_params.t2);
    fp2_copy(&(global_params.t1), global_params.fp2res);
    fp2_mul(&(global_params.t5), global_params.t3, Ap24_Am24.X);
    fp2_mul(&(global_params.fp2res), global_params.t5, global_params.t3);
    fp2_copy(&(global_params.t3), global_params.fp2res);
    fp2_mul(&(global_params.t6), global_params.t2, Ap24_Am24.Z);

    fp2_mul(&(global_params.fp2res), global_params.t2, global_params.t6);
    fp2_copy(&(global_params.t2), global_params.fp2res);
    fp2_sub(&(global_params.fp2res), global_params.t2, global_params.t3);
    fp2_copy(&(global_params.t3), global_params.fp2res);
    fp2_sub(&(global_params.t2), global_params.t5, global_params.t6);
    fp2_mul(&(global_params.fp2res), global_params.t2, global_params.t1);
    fp2_copy(&(global_params.t1), global_params.fp2res);
    fp2_add(&(global_params.t2), global_params.t3, global_params.t1);
    fp2_sqr(&(global_params.fp2res), global_params.t2);
    fp2_copy(&(global_params.t2), global_params.fp2res);

    fp2_mul(&(Res->X), global_params.t2, global_params.t4);
    fp2_sub(&(global_params.fp2res), global_params.t3, global_params.t1);
    fp2_copy(&(global_params.t1), global_params.fp2res);
    fp2_sqr(&(global_params.fp2res), global_params.t1);
    fp2_copy(&(global_params.t1), global_params.fp2res);
    fp2_mul(&(Res->Z), global_params.t1, global_params.t0);
    return;
}

void xTPLe(MP *Res, const MP P, const MP Ap24_Am24, int e) {
    MP_copy(Res, P);
    MP MP_temp;
    MP_init(&MP_temp);
    for (int i = 1; i <= e; i++) {
        xTPL(&MP_temp, *Res, Ap24_Am24);
        MP_copy(Res, MP_temp);
    }
    MP_clear(&MP_temp);
    return;
}

void Ladder3pt(MP *Res, const mpz_t m, const MP P, const MP Q, const MP R,
               const MP A_1) {
    MP p0, p1, p2, Ap24_1, MP_temp0, MP_temp1, MP_temp2;
    MP_init(&p0);
    MP_init(&p1);
    MP_init(&p2);
    MP_init(&Ap24_1);
    MP_init(&MP_temp0);
    MP_init(&MP_temp1);
    MP_init(&MP_temp2);
    MP_copy(&p0, Q);
    MP_copy(&p1, P);
    MP_copy(&p2, R);

    fp2_set_si(&(Ap24_1.Z), 0, 1);
    fp2_add_si(&(Ap24_1.X), A_1.X, 2);
    fp2_set_si(&(global_params.t0), 0, 4);
    fp2_div(&(global_params.fp2res), Ap24_1.X, global_params.t0);
    fp2_copy(&(Ap24_1.X), global_params.fp2res);

    for (long int i = 0; i < mpz_sizeinbase(m, 2); i++) {
        if (mpz_tstbit(m, i) == 1) {
            xDBLADD(&MP_temp0, &MP_temp1, p0, p1, p2, Ap24_1);
            MP_copy(&p0, MP_temp0);
            MP_copy(&p1, MP_temp1);
        } else {
            xDBLADD(&MP_temp0, &MP_temp2, p0, p2, p1, Ap24_1);
            MP_copy(&p0, MP_temp0);
            MP_copy(&p2, MP_temp2);
        }
    }
    MP_copy(Res, p1);
    MP_clear(&p0);
    MP_clear(&p1);
    MP_clear(&p2);
    MP_clear(&Ap24_1);
    MP_clear(&MP_temp0);
    MP_clear(&MP_temp1);
    MP_clear(&MP_temp2);
    return;
}

void jInvariant(fp2 *j, const MP A_C) {
    fp2_sqr(j, A_C.X);

    fp2_sqr(&(global_params.t1), A_C.Z);

    fp2_add(&(global_params.t0), global_params.t1, global_params.t1);

    fp2_sub(&(global_params.fp2res), *j, global_params.t0);
    fp2_copy(&(global_params.t0), global_params.fp2res);

    fp2_sub(&(global_params.fp2res), global_params.t0, global_params.t1);
    fp2_copy(&(global_params.t0), global_params.fp2res);

    fp2_sub(j, global_params.t0, global_params.t1);

    fp2_sqr(&(global_params.fp2res), global_params.t1);
    fp2_copy(&(global_params.t1), global_params.fp2res);

    fp2_mul(&(global_params.fp2res), *j, global_params.t1);
    fp2_copy(j, global_params.fp2res);

    fp2_add(&(global_params.fp2res), global_params.t0, global_params.t0);
    fp2_copy(&(global_params.t0), global_params.fp2res);

    fp2_add(&(global_params.fp2res), global_params.t0, global_params.t0);
    fp2_copy(&(global_params.t0), global_params.fp2res);

    fp2_sqr(&(global_params.t1), global_params.t0);

    fp2_mul(&(global_params.fp2res), global_params.t0, global_params.t1);
    fp2_copy(&(global_params.t0), global_params.fp2res);

    fp2_add(&(global_params.fp2res), global_params.t0, global_params.t0);
    fp2_copy(&(global_params.t0), global_params.fp2res);

    fp2_add(&(global_params.fp2res), global_params.t0, global_params.t0);
    fp2_copy(&(global_params.t0), global_params.fp2res);

    fp2_inv(&(global_params.fp2res), *j);
    fp2_copy(j, global_params.fp2res);

    fp2_mul(&(global_params.fp2res), global_params.t0, *j);
    fp2_copy(j, global_params.fp2res);
    return;
}

void get_A(fp2 *A, const MP P, const MP Q, const MP R) {
    fp2_add(&(global_params.t1), P.X, Q.X);
    fp2_mul(&(global_params.t0), P.X, Q.X);
    fp2_mul(A, R.X, global_params.t1);
    fp2_add(&(global_params.fp2res), *A, global_params.t0);
    fp2_copy(A, global_params.fp2res);

    fp2_mul(&(global_params.fp2res), global_params.t0, R.X);
    fp2_copy(&(global_params.t0), global_params.fp2res);
    fp2_sub_si(&(global_params.fp2res), *A, 1);
    fp2_copy(A, global_params.fp2res);
    fp2_add(&(global_params.fp2res), global_params.t0, global_params.t0);
    fp2_copy(&(global_params.t0), global_params.fp2res);
    fp2_add(&(global_params.fp2res), global_params.t1, R.X);
    fp2_copy(&(global_params.t1), global_params.fp2res);

    fp2_add(&(global_params.fp2res), global_params.t0, global_params.t0);
    fp2_copy(&(global_params.t0), global_params.fp2res);
    fp2_sqr(&(global_params.fp2res), *A);
    fp2_copy(A, global_params.fp2res);
    fp2_inv(&(global_params.fp2res), global_params.t0);
    fp2_copy(&(global_params.t0), global_params.fp2res);
    fp2_mul(&(global_params.fp2res), *A, global_params.t0);
    fp2_copy(A, global_params.fp2res);

    fp2_sub(&(global_params.fp2res), *A, global_params.t1);
    fp2_copy(A, global_params.fp2res);
    return;
}
