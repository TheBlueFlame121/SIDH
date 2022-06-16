#include "sidh.h"
#include <gmp.h>

void isogen2(fp2 *x1, fp2 *x2, fp2 *x3, mpz_t sk2, int e2, int e3, mpz_t p,
             fp2 xP2, fp2 xQ2, fp2 xR2, fp2 xP3, fp2 xQ3, fp2 xR3) {
    MP A_C, Ap24_C24, S, P1, P2, P3, P, Q, R, test1;
    MP_init(&A_C);
    MP_init(&Ap24_C24);
    MP_init(&S);
    MP_init(&P1);
    MP_init(&P2);
    MP_init(&P3);
    MP_init(&P);
    MP_init(&Q);
    MP_init(&R);
    MP_init(&test1);

    MP_set_si(&A_C, 0, 6, 0, 1);
    MP_set_si(&Ap24_C24, 0, 8, 0, 4);

    MP_set_x(&P, xP2);
    MP_set_x(&Q, xQ2);
    MP_set_x(&R, xR2);

    MP_set_x(&P1, xP3);
    MP_set_x(&P2, xQ3);
    MP_set_x(&P3, xR3);

    Ladder3pt(&S, sk2, P, Q, R, A_C);

    e_2_iso_optional(&test1, &P, &Q, &R, Ap24_C24, S, P1, P2, P3, e2);
    MP_copy(&Ap24_C24, test1);

    MP_norm(&P1, P);
    MP_norm(&P2, Q);
    MP_norm(&P3, R);

    fp2_copy(x1, P1.X);
    fp2_copy(x2, P2.X);
    fp2_copy(x3, P3.X);

    MP_clear(&A_C);
    MP_clear(&Ap24_C24);
    MP_clear(&S);
    MP_clear(&P1);
    MP_clear(&P2);
    MP_clear(&P3);
    MP_clear(&P);
    MP_clear(&Q);
    MP_clear(&R);
    MP_clear(&test1);

    return;
}

void isogen3(fp2 *x1, fp2 *x2, fp2 *x3, mpz_t sk3, int e2, int e3, mpz_t p, fp2 xP2, fp2 xQ2, fp2 xR2, fp2 xP3, fp2 xQ3, fp2 xR3) {
    MP A_C, Ap24_Am24, S, P1, P2, P3, P, Q, R, test1;
    MP_init(&A_C);
    MP_init(&Ap24_Am24);
    MP_init(&S);
    MP_init(&P1);
    MP_init(&P2);
    MP_init(&P3);
    MP_init(&P);
    MP_init(&Q);
    MP_init(&R);
    MP_init(&test1);

    MP_set_si(&A_C, 0, 6, 0, 1);
    MP_set_si(&Ap24_Am24, 0, 8, 0, 4);

    MP_set_x(&P, xP3);
    MP_set_x(&Q, xQ3);
    MP_set_x(&R, xR3);

    MP_set_x(&P1, xP2);
    MP_set_x(&P2, xQ2);
    MP_set_x(&P3, xR2);

    Ladder3pt(&S, sk3, P, Q, R, A_C);

    e_3_iso_optional(&test1, &P, &Q, &R, Ap24_Am24, S, P1, P2, P3, e3);
    MP_copy(&Ap24_Am24, test1);

    MP_norm(&P1, P);
    MP_norm(&P2, Q);
    MP_norm(&P3, R);

    fp2_copy(x1, P1.X);
    fp2_copy(x2, P2.X);
    fp2_copy(x3, P3.X);

    MP_clear(&A_C);
    MP_clear(&Ap24_Am24);
    MP_clear(&S);
    MP_clear(&P1);
    MP_clear(&P2);
    MP_clear(&P3);
    MP_clear(&P);
    MP_clear(&Q);
    MP_clear(&R);
    MP_clear(&test1);

    return;
}

void isoex2(fp2 *j, const mpz_t sk2, const fp2 x1, const fp2 x2, const fp2 x3,
            const int e2) {
    MP A_C, P1, P2, P3, S, Ap24_C24, test1;
    MP_init(&A_C);
    MP_init(&P1);
    MP_init(&P2);
    MP_init(&P3);
    MP_init(&S);
    MP_init(&Ap24_C24);
    MP_init(&test1);

    MP_set_x(&P1, x1);
    MP_set_x(&P2, x2);
    MP_set_x(&P3, x3);

    fp2 A;
    fp2_init(&A);

    get_A(&A, P1, P2, P3);
    MP_set_x(&A_C, A);

    Ladder3pt(&S, sk2, P1, P2, P3, A_C);

    fp2_add_si(&(Ap24_C24.X), A, 2);
    fp2_set_si(&(Ap24_C24.Z), 0, 4);

    e_2_iso(&test1, Ap24_C24, S, e2);
    MP_copy(&Ap24_C24, test1);

    fp2_copy(&(A_C.Z), Ap24_C24.Z);
    fp2_add(&(global_params.fp2res), Ap24_C24.Z, Ap24_C24.Z);
    fp2_copy(&(Ap24_C24.Z), global_params.fp2res);
    fp2_add(&(global_params.fp2res), Ap24_C24.X, Ap24_C24.X);
    fp2_copy(&(Ap24_C24.X), global_params.fp2res);
    fp2_add(&(global_params.fp2res), Ap24_C24.X, Ap24_C24.X);
    fp2_copy(&(Ap24_C24.X), global_params.fp2res);
    fp2_sub(&(A_C.X), Ap24_C24.X, Ap24_C24.Z);

    jInvariant(j, A_C);

    MP_clear(&A_C);
    MP_clear(&P1);
    MP_clear(&P2);
    MP_clear(&P3);
    MP_clear(&S);
    MP_clear(&Ap24_C24);
    MP_clear(&test1);

    return;
}

void isoex3(fp2 *j, const mpz_t sk3, const fp2 x1, const fp2 x2, const fp2 x3,
            const int e3) {
    MP A_C, Ap24_Am24, P1, P2, P3, S, test1;
    MP_init(&A_C);
    MP_init(&Ap24_Am24);
    MP_init(&P1);
    MP_init(&P2);
    MP_init(&P3);
    MP_init(&S);
    MP_init(&test1);

    MP_set_x(&P1, x1);
    MP_set_x(&P2, x2);
    MP_set_x(&P3, x3);

    fp2 A;
    fp2_init(&A);

    get_A(&A, P1, P2, P3);
    MP_set_x(&A_C, A);

    Ladder3pt(&S, sk3, P1, P2, P3, A_C);

    fp2_add_si(&(Ap24_Am24.X), A, 2);
    fp2_sub_si(&(Ap24_Am24.Z), A, 2);

    e_3_iso(&test1, Ap24_Am24, S, e3);
    MP_copy(&Ap24_Am24, test1);

    fp2_sub(&(A_C.Z), Ap24_Am24.X, Ap24_Am24.Z);
    fp2_add(&(A_C.X), Ap24_Am24.Z, Ap24_Am24.X);
    fp2_add(&(global_params.fp2res), A_C.X, A_C.X);
    fp2_copy(&(A_C.X), global_params.fp2res);

    jInvariant(j, A_C);

    MP_clear(&A_C);
    MP_clear(&Ap24_Am24);
    MP_clear(&P1);
    MP_clear(&P2);
    MP_clear(&P3);
    MP_clear(&S);
    MP_clear(&test1);

    return;
}
