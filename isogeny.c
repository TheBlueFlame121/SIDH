#include "isogeny.h"
#include <gmp.h>

void iso_2_curve(MP *Ap24_C24, const MP P) {
    fp2_sqr(&(Ap24_C24->X), P.X);
    fp2_sqr(&(Ap24_C24->Z), P.Z);
    fp2_sub(&(global_params.fp2res), Ap24_C24->Z, Ap24_C24->X);
    fp2_copy(&(Ap24_C24->X), global_params.fp2res);
    return;
}

void iso_2_eval(MP *Res, const MP P2, const MP Q) {
    fp2_add(&(global_params.t0), P2.X, P2.Z);
    fp2_sub(&(global_params.t1), P2.X, P2.Z);
    fp2_add(&(global_params.t2), Q.X, Q.Z);

    fp2_sub(&(global_params.t3), Q.X, Q.Z);
    fp2_mul(&(global_params.fp2res), global_params.t0, global_params.t3);
    fp2_copy(&(global_params.t0), global_params.fp2res);
    fp2_mul(&(global_params.fp2res), global_params.t1, global_params.t2);
    fp2_copy(&(global_params.t1), global_params.fp2res);

    fp2_add(&(global_params.t2), global_params.t0, global_params.t1);
    fp2_sub(&(global_params.t3), global_params.t0, global_params.t1);
    fp2_mul(&(Res->X), Q.X, global_params.t2);
    fp2_mul(&(Res->Z), Q.Z, global_params.t3);
    return;
}

void iso_4_curve(MP *Ap24_C24, fp2 *K1, fp2 *K2, fp2 *K3, const MP P4) {
    fp2_sub(K2, P4.X, P4.Z);
    fp2_add(K3, P4.X, P4.Z);
    fp2_sqr(K1, P4.Z);
    fp2_add(&(global_params.fp2res), *K1, *K1);
    fp2_copy(K1, global_params.fp2res);
    fp2_sqr(&(Ap24_C24->Z), *K1);
    fp2_add(&(global_params.fp2res), *K1, *K1);
    fp2_copy(K1, global_params.fp2res);
    fp2_sqr(&(Ap24_C24->X), P4.X);
    fp2_add(&(global_params.fp2res), Ap24_C24->X, Ap24_C24->X);
    fp2_copy(&(Ap24_C24->X), global_params.fp2res);
    fp2_sqr(&(global_params.fp2res), Ap24_C24->X);
    fp2_copy(&(Ap24_C24->X), global_params.fp2res);
    return;
}

void iso_4_eval(MP *Res, const fp2 K1, const fp2 K2, const fp2 K3, const MP Q) {
    fp2_add(&(global_params.t0), Q.X, Q.Z);
    fp2_sub(&(global_params.t1), Q.X, Q.Z);
    fp2_mul(&(global_params.t2), global_params.t0, K2);
    fp2_mul(&(global_params.t3), global_params.t1, K3);
    fp2_mul(&(global_params.fp2res), global_params.t0, global_params.t1);
    fp2_copy(&(global_params.t0), global_params.fp2res);
    fp2_mul(&(global_params.fp2res), global_params.t0, K1);
    fp2_copy(&(global_params.t0), global_params.fp2res);
    fp2_add(&(global_params.t1), global_params.t2, global_params.t3);
    fp2_sub(&(global_params.fp2res), global_params.t2, global_params.t3);
    fp2_copy(&(global_params.t3), global_params.fp2res);
    fp2_sqr(&(global_params.fp2res), global_params.t1);
    fp2_copy(&(global_params.t1), global_params.fp2res);
    fp2_sqr(&(global_params.fp2res), global_params.t3);
    fp2_copy(&(global_params.t3), global_params.fp2res);
    fp2_add(&(global_params.t2), global_params.t0, global_params.t1);
    fp2_sub(&(global_params.fp2res), global_params.t3, global_params.t0);
    fp2_copy(&(global_params.t0), global_params.fp2res);
    fp2_mul(&(Res->X), global_params.t2, global_params.t1);
    fp2_mul(&(Res->Z), global_params.t3, global_params.t0);
    return;
}

void iso_3_curve(MP *Ap24_Am24, fp2 *K1, fp2 *K2, const MP P3) {
    fp2_sub(K1, P3.X, P3.Z);
    fp2_sqr(&(global_params.t0), *K1);
    fp2_add(K2, P3.X, P3.Z);
    fp2_sqr(&(global_params.t1), *K2);
    fp2_add(&(global_params.t2), global_params.t0, global_params.t1);

    fp2_add(&(global_params.t3), *K1, *K2);
    fp2_sqr(&(global_params.fp2res), global_params.t3);
    fp2_copy(&(global_params.t3), global_params.fp2res);
    fp2_sub(&(global_params.fp2res), global_params.t3, global_params.t2);
    fp2_copy(&(global_params.t3), global_params.fp2res);
    fp2_add(&(global_params.t2), global_params.t1, global_params.t3);
    fp2_add(&(global_params.fp2res), global_params.t3, global_params.t0);
    fp2_copy(&(global_params.t3), global_params.fp2res);

    fp2_add(&(global_params.t4), global_params.t3, global_params.t0);
    fp2_add(&(global_params.fp2res), global_params.t4, global_params.t4);
    fp2_copy(&(global_params.t4), global_params.fp2res);
    fp2_add(&(global_params.fp2res), global_params.t1, global_params.t4);
    fp2_copy(&(global_params.t4), global_params.fp2res);
    fp2_mul(&(Ap24_Am24->Z), global_params.t2, global_params.t4);
    fp2_add(&(global_params.t4), global_params.t1, global_params.t2);

    fp2_add(&(global_params.fp2res), global_params.t4, global_params.t4);
    fp2_copy(&(global_params.t4), global_params.fp2res);
    fp2_add(&(global_params.fp2res), global_params.t0, global_params.t4);
    fp2_copy(&(global_params.t4), global_params.fp2res);
    fp2_mul(&(Ap24_Am24->X), global_params.t3, global_params.t4);
    return;
}

void iso_e_3val(MP *Res, const fp2 K1, const fp2 K2, const MP Q) {
    fp2_add(&(global_params.t0), Q.X, Q.Z);
    fp2_sub(&(global_params.t1), Q.X, Q.Z);
    fp2_mul(&(global_params.fp2res), K1, global_params.t0);
    fp2_copy(&(global_params.t0), global_params.fp2res);

    fp2_mul(&(global_params.fp2res), K2, global_params.t1);
    fp2_copy(&(global_params.t1), global_params.fp2res);
    fp2_add(&(global_params.t2), global_params.t0, global_params.t1);
    fp2_sub(&(global_params.fp2res), global_params.t1, global_params.t0);
    fp2_copy(&(global_params.t0), global_params.fp2res);

    fp2_sqr(&(global_params.fp2res), global_params.t2);
    fp2_copy(&(global_params.t2), global_params.fp2res);
    fp2_sqr(&(global_params.fp2res), global_params.t0);
    fp2_copy(&(global_params.t0), global_params.fp2res);
    fp2_mul(&(Res->X), Q.X, global_params.t2);

    fp2_mul(&(Res->Z), Q.Z, global_params.t0);
    return;
}

void e_2_iso_optional(MP *ACp, MP *Res1, MP *Res2, MP *Res3, const MP AC,
                      const MP S, const MP P1, const MP P2, const MP P3,
                      int e2) {
    MP T, Sp, MP_temp;
    MP_init(&T);
    MP_init(&Sp);
    MP_init(&MP_temp);
    MP_copy(&Sp, S);
    MP_copy(ACp, AC);
    MP_copy(Res1, P1);
    MP_copy(Res2, P2);
    MP_copy(Res3, P3);

    fp2 K1, K2, K3;
    fp2_init(&K1);
    fp2_init(&K2);
    fp2_init(&K3);

    for (int e = e2 - 2; e >= 0; e -= 2) {
        xDBLe(&T, Sp, *ACp, e);
        iso_4_curve(ACp, &K1, &K2, &K3, T);
        if (e != 0) {
            iso_4_eval(&MP_temp, K1, K2, K3, Sp);
            MP_copy(&Sp, MP_temp);
        }
        iso_4_eval(&MP_temp, K1, K2, K3, *Res1);
        MP_copy(Res1, MP_temp);
        iso_4_eval(&MP_temp, K1, K2, K3, *Res2);
        MP_copy(Res2, MP_temp);
        iso_4_eval(&MP_temp, K1, K2, K3, *Res3);
        MP_copy(Res3, MP_temp);
    }

    MP_clear(&T);
    MP_clear(&Sp);
    MP_clear(&MP_temp);

    fp2_clear(&K1);
    fp2_clear(&K2);
    fp2_clear(&K3);

    return;
}

void e_2_iso(MP *ACp, const MP AC, const MP S, int e2) {
    MP T, Sp, MP_temp;
    MP_init(&T);
    MP_init(&Sp);
    MP_init(&MP_temp);
    MP_copy(&Sp, S);
    MP_copy(ACp, AC);

    fp2 K1, K2, K3;
    fp2_init(&K1);
    fp2_init(&K2);
    fp2_init(&K3);

    for (int e = e2 - 2; e >= 0; e -= 2) {
        xDBLe(&T, Sp, *ACp, e);
        iso_4_curve(ACp, &K1, &K2, &K3, T);
        if (e != 0) {
            iso_4_eval(&MP_temp, K1, K2, K3, Sp);
            MP_copy(&Sp, MP_temp);
        }
    }

    MP_clear(&T);
    MP_clear(&Sp);
    MP_clear(&MP_temp);

    fp2_clear(&K1);
    fp2_clear(&K2);
    fp2_clear(&K3);

    return;
}

void e_3_iso_optional(MP *AAp, MP *Res1, MP *Res2, MP *Res3, const MP AA, const MP S, const MP P1, const MP P2, const MP P3, int e3) {
    MP T, Sp, MP_temp;
    MP_init(&T);
    MP_init(&Sp);
    MP_init(&MP_temp);
    MP_copy(&Sp, S);
    MP_copy(AAp, AA);
    MP_copy(Res1, P1);
    MP_copy(Res2, P2);
    MP_copy(Res3, P3);

    fp2 K1, K2;
    fp2_init(&K1);
    fp2_init(&K2);

    for (int e = e3 - 1; e >= 0; e -= 1) {
        xTPLe(&T, Sp, *AAp, e);
        iso_3_curve(AAp, &K1, &K2, T);
        if (e != 0) {
            iso_e_3val(&MP_temp, K1, K2, Sp);
            MP_copy(&Sp, MP_temp);
        }
        iso_e_3val(&MP_temp, K1, K2, *Res1);
        MP_copy(Res1, MP_temp);
        iso_e_3val(&MP_temp, K1, K2, *Res2);
        MP_copy(Res2, MP_temp);
        iso_e_3val(&MP_temp, K1, K2, *Res3);
        MP_copy(Res3, MP_temp);
    }

    MP_clear(&T);
    MP_clear(&Sp);
    MP_clear(&MP_temp);

    fp2_clear(&K1);
    fp2_clear(&K2);

    return;
}

void e_3_iso(MP *AAp, const MP AA, const MP S, int e3) {
    MP T, Sp, MP_temp;
    MP_init(&T);
    MP_init(&Sp);
    MP_init(&MP_temp);
    MP_copy(&Sp, S);
    MP_copy(AAp, AA);

    fp2 K1, K2;
    fp2_init(&K1);
    fp2_init(&K2);

    for (int e = e3 - 1; e >= 0; e -= 1) {
        xTPLe(&T, Sp, *AAp, e);
        iso_3_curve(AAp, &K1, &K2, T);
        if (e != 0) {
            iso_e_3val(&MP_temp, K1, K2, Sp);
            MP_copy(&Sp, MP_temp);
        }
    }

    MP_clear(&T);
    MP_clear(&Sp);
    MP_clear(&MP_temp);

    fp2_clear(&K1);
    fp2_clear(&K2);

    return;
}
