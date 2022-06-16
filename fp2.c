#include "fp2.h"
#include <time.h>

fp2_params global_params;

void fp2_init(fp2 *Res) {
    mpz_init(Res->a);
    mpz_init(Res->b);
}

void fp2_setup(mpz_t prime) {
    // init prime and few mpz registers
    mpz_init(global_params.p);
    mpz_init(global_params.tmp1);
    mpz_init(global_params.tmp2);
    mpz_init(global_params.tmp3);

    mpz_set(global_params.p, prime);

    // init random
    gmp_randinit_mt(global_params.state);
    gmp_randseed_ui(global_params.state, time(0));

    // init some fp2 registers
    fp2_init(&(global_params.fp2tmp1));
    fp2_init(&(global_params.t0));
    fp2_init(&(global_params.t1));
    fp2_init(&(global_params.t2));
    fp2_init(&(global_params.t3));
    fp2_init(&(global_params.t4));
    fp2_init(&(global_params.t5));
    fp2_init(&(global_params.t6));
    fp2_init(&(global_params.fp2res));
}

void fp2_clear(fp2 *Res) {
    mpz_clear(Res->a);
    mpz_clear(Res->b);
}

void fp2_set_str(fp2 *Res, const char *a, const char *b) {
    mpz_set_str(global_params.tmp1, a, 0);
    mpz_mod(Res->a, global_params.tmp1, global_params.p);
    mpz_set_str(global_params.tmp1, b, 0);
    mpz_mod(Res->b, global_params.tmp1, global_params.p);
}

void fp2_set_str16(fp2 *Res, const char *a, const char *b) {
    mpz_set_str(global_params.tmp1, a, 16);
    mpz_mod(Res->a, global_params.tmp1, global_params.p);
    mpz_set_str(global_params.tmp1, b, 16);
    mpz_mod(Res->b, global_params.tmp1, global_params.p);
}

void fp2_set_mpz(fp2 *Res, const mpz_t a, const mpz_t b) {
    mpz_mod(Res->a, a, global_params.p);
    mpz_mod(Res->b, b, global_params.p);
}

void fp2_set_si(fp2 *Res, const long int a, const long int b) {
    mpz_set_si(global_params.tmp1, a);
    mpz_mod(Res->a, global_params.tmp1, global_params.p);
    mpz_set_si(global_params.tmp1, b);
    mpz_mod(Res->b, global_params.tmp1, global_params.p);
}

void fp2_to_string(const fp2 x, char *a, char *b) {
    gmp_sprintf(a, "%Zd", x.a);
    gmp_sprintf(b, "%Zd", x.b);
}

int fp2_iszero(const fp2 x) {
    return (mpz_cmp_si(x.a, 0) == 0) && (mpz_cmp_si(x.b, 0) == 0);
}

void fp2_print(const fp2 x, char *a) {
    if (fp2_iszero(x)) {
        gmp_sprintf(a, "0");
        return;
    } else if (mpz_cmp_si(x.a, 0) == 0) {
        gmp_sprintf(a, "%Zd", x.b);
        return;
    } else if (mpz_cmp_si(x.b, 0) == 0) {
        gmp_sprintf(a, "%Zdi", x.a);
        return;
    }
    gmp_sprintf(a, "%Zdi + %Zd", x.a, x.b);
}

void fp2_copy(fp2 *Res, const fp2 x) {
    mpz_set(Res->a, x.a);
    mpz_set(Res->b, x.b);
}

int fp2_cmp(const fp2 x, const fp2 y) {
    int flag = mpz_cmp(x.a, y.a);
    if (!flag)
        flag = mpz_cmp(x.b, y.b);
    return flag;
}

int fp2_isone(const fp2 x) {
    return (mpz_cmp_si(x.a, 0) == 0) && (mpz_cmp_si(x.b, 1) == 0);
}

void fp2_random(fp2 *Res) {
    mpz_urandomm(Res->a, global_params.state, global_params.p);
    mpz_urandomm(Res->b, global_params.state, global_params.p);
}

int fp2_isequal(const fp2 x, const fp2 y) { return (fp2_cmp(x, y) == 0); }

void fp2_add(fp2 *Res, const fp2 x, const fp2 y) {
    mpz_add(global_params.tmp1, x.a, y.a);
    mpz_mod(Res->a, global_params.tmp1, global_params.p);
    mpz_add(global_params.tmp1, x.b, y.b);
    mpz_mod(Res->b, global_params.tmp1, global_params.p);
}

void fp2_add_si(fp2 *Res, const fp2 x, const long int b) {
    mpz_add_ui(global_params.tmp1, x.b, b);
    mpz_mod(Res->b, global_params.tmp1, global_params.p);
    mpz_set(Res->a, x.a);
}

void fp2_sub(fp2 *Res, const fp2 x, const fp2 y) {
    mpz_sub(global_params.tmp1, x.a, y.a);
    mpz_mod(Res->a, global_params.tmp1, global_params.p);
    mpz_sub(global_params.tmp1, x.b, y.b);
    mpz_mod(Res->b, global_params.tmp1, global_params.p);
}

void fp2_sub_si(fp2 *Res, const fp2 x, const long int b) {
    mpz_sub_ui(global_params.tmp1, x.b, b);
    mpz_mod(Res->b, global_params.tmp1, global_params.p);
    mpz_set(Res->a, x.a);
}

void fp2_neg(fp2 *Res, const fp2 x) {
    if (mpz_cmp_si(x.a, 0) == 0)
        mpz_set(Res->a, x.a);
    else
        mpz_sub(Res->a, global_params.p, x.a);
    if (mpz_cmp_si(x.b, 0) == 0)
        mpz_set(Res->b, x.b);
    else
        mpz_sub(Res->b, global_params.p, x.b);
}

void fp2_scalar(fp2 *Res, const fp2 x, const mpz_t s) {
    mpz_mul(global_params.tmp1, x.a, s);
    mpz_mod(Res->a, global_params.tmp1, global_params.p);
    mpz_mul(global_params.tmp1, x.b, s);
    mpz_mod(Res->b, global_params.tmp1, global_params.p);
}

void fp2_scalar_si(fp2 *Res, const fp2 x, const long int s) {
    mpz_mul_si(global_params.tmp1, x.a, s);
    mpz_mod(Res->a, global_params.tmp1, global_params.p);
    mpz_mul_si(global_params.tmp1, x.b, s);
    mpz_mod(Res->b, global_params.tmp1, global_params.p);
}

void fp2_mul(fp2 *Res, const fp2 x, const fp2 y) {
    mpz_mul(global_params.tmp1, x.a, y.b);
    mpz_mul(global_params.tmp2, y.a, x.b);
    mpz_add(global_params.tmp3, global_params.tmp1, global_params.tmp2);
    mpz_mod(Res->a, global_params.tmp3, global_params.p);

    mpz_mul(global_params.tmp1, x.b, y.b);
    mpz_mul(global_params.tmp2, x.a, y.a);
    mpz_sub(global_params.tmp3, global_params.tmp1, global_params.tmp2);
    mpz_mod(Res->b, global_params.tmp3, global_params.p);
}

void fp2_sqr(fp2 *Res, const fp2 x) {
    mpz_mul(global_params.tmp1, x.a, x.b);
    mpz_add(global_params.tmp2, global_params.tmp1, global_params.tmp1);
    mpz_mod(Res->a, global_params.tmp2, global_params.p);

    mpz_mul(global_params.tmp1, x.a, x.a);
    mpz_mul(global_params.tmp2, x.b, x.b);
    mpz_sub(global_params.tmp3, global_params.tmp2, global_params.tmp1);
    mpz_mod(Res->b, global_params.tmp3, global_params.p);
}

//(a*i + b)^-1 = (b-i*a)/(b^2+a^2).
void fp2_inv(fp2 *Res, const fp2 x) {
    mpz_mul(global_params.tmp1, x.b, x.b);
    mpz_mul(global_params.tmp2, x.a, x.a);
    mpz_add(global_params.tmp3, global_params.tmp1, global_params.tmp2);
    mpz_sub_ui(global_params.tmp1, global_params.p, 2);
    mpz_powm(global_params.tmp2, global_params.tmp3, global_params.tmp1,
             global_params.p);

    mpz_mul(global_params.tmp1, x.b, global_params.tmp2);
    mpz_mod(Res->b, global_params.tmp1, global_params.p);

    mpz_mul(global_params.tmp1, x.a, global_params.tmp2);
    mpz_sub(global_params.tmp3, global_params.p, global_params.tmp1);
    mpz_mod(Res->a, global_params.tmp3, global_params.p);
}

void fp2_div(fp2 *Res, const fp2 x, const fp2 y) {
    fp2_inv(&(global_params.fp2tmp1), y);
    fp2_mul(Res, x, global_params.fp2tmp1);
}
