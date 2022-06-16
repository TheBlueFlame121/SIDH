#include <gmp.h>

typedef struct{
	mpz_t a, b;
}fp2;

struct fp2_params;
typedef struct fp2_params fp2_params;

struct fp2_params{
	mpz_t p, tmp1, tmp2, tmp3;
	gmp_randstate_t state;
	fp2 fp2tmp1, t0, t1, t2, t3, t4, t5, t6, fp2res;
};

extern fp2_params global_params;

void fp2_init(fp2 *Res);

void fp2_setup(mpz_t prime);

void fp2_clear(fp2 *Res);

void fp2_set_str(fp2 *Res, const char *a, const char *b);

void fp2_set_str16(fp2 *Res, const char *a, const char *b);

void fp2_set_mpz(fp2 *Res, const mpz_t a, const mpz_t b);

void fp2_set_si(fp2 *Res, const long int a, const long int b);

void fp2_to_string(const fp2 x, char *a, char *b);

int fp2_iszero(const fp2 x);

void fp2_print(const fp2 x, char *a);

void fp2_copy(fp2 *Res, const fp2 x);

int fp2_cmp(const fp2 x, const fp2 y);

int fp2_isone(const fp2 x);

void fp2_random(fp2 *Res);

int fp2_isequal(const fp2 x, const fp2 y);

void fp2_add(fp2 *Res, const fp2 x, const fp2 y);

void fp2_add_si(fp2 *Res, const fp2 x, const long int b);

void fp2_sub(fp2 *Res, const fp2 x, const fp2 y);

void fp2_sub_si(fp2 *Res, const fp2 x, const long int b);

void fp2_neg(fp2 *Res, const fp2 x);

void fp2_scalar(fp2 *Res, const fp2 x, const mpz_t s);

void fp2_scalar_si(fp2 *Res, const fp2 x, const long int s);

void fp2_mul(fp2 *Res, const fp2 x, const fp2 y);

void fp2_sqr(fp2 *Res, const fp2 x);

//(a*i + b)^-1 = (b-i*a)/(b^2+a^2).
void fp2_inv(fp2 *Res, const fp2 x);

void fp2_div(fp2 *Res, const fp2 x, const fp2 y);
