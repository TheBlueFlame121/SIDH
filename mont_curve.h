#include <gmp.h>
#include "fp2.h"

typedef struct {
    fp2 A;
}MC;

typedef struct {
    fp2 X, Z;
}MP;

void MC_init(MC *Res);

void MP_init(MP *Res);

void MP_set(MP *Res, const fp2 x, const fp2 z);

void MP_set_mpz(MP *Res, const mpz_t x1, const mpz_t x0, const mpz_t z1, const mpz_t z0);

void MP_set_si(MP *Res, const long int x1, const long int x0, const long int z1, const long int z0);

void MP_set_x(MP *Res, const fp2 x);

void MP_copy(MP *Res, const MP P);

void MP_clear(MP *Res);

void MP_norm(MP *Res, const MP P);

void xDBL(MP *Res, const MP P, const MP Ap24_C24); 

void xDBLe(MP *Res, const MP P, const MP Ap24_C24, int e);

// Be careful with the parmams
// In Public Key:= R = P-Q
// In function definition:= R = Q-P
// Might not make a difference as we are omitting y.
void xDBLADD(MP *Dub, MP *Add, const MP P, const MP Q, const MP R, const MP Ap24_1);

void xTPL(MP *Res, const MP P, const MP Ap24_Am24);

void xTPLe(MP *Res, const MP P, const MP Ap24i_Am24, int e);
    
// Be careful with the parmams
// In Public Key:= R = P-Q
// In function definition:= R = Q-P
// Might not make a difference as we are omitting y.
void Ladder3pt(MP *Res,const mpz_t m, const MP P, const MP Q, const MP R, const MP A_1);

void jInvariant(fp2 *j, const MP A_C);

void get_A(fp2 *A, const MP P, const MP Q, const MP R);

