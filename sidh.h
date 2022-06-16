#include <gmp.h>
#include "isogeny.h"

void isogen2(fp2 *x1, fp2 *x2, fp2 *x3, mpz_t sk2, int e2, int e3, mpz_t p,
             fp2 xP2, fp2 xQ2, fp2 xR2, fp2 xP3, fp2 xQ3, fp2 xR3);

void isogen3(fp2 *x1, fp2 *x2, fp2 *x3, mpz_t sk3, int e2, int e3, mpz_t p,
             fp2 xP2, fp2 xQ2, fp2 xR2, fp2 xP3, fp2 xQ3, fp2 xR3);
    
void isoex2(fp2 *j, const mpz_t sk2, const fp2 x1, const fp2 x2, const fp2 x3,
            const int e2);
   
void isoex3(fp2 *j, const mpz_t sk3, const fp2 x1, const fp2 x2, const fp2 x3,
            const int e3);
  
