#include <gmp.h>
#include "mont_curve.h"

void iso_2_curve(MP *Ap24_C24, const MP P);

void iso_2_eval(MP *Res, const MP P2, const MP Q);

void iso_4_curve(MP *Ap24_C24, fp2 *K1, fp2 *K2, fp2 *K3, const MP P4);

void iso_4_eval(MP *Res, const fp2 K1, const fp2 K2, const fp2 K3, const MP Q);

void iso_3_curve(MP *Ap24_Am24, fp2 *K1, fp2 *K2, const MP P3);

void iso_e_3val(MP *Res, const fp2 K1, const fp2 K2, const MP Q);

void e_2_iso_optional(MP *ACp, MP *Res1, MP *Res2, MP *Res3, const MP AC, const MP S, const MP P1, const MP P2, const MP P3, int e2); 

void e_2_iso(MP *ACp, const MP AC, const MP S, int e2); 

void e_3_iso_optional(MP *AAp, MP *Res1, MP *Res2, MP *Res3, const MP AA, const MP S, const MP P1, const MP P2, const MP P3, int e3); 

void e_3_iso(MP *AAp, const MP AA, const MP S, int e3); 
