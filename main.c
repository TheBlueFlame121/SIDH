#include "sidh.h"
#include "utils.h"
#include <gmp.h>
#include <stdio.h>
#include <string.h>

#define IN_LEN 1000

int main(int argc, char *argv[]) {
  char sp[IN_LEN], seA[IN_LEN], seB[IN_LEN], sPA0[IN_LEN], sPA1[IN_LEN],
      sQA0[IN_LEN], sQA1[IN_LEN], sRA0[IN_LEN], sRA1[IN_LEN], sPB0[IN_LEN],
      sPB1[IN_LEN], sQB0[IN_LEN], sQB1[IN_LEN], sRB0[IN_LEN], sRB1[IN_LEN],
      file[100] = "params/";

  if (argc == 1) {
    printf("Usage: ./sidh_test <filename>\nOptional Second argument for "
           "verbose\n");
    printf("\nExample: ./sidh_test sikep16.txt v\n");
    exit(0);
  }
  strcat(file, argv[1]);
  int VERBOSE = 0;
  if (argc == 3)
    VERBOSE = 1;

  params_from_file(sp, seA, seB, sPA0, sPA1, sQA0, sQA1, sRA0, sRA1, sPB0, sPB1,
                   sQB0, sQB1, sRB0, sRB1, file);

  if (VERBOSE) {
    printf("----PARAMS READ----\n");
    printf("P = %s\n", sp);
    printf("eA = %s\n", seA);
    printf("eB = %s\n", seB);
    printf("PA0 = %s\n", sPA0);
    printf("PA1 = %s\n", sPA1);
    printf("QA0 = %s\n", sQA0);
    printf("QA1 = %s\n", sQA1);
    printf("RA0 = %s\n", sRA0);
    printf("RA1 = %s\n", sRA1);
    printf("PB0 = %s\n", sPB0);
    printf("PB1 = %s\n", sPB1);
    printf("QB0 = %s\n", sQB0);
    printf("QB1 = %s\n", sQB1);
    printf("RB0 = %s\n", sRB0);
    printf("RB1 = %s\n", sRB1);
    printf("-------------------\n\n");
  }

  mpz_t p, temp;
  mpz_init(p);
  mpz_set_str(p, sp, 0);
  fp2_setup(p);
  if (VERBOSE)
    gmp_printf("Field setup with prime: %#Zx\n", p);

  int eA, eB;
  mpz_init(temp);
  mpz_set_str(temp, seA, 0);
  eA = mpz_get_si(temp);
  mpz_set_str(temp, seB, 0);
  eB = mpz_get_si(temp);
  // sscanf(seA, "%d", &eA);
  // sscanf(seB, "%d", &eB);

  if (VERBOSE)
    printf("Public sizes: 2^%d and 3^%d\n", eA, eB);

  mpz_t temp1, temp2, sk2, sk3;
  mpz_init(temp1);
  mpz_init(temp2);
  mpz_init(sk2);
  mpz_init(sk3);

  mpz_set_si(temp1, 2);
  mpz_set_si(temp2, 3);
  mpz_pow_ui(temp1, temp1, eA);
  mpz_pow_ui(temp2, temp2, eB);
  mpz_urandomm(sk2, global_params.state, temp1);
  mpz_urandomm(sk3, global_params.state, temp2);

  // gmp_printf("\n\nSanity: %#Zx %#Zx\n", temp1, temp2);
  // gmp_printf("Private Keys: %#Zx %#Zx\n", sk2, sk3);

  // Alice Side parameter definition
  fp2 App, Aqp, Arp, Ap, Aq, Ar, sec_A;
  fp2_init(&App);
  fp2_init(&Aqp);
  fp2_init(&Arp);
  fp2_init(&Ap);
  fp2_init(&Aq);
  fp2_init(&Ar);
  fp2_init(&sec_A);

  fp2_set_str(&Ap, sPA1, sPA0);
  fp2_set_str(&Aq, sQA1, sQA0);
  fp2_set_str(&Ar, sRA1, sRA0);

  if (VERBOSE) {
    gmp_printf("\n\nAlice Gen points: xPA||xQA||xRA\n");
    gmp_printf("%#Zx*u + %#Zx\n", Ap.a, Ap.b);
    gmp_printf("%#Zx*u + %#Zx\n", Aq.a, Aq.b);
    gmp_printf("%#Zx*u + %#Zx\n", Ar.a, Ar.b);
  }

  // Bob Side parameter definition
  fp2 Bpp, Bqp, Brp, Bp, Bq, Br, sec_B;
  fp2_init(&Bpp);
  fp2_init(&Bqp);
  fp2_init(&Brp);
  fp2_init(&Bp);
  fp2_init(&Bq);
  fp2_init(&Br);
  fp2_init(&sec_B);

  fp2_set_str(&Bp, sPB1, sPB0);
  fp2_set_str(&Bq, sQB1, sQB0);
  fp2_set_str(&Br, sRB1, sRB0);

  if (VERBOSE) {
    gmp_printf("\n\nBob Gen points:   xPB||xQB||xRB\n");
    gmp_printf("%#Zx*u + %#Zx\n", Bp.a, Bp.b);
    gmp_printf("%#Zx*u + %#Zx\n", Bq.a, Bq.b);
    gmp_printf("%#Zx*u + %#Zx\n", Br.a, Br.b);
  }

  // Alice Round 1
  isogen2(&Bpp, &Bqp, &Brp, sk2, eA, eB, p, Ap, Aq, Ar, Bp, Bq, Br);
  if (VERBOSE) {
    gmp_printf("\nAlice Public Key: xPBp||xQBp||xRBp\n");
    gmp_printf("%#Zx*u + %#Zx\n", Bpp.a, Bpp.b);
    gmp_printf("%#Zx*u + %#Zx\n", Bqp.a, Bqp.b);
    gmp_printf("%#Zx*u + %#Zx\n", Brp.a, Brp.b);
  }

  // Bob Round 1
  isogen3(&App, &Aqp, &Arp, sk3, eA, eB, p, Ap, Aq, Ar, Bp, Bq, Br);
  if (VERBOSE) {
    gmp_printf("\nBob Public Key:   xPAp||xQAp||xRAp\n");
    gmp_printf("%#Zx*u + %#Zx\n", App.a, App.b);
    gmp_printf("%#Zx*u + %#Zx\n", Aqp.a, Aqp.b);
    gmp_printf("%#Zx*u + %#Zx\n", Arp.a, Arp.b);
  }

  // Alice Round 2
  isoex2(&sec_A, sk2, App, Aqp, Arp, eA);
  if (VERBOSE)
    gmp_printf("\n\nAlice shared secret:\n%#Zx*u + %#Zx\n", sec_A.a, sec_A.b);

  // Bob Round 2
  isoex3(&sec_B, sk3, Bpp, Bqp, Brp, eB);
  if (VERBOSE)
    gmp_printf("Bob shared secret:\n%#Zx*u + %#Zx\n", sec_B.a, sec_B.b);

  if (fp2_isequal(sec_A, sec_B))
    printf("\nSuccess\n");
  else
    printf("\nFailure\n");

  return 0;
}
