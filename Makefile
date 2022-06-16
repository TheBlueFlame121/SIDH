sidh_test: main.c utils.c fp2.c mont_curve.c isogeny.c sidh.c
	gcc -o sidh_test main.c utils.c fp2.c mont_curve.c isogeny.c sidh.c -lgmp
