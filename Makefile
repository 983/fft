all:
	gcc fft_complex.c fft.c test.c -o test -Wall -Wextra -ansi -pedantic -O2
