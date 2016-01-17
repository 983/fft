#include "fft.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#define MAX_ALLOWED_ERROR 1e-10

void test_fft2(){
    struct fft f[1];

    fft_complex a[2];
    fft_complex b[2];

    fft_init(f, 2);

    a[0] = b[0] = fftc(1, 3);
    a[1] = b[1] = fftc(2, 4);

    fft_fft(f, b);

    assert(fftc_abs(fftc_sub(b[0], fftc( 3,  7))) < MAX_ALLOWED_ERROR);
    assert(fftc_abs(fftc_sub(b[1], fftc(-1, -1))) < MAX_ALLOWED_ERROR);

    fft_ifft(f, b);

    assert(fftc_abs(fftc_sub(b[0], a[0])) < MAX_ALLOWED_ERROR);
    assert(fftc_abs(fftc_sub(b[1], a[1])) < MAX_ALLOWED_ERROR);

    fft_free(f);
}

void test_fft_at_least_4(size_t n, size_t frequency){
    size_t i;
    struct fft f[1];
    fft_complex *a = (fft_complex*)malloc(sizeof(*a)*n);
    fft_complex *b = (fft_complex*)malloc(sizeof(*a)*n);

    /* fill the arrays `a` and `b` with a cosine wave of `frequency` Hertz` */
    for (i = 0; i < n; i++){
        a[i] = b[i] = fftc(cos(i*2.0/n*FFT_PI*frequency), 0.0);
    }

    fft_init(f, n);

    fft_fft(f, b);

    /* the fourier transform of a cosine wave has two non-zero entries */
    /* those are at the index `frequency` and `n - frequency` */
    /* they have a magnitude of n/2 */

    assert(fftc_abs(fftc_sub(b[    frequency], fftc(n/2, 0.0))) < MAX_ALLOWED_ERROR);
    assert(fftc_abs(fftc_sub(b[n - frequency], fftc(n/2, 0.0))) < MAX_ALLOWED_ERROR);
    for (i =                 0; i <     frequency; i++) assert(fftc_abs(b[i]) < MAX_ALLOWED_ERROR);
    for (i =     frequency + 1; i < n - frequency; i++) assert(fftc_abs(b[i]) < MAX_ALLOWED_ERROR);
    for (i = n - frequency + 1; i < n            ; i++) assert(fftc_abs(b[i]) < MAX_ALLOWED_ERROR);

    fft_ifft(f, b);

    for (i = 0; i < n; i++) assert(fftc_abs(fftc_sub(a[i], b[i])) < MAX_ALLOWED_ERROR);

    fft_free(f);
    free(a);
    free(b);
}

int main(){
    size_t frequency, n;

    test_fft2();

    for (n = 4; n <= 512; n *= 2){
        for (frequency = 1; frequency < n/2; frequency++){
            test_fft_at_least_4(n, frequency);
        }
    }

    printf("all tests passed!\n");

    return 0;
}
