#include "fft.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

float my_abs(float x){
    return x < 0.0f ? -x : x;
}

float randf(){
    return rand()/(float)RAND_MAX;
}

int calculate_shift(float *re, float *im, float *re_shifted, float *im_shifted, int n){
    struct FFT fft[1];
    fft_init(fft, n);

    fft_fft(fft, re, im);
    fft_fft(fft, re_shifted, im_shifted);

    int i;
    for (i = 0; i < n; i++){
        re[i] *= re_shifted[i];
        im[i] *= im_shifted[i];
    }

    fft_ifft(fft, re, im);

    int index_max = -1;
    float max_value = -1.0f;
    for (i = 0; i < n/2; i++){
        float abs_value = my_abs(re[i]);

        if (abs_value > max_value){
            max_value = abs_value;
            index_max = i;
        }
    }

    fft_free(fft);

    return index_max;
}

void test_shift(int shift){
    int i, n = 512;

    float *re         = (float*)malloc(sizeof(*re)*n*2);
    float *im         = (float*)malloc(sizeof(*im)*n*2);
    float *re_shifted = (float*)malloc(sizeof(*re)*n*2);
    float *im_shifted = (float*)malloc(sizeof(*im)*n*2);

    /* initialize an array with random values and padded zeros */
    for (i = 0; i < n; i++){
        re[i] = randf();
        im[i] = 0.0f;

        re[i + n] = 0.0f;
        im[i + n] = 0.0f;
    }

    /* create a shifted signal */
    for (i = 0; i < n; i++){
        int j = (i + shift) % n;
        re_shifted[i] = re[j];
        im_shifted[i] = im[j];

        re_shifted[i + n] = 0.0f;
        im_shifted[i + n] = 0.0f;
    }

    int calculated_shift = calculate_shift(re, im, re_shifted, im_shifted, n*2);

    assert(shift == calculated_shift);

    free(re);
    free(im);
    free(re_shifted);
    free(im_shifted);
}

void test_frequency(int frequency, int n){
    const float pi = 3.14159265358979f;

    float *re = (float*)malloc(sizeof(*re)*n);
    float *im = (float*)malloc(sizeof(*im)*n);

    int i;
    for (i = 0; i < n; i++){
        re[i] = cos(2.0f*pi*frequency/n*i);
        im[i] = 0.0f;
    }

    struct FFT fft[1];
    fft_init(fft, n);
    fft_fft(fft, re, im);

    /* The bin re[frequency] and re[n - frequency] should be n/2. */
    /* All other bins should be 0. */
    float max_error = 0.01f;

    for (i = 0; i < n; i++){
        assert(my_abs(im[i]) < max_error);
    }

    for (i = 0; i < frequency; i++){
        assert(my_abs(re[i]) < max_error);
    }

    assert(my_abs(re[frequency] - n/2) < max_error);

    for (i = frequency + 1; i < n - frequency; i++){
        assert(my_abs(re[i]) < max_error);
    }

    assert(my_abs(re[n - frequency] - n/2) < max_error);

    for (i = n - frequency + 1; i < n; i++){
        assert(my_abs(re[i]) < max_error);
    }

    fft_ifft(fft, re, im);

    /* check that x = fft(ifft(x)) */

    for (i = 0; i < n; i++){
        float x = cos(2.0f*pi*frequency/n*i);
        assert(my_abs(re[i] - x) < max_error);
        assert(my_abs(im[i]) < max_error);
    }

    fft_free(fft);
    free(re);
    free(im);
}

int main(){

    for (int shift = 0; shift < 100; shift++){
        test_shift(shift);
    }

    for (int n = 16; n <= 256; n *= 2){
        for (int frequency = 1; frequency < n/2; frequency++){
            test_frequency(frequency, n);
        }
    }

    return 0;
}
