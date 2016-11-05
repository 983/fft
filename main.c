#include "fft.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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
    for (i = 1; i < n; i++){
        float abs_value = re[i];

        if (abs_value < 0.0f){
            abs_value = -abs_value;
        }

        if (abs_value > max_value){
            max_value = abs_value;
            index_max = i;
        }
    }

    fft_free(fft);

    return index_max;
}

float randf(){
    return rand()/(float)RAND_MAX;
}

int main(){
    int i, n = 512;

    float re[n*2];
    float im[n*2];
    float re_shifted[n*2];
    float im_shifted[n*2];

    /* initialize an array with random values and padded zeros */
    for (i = 0; i < n; i++){
        re[i] = randf();
        im[i] = 0.0f;

        re[i + n] = 0.0f;
        im[i + n] = 0.0f;
    }

    int shift = 13;

    /* create a shifted signal */
    for (i = 0; i < n; i++){
        int j = (i + shift) % n;
        re_shifted[i] = re[j];
        im_shifted[i] = im[j];

        re_shifted[i + n] = 0.0f;
        im_shifted[i + n] = 0.0f;
    }

    printf("signal is shifted by %i\n", calculate_shift(re, im, re_shifted, im_shifted, n*2));

    return 0;
}
