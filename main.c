#include "fft.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

float randf(){
    return rand()/(float)RAND_MAX;
}

int main(){
    size_t i, n = 512;

    struct FFT fft[1];
    fft_init(fft, n*2);

    float re[n*2];
    float im[n*2];
    float re_shifted[n*2];
    float im_shifted[n*2];

    for (i = 0; i < n; i++){
        re[i] = randf();
        im[i] = 0.0f;

        re[i + n] = 0.0f;
        im[i + n] = 0.0f;
    }

    size_t shift = 13;

    for (i = 0; i < n; i++){
        size_t j = (i + shift) % n;
        re_shifted[i] = re[j];
        im_shifted[i] = im[j];

        re_shifted[i + n] = 0.0f;
        im_shifted[i + n] = 0.0f;
    }

    fft_fft(fft, re, im);
    fft_fft(fft, re_shifted, im_shifted);

    for (i = 0; i < n*2; i++){
        re[i] *= re_shifted[i];
        im[i] *= im_shifted[i];
    }

    fft_ifft(fft, re, im);

    size_t index_max = -1;
    float max_value = -1.0f;
    for (i = 1; i < n*2; i++){
        float abs_value = re[i];

        if (abs_value < 0.0f){
            abs_value = -abs_value;
        }

        if (abs_value > max_value){
            max_value = abs_value;
            index_max = i;
        }
    }

    printf("index: %i\n", (int)index_max);

    fft_free(fft);

    return 0;
}
