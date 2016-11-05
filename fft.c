#include "fft.h"

#include <math.h>
#include <stdlib.h>

static int fft_bitreverse(int x, int n){
    int i, result = 0;
    for (i = 0; i < n; i++){
        result = (result << 1) | ((x >> i) & 1);
    }
    return result;
}

static int fft_bitlength(int x){
    int n = 0;
    for (; x; x >>= 1) n++;
    return n;
}

void fft_init(struct FFT *fft, int n){
    fft->c = (float*)malloc(n*sizeof(*fft->c));
    fft->s = (float*)malloc(n*sizeof(*fft->s));
    fft->reversed = (int*)malloc(n*sizeof(*fft->reversed));
    fft->n = n;

    const float pi = 3.14159265358979f;

    int i, log2_n = fft_bitlength(n) - 1;
    for (i = 0; i < n; i++){
        float t = -2.0f*pi/n*i;
        fft->c[i] = cos(t);
        fft->s[i] = sin(t);
        fft->reversed[i] = fft_bitreverse(i, log2_n);
    }
}

void fft_free(struct FFT *fft){
    free(fft->reversed);
    free(fft->c);
    free(fft->s);
}

static void fft_butterfly_swap(struct FFT *fft, float *re, float *im){
    int i;
    for (i = 0; i < fft->n; i++){
        int j = fft->reversed[i];
        if (i < j){
            float temp;

            temp = re[i];
            re[i] = re[j];
            re[j] = temp;

            temp = im[i];
            im[i] = im[j];
            im[j] = temp;
        }
    }
}

void fft_fft(struct FFT *fft, float *re, float *im){
    int n = fft->n;
    float *c = fft->c;
    float *s = fft->s;
    fft_butterfly_swap(fft, re, im);
    int block_size, block_offset, i;
    for (block_size = 2; block_size <= n; block_size *= 2){
        int cos_sin_stride = n/block_size;
        for (block_offset = 0; block_offset < n; block_offset += block_size){
            for (i = 0; i < block_size/2; i++){
                float ar = re[block_offset + i];
                float ai = im[block_offset + i];
                float br = re[block_offset + i + block_size/2];
                float bi = im[block_offset + i + block_size/2];
                float cr = c[i*cos_sin_stride];
                float ci = s[i*cos_sin_stride];
                float dr = br*cr - bi*ci;
                float di = br*ci + bi*cr;
                re[block_offset + i               ] = ar + dr;
                im[block_offset + i               ] = ai + di;
                re[block_offset + i + block_size/2] = ar - dr;
                im[block_offset + i + block_size/2] = ai - di;
            }
        }
    }
}

void fft_ifft(struct FFT *fft, float *re, float *im){
    int n = fft->n;
    float *c = fft->c;
    float *s = fft->s;
    fft_butterfly_swap(fft, re, im);
    int block_size, block_offset, i;
    for (block_size = 2; block_size <= n; block_size *= 2){
        int cos_sin_stride = n/block_size;
        for (block_offset = 0; block_offset < n; block_offset += block_size){
            for (i = 0; i < block_size/2; i++){
                float ar = re[block_offset + i];
                float ai = im[block_offset + i];
                float br = re[block_offset + i + block_size/2];
                float bi = im[block_offset + i + block_size/2];
                float cr = c[i*cos_sin_stride];
                float ci = -s[i*cos_sin_stride];
                float dr = br*cr - bi*ci;
                float di = br*ci + bi*cr;
                re[block_offset + i               ] = ar + dr;
                im[block_offset + i               ] = ai + di;
                re[block_offset + i + block_size/2] = ar - dr;
                im[block_offset + i + block_size/2] = ai - di;
            }
        }
    }
}
