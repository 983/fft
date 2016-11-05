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

static void fft_swap_and_first_pass(struct FFT *fft, float *re, float *im){
    int i, n = fft->n;
    for (i = 0; i < n; i++){
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

    if (n >= 2){
        for (i = 0; i < n; i += 2){
            float ar = re[i + 0];
            float ai = im[i + 0];
            float br = re[i + 1];
            float bi = im[i + 1];
            re[i + 0] = ar + br;
            im[i + 0] = ai + bi;
            re[i + 1] = ar - br;
            im[i + 1] = ai - bi;
        }
    }
}

void fft_fft(struct FFT *fft, float *re, float *im){
    int block_size, block_offset, i, n = fft->n;
    float *c = fft->c;
    float *s = fft->s;

    fft_swap_and_first_pass(fft, re, im);

    if (n >= 4){
        for (i = 0; i < n; i += 4){
            float ar = re[i + 0];
            float ai = im[i + 0];
            float cr = re[i + 1];
            float ci = im[i + 1];
            float br = re[i + 2];
            float bi = im[i + 2];
            float di = re[i + 3];
            float dr = im[i + 3];
            re[i + 0] = ar + br;
            im[i + 0] = ai + bi;
            re[i + 1] = cr - dr;
            im[i + 1] = ci - di;
            re[i + 2] = ar - br;
            im[i + 2] = ai - bi;
            re[i + 3] = cr + dr;
            im[i + 3] = ci + di;
        }
    }

    for (block_size = 8; block_size <= n; block_size *= 2){
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
    int block_size, block_offset, i, n = fft->n;
    float *c = fft->c;
    float *s = fft->s;

    fft_swap_and_first_pass(fft, re, im);

    if (n >= 4){
        for (i = 0; i < n; i += 4){
            float ar = re[i + 0];
            float ai = im[i + 0];
            float cr = re[i + 1];
            float ci = im[i + 1];
            float br = re[i + 2];
            float bi = im[i + 2];
            float di = -re[i + 3];
            float dr = im[i + 3];
            re[i + 0] = ar + br;
            im[i + 0] = ai + bi;
            re[i + 1] = cr - dr;
            im[i + 1] = ci - di;
            re[i + 2] = ar - br;
            im[i + 2] = ai - bi;
            re[i + 3] = cr + dr;
            im[i + 3] = ci + di;
        }
    }

    for (block_size = 8; block_size <= n; block_size *= 2){
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

    float scale = 1.0f/n;
    for (i = 0; i < n; i++){
        re[i] *= scale;
        im[i] *= scale;
    }
}
