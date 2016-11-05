#ifndef FFT_H
#define FFT_H

#include "fft_complex.h"

#include <stddef.h>

struct fft {
    fft_complex *cos_sin;
    size_t *reversed;
    size_t n;
};

void fft_init(struct fft *f, size_t n);
void fft_free(struct fft *f);
void fft_fft (struct fft *f, fft_complex *values);
void fft_ifft(struct fft *f, fft_complex *values);

#endif
