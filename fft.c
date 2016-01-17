#include "fft.h"

#include <assert.h>
#include <stdlib.h>

/* get index of highest set bit in x */
static size_t fft_bitlength(size_t x){
    size_t n = 0;
    for (; x; x >>= 1) n++;
    return n;
}

static size_t fft_bitreverse(size_t x, size_t n){
    size_t i, result = 0;
    for (i = 0; i < n; i++){
        result = (result << 1) | (x & 1);
        x >>= 1;
    }
    return result;
}

void fft_init(struct fft *f, size_t n){
    size_t i, log2_n = fft_bitlength(n) - 1;

    f->n = n;
    f->cos_sin = (fft_complex*)malloc(sizeof(*f->cos_sin)*n/2);
    f->reversed = (size_t*)malloc(sizeof(*f->reversed)*n);

    assert("n must be exact power of two" && n == ((size_t)1) << log2_n);

    for (i = 0; i < n/2; i++){
        f->cos_sin[i] = fftc_polar(-2.0*FFT_PI/n*i);
    }

    for (i = 0; i < n; i++){
        f->reversed[i] = fft_bitreverse(i, log2_n);
    }
}

void fft_free(struct fft *f){
    free(f->cos_sin);
    free(f->reversed);
}

static void fft_butterfly_swap(struct fft *f, fft_complex *values){
    size_t i, j, n = f->n;

    for (i = 0; i < n; i++){
        j = f->reversed[i];
        if (i < j){
            fft_complex temp = values[i];
            values[i] = values[j];
            values[j] = temp;
        }
    }
}

void fft_fft(struct fft *f, fft_complex *values){
    size_t i, block_size, block_offset, n = f->n;

    /* fourier transform of a single element is the same element */
    if (n < 2) return;

    fft_butterfly_swap(f, values);

    for (i = 0; i < n; i += 2){
        fft_complex a = values[i + 0];
        fft_complex b = values[i + 1];
        values[i + 0] = fftc_add(a, b);
        values[i + 1] = fftc_sub(a, b);
    }

    if (n < 4) return;

    for (i = 0; i < n; i += 4){
        fft_complex a = values[i + 0];
        fft_complex c = values[i + 1];
        fft_complex b = values[i + 2];
        fft_complex d = values[i + 3];
        d = fftc(d.imag, d.real);
        values[i + 0] = fftc_add(a, b);
        values[i + 1] = fftc_sub(c, d);
        values[i + 2] = fftc_sub(a, b);
        values[i + 3] = fftc_add(c, d);
    }

    for (block_size = 8; block_size <= n; block_size *= 2){
        size_t cos_sin_stride = n/block_size;
        for (block_offset = 0; block_offset < n; block_offset += block_size){
            for (i = 0; i < block_size/2; i++){
                fft_complex c = f->cos_sin[i*cos_sin_stride];
                fft_complex a = values[block_offset + i];
                fft_complex b = values[block_offset + i + block_size/2];
                fft_complex d = fftc_mul(b, c);
                values[block_offset + i               ] = fftc_add(a, d);
                values[block_offset + i + block_size/2] = fftc_sub(a, d);
            }
        }
    }
}

/* there are only two differences to fft_fft */
/* first, cos_sin is complex conjugated      */
/* second, the result is scaled by `1.0/n`   */
void fft_ifft(struct fft *f, fft_complex *values){
    size_t i, block_size, block_offset, n = f->n;

    if (n < 2) return;

    fft_butterfly_swap(f, values);

    for (i = 0; i < n; i++) values[i] = fftc_muls(values[i], 1.0/n);

    for (i = 0; i < n; i += 2){
        fft_complex a = values[i + 0];
        fft_complex b = values[i + 1];
        values[i + 0] = fftc_add(a, b);
        values[i + 1] = fftc_sub(a, b);
    }

    if (n < 4) return;

    for (i = 0; i < n; i += 4){
        fft_complex a = values[i + 0];
        fft_complex c = values[i + 1];
        fft_complex b = values[i + 2];
        fft_complex d = values[i + 3];
        d = fftc(d.imag, -d.real);
        values[i + 0] = fftc_add(a, b);
        values[i + 1] = fftc_sub(c, d);
        values[i + 2] = fftc_sub(a, b);
        values[i + 3] = fftc_add(c, d);
    }

    for (block_size = 8; block_size <= n; block_size *= 2){
        size_t cos_sin_stride = n/block_size;
        for (block_offset = 0; block_offset < n; block_offset += block_size){
            for (i = 0; i < block_size/2; i++){
                fft_complex c = fftc_conj(f->cos_sin[i*cos_sin_stride]);
                fft_complex a = values[block_offset + i];
                fft_complex b = values[block_offset + i + block_size/2];
                fft_complex d = fftc_mul(b, c);
                values[block_offset + i               ] = fftc_add(a, d);
                values[block_offset + i + block_size/2] = fftc_sub(a, d);
            }
        }
    }
}
