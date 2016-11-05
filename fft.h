#ifndef FFT_INCLUDED
#define FFT_INCLUDED

struct FFT {
    int n, *reversed;
    float *c, *s;
};

void fft_init(struct FFT *fft, int n);
void fft_free(struct FFT *fft);
void fft_fft(struct FFT *fft, float *re, float *im);
void fft_ifft(struct FFT *fft, float *re, float *im);

#endif
