#ifndef FFT_COMPLEX
#define FFT_COMPLEX

#define FFT_PI 3.1415926535897932384626433832795

typedef double fft_float;

typedef struct {
    fft_float real;
    fft_float imag;
} fft_complex;

fft_complex fftc(fft_float real, fft_float imag);
fft_complex fftc_add(fft_complex a, fft_complex b);
fft_complex fftc_sub(fft_complex a, fft_complex b);
fft_complex fftc_mul(fft_complex a, fft_complex b);
fft_complex fftc_muls(fft_complex a, fft_float b);
fft_complex fftc_conj(fft_complex a);
fft_complex fftc_polar(fft_float radians);
fft_float   fftc_abs(fft_complex a);
fft_float   fftc_angle(fft_complex a);

#endif
