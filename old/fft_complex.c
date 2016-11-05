#include "fft_complex.h"

#include <math.h>

fft_complex fftc(fft_float real, fft_float imag){
    fft_complex c;
    c.real = real;
    c.imag = imag;
    return c;
}

fft_complex fftc_add(fft_complex a, fft_complex b){
    fft_complex c;
    c.real = a.real + b.real;
    c.imag = a.imag + b.imag;
    return c;
}

fft_complex fftc_sub(fft_complex a, fft_complex b){
    fft_complex c;
    c.real = a.real - b.real;
    c.imag = a.imag - b.imag;
    return c;
}

fft_complex fftc_mul(fft_complex a, fft_complex b){
    fft_complex c;
    c.real = a.real * b.real - a.imag * b.imag;
    c.imag = a.real * b.imag + a.imag * b.real;
    return c;
}

fft_complex fftc_muls(fft_complex a, fft_float b){
    fft_complex c;
    c.real = a.real * b;
    c.imag = a.imag * b;
    return c;
}

fft_complex fftc_conj(fft_complex a){
    return fftc(a.real, -a.imag);
}

fft_complex fftc_polar(fft_float radians){
    return fftc(cos(radians), sin(radians));
}

fft_float fftc_abs(fft_complex a){
    return sqrt(a.real*a.real + a.imag*a.imag);
}

fft_float fftc_angle(fft_complex a){
    return atan2(a.imag, a.real);
}
