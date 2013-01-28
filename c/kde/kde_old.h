#ifndef KDE_OLD_H_
#define KDE_OLD_H_

#if defined(__cplusplus) && !defined(WIN32)
extern "C" {
#endif


#include <complex.h>



int fft(double *data, int length, double complex *fft_data);

int ifft(double *data, int length, double *ifft_data);

int dct1d(double *data, int length, double *dct_data);

int idct1d(double *data, int length, double *dct_data);


#if defined(__cplusplus) && !defined(WIN32)
}
#endif

#endif //KDE_OLD_H_
