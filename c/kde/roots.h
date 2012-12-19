#ifndef ROOTS_H_
#define ROOTS_H_

#if defined(__cplusplus) && !defined(WIN32)
extern "C" {
#endif

int fzero(double* t_star, double N, double* It, double* a2, int n);

#if defined(__cplusplus) && !defined(WIN32)
}
#endif

#endif	//ROOTS_H_
