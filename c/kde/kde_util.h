#ifndef _KDE_UTIL
#define _KDE_UTIL

#include <stdlib.h>
#include <math.h>
#include <malloc/malloc.h>

#define VERBOSE
#ifdef VERBOSE
#define xml_out(...) fprintf ( stdout , __VA_ARGS__ )
#else
#define xml_out(...) 
#endif

#define XML_IN xml_out("<%s>\n",  __FUNCTION__)
#define XML_IN_T(str) xml_out("<%s title=\"%s\">\n", __FUNCTION__ , str)
#define XML_OUT xml_out("</%s>\n", __FUNCTION__)

int print_vec(double *v, char* title, int start, int end);
int file_read_into_array_doubles(const char *filename , double **out_data, int *length);
void kde_dct_fftw(double *in, int n, double* out);

#endif //_KDE_UTIL
