/* 
 * File:   utils.h
 * Author: mat
 *
 * Created on February 13, 2013, 8:33 PM
 */

#ifndef UTILS_H
#define	UTILS_H

#include <stddef.h>

int file_read_into_array_doubles_mat(const char *filename, double **out_data, int *xsize, int *ysize);
int double_vector_print(int l, double* v);
int int_vector_print(int l, int* v);
double* unique(int l, double* v, int* l_out);
int histogram(int* counts, int len, double* data, double* bins);

#endif	/* UTILS_H */

