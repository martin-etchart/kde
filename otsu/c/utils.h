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
int double_vector_save_to_file(char* filename, int l, double* v);
int int_vector_save_to_file(char* filename, int l, int* v);
int unique(int l, double* v, int* l_out, double** v_out);
int histogram(int* counts, int len, int nbins, double* data, double* bins);

int cumsum(double* b, double* a, int N);
int vector_pow(double* b, double* a, int power, int N);
int divide_vectors(double* c, double* a, double* b, int N);

#endif	/* UTILS_H */
