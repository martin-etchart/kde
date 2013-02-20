/* 
 * File:   main.cpp
 * Author: mat
 *
 * Created on February 13, 2013, 8:26 PM
 */

#include <cstdlib>
#include <iostream>
#include "utils.h"
#include <math.h>
#include <vector>

using namespace std;

int LEVELS = 256;

/*
 * 
 */

int main(int argc, char** argv)
{

    // Read image from text file
    int xsize;
    int ysize;
    double *data = NULL;
    const char * full_fname = "../../pics/trimodal2gaussian.txt";
    //const char * full_fname = "../../pics/peppers.txt";
    file_read_into_array_doubles_mat(full_fname, &data, &xsize, &ysize);
    // Unique
    int unI_size;
    double* unI;
    unique(xsize*ysize, data, &unI_size, &unI);
    std::cout << std::endl << "unI_c = ";
    double_vector_print(unI_size, unI);

    // Histogram
    int* counts;
    double* bins;

    int nbins;
    if (unI_size < LEVELS)
	nbins = unI_size;
    else
	nbins = LEVELS;

    if (nbins < LEVELS)
    {
	counts = new int[unI_size];
	for (int i = 0; i < unI_size; i++)
	    counts[i] = 0;
	bins = new double[unI_size];
	for (int i = 0; i < unI_size; i++)
	    bins[i] = unI[i];
    } else if (nbins == LEVELS)
    {
	counts = new int[LEVELS];
	for (int i = 0; i < LEVELS; i++)
	    counts[i] = 0;
	bins = new double[LEVELS];
	double step = (LEVELS - 1) / LEVELS / 2;
	bins[0] = step;
	for (int i = 1; i < LEVELS; i++)
	    bins[i] = bins[i - 1] + step;
    }

    // free aca ya puedo liberar la memoria de unI

    histogram(counts, xsize*ysize, nbins, data, bins);

    //  Display vector for matlab comparison
    std::cout << std::endl << "cBins = ";
    double_vector_print(nbins, bins);
    double_vector_save_to_file("cBins.txt", nbins, bins);
    std::cout << std::endl << "cHist = ";
    int_vector_print(nbins, counts);
    int_vector_save_to_file("cHist.txt", nbins, counts);

    double P[nbins];
    for (int i = 0; i < nbins; i++)
	P[i] = (double) counts[i] / (xsize * ysize);

    // cumsum
    double w[nbins], mu[nbins], Pi[nbins];
    Pi[0] = P[0];
    for (int i = 1; i < nbins; i++)
	Pi[i] = (i + 1) * P[i];

    cumsum(w, P, nbins);
    cumsum(mu, Pi, nbins);

    double_vector_save_to_file("P.txt", nbins, P);
    double_vector_save_to_file("w.txt", nbins, w);
    double_vector_save_to_file("mu.txt", nbins, mu);

    double w_aux[nbins - 2], mu_aux[nbins - 2], aux1[nbins - 2], aux2[nbins - 2], sigma2B[nbins - 2];
    for (int i = 0; i < nbins - 2; i++)
    {
	w_aux[i] = w[i + 1];
	mu_aux[i] = mu[i + 1];
	aux1[i] = mu[nbins - 1] * w_aux[i] - mu_aux[i];
    }

    vector_pow(aux2, aux1, 2, nbins - 2);
    divide_vectors(aux1, aux2, w_aux, nbins - 2);
    for (int i = 0; i < nbins - 2; i++)
	aux2[i] = 1 - w_aux[i];
    divide_vectors(sigma2B, aux1, aux2, nbins - 2);

    double_vector_save_to_file("sigma2B.txt", nbins - 2, sigma2B);

    double sigma2Bmax;
    int ind_sigma2Bmax;
    vector_max(&sigma2Bmax, &ind_sigma2Bmax, sigma2B, nbins - 2);
    std::cout << "Maximo: " << sigma2Bmax << std::endl << "Posicion del maximo: " << ind_sigma2Bmax << std::endl;
    std::cout << "pixval: " << bins[ind_sigma2Bmax + 1] << std::endl;
    double data_out[xsize * ysize];
    for (int i = 0; i < xsize * ysize; i++)
	if (data[i] <= bins[ind_sigma2Bmax + 1])
	    data_out[i] = 0;
	else
	    data_out[i] = 1;

    double_vector_save_to_file("data_out.txt", xsize*ysize, data_out);

    return 0;
}