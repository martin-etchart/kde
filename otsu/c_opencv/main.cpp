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

int otsu(double** Iseg, double** thr, double* data, int xsize, int ysize, int N);

int main(int argc, char** argv)
{
    // Read image from text file
    int xsize;
    int ysize;
    double *data = NULL;
    //        const char * full_fname = "../../pics/trimodal2gaussian.txt";
    const char * full_fname = "../../pics/peppers.txt";
    file_read_into_array_doubles_mat(full_fname, &data, &xsize, &ysize);
    int modes = 3;

    double* Iseg;
    Iseg = new double[xsize * ysize];
    double* thr;
    otsu(&Iseg, &thr, data, xsize, ysize, modes);
    double_vector_save_to_file("data_out.txt", xsize*ysize, Iseg);

    if (modes == 3)
    {
	double_vector_save_to_file("thr.txt", 2, thr);
	std::cout << "thr = ";
	double_vector_print(2, thr);
    } else
    {
	double_vector_save_to_file("thr.txt", 1, thr);
	std::cout << "thr = ";
	double_vector_print(1, thr);
    }
}

int otsu(double** Iseg, double** thr, double* data, int xsize, int ysize, int N)
{
    int LEVELS = 256;
    double data_out[xsize * ysize];

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

    double* w_aux;
    double* mu_aux;
    double* aux1;
    double* aux2;
    double* sigma2B;

    double sigma2Bmax;
    int ind_sigma2Bmax;

    if (N == 2)
    {
	w_aux = new double[nbins - 2];
	mu_aux = new double[nbins - 2];
	aux1 = new double[nbins - 2];
	aux2 = new double[nbins - 2];
	sigma2B = new double[nbins - 2];

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

	vector_max(&sigma2Bmax, &ind_sigma2Bmax, sigma2B, nbins - 2);
	std::cout << "Maximo: " << sigma2Bmax << std::endl << "Posicion del maximo: " << ind_sigma2Bmax << std::endl;
	std::cout << "pixval: " << bins[ind_sigma2Bmax + 1] << std::endl;

	for (int i = 0; i < xsize * ysize; i++)
	    if (data[i] <= bins[ind_sigma2Bmax + 1])
		data_out[i] = 0.0;
	    else
		data_out[i] = 1.0;

	double* th_aux;
	th_aux = new double[1];
	th_aux[0] = bins[ind_sigma2Bmax + 1];
	*thr = th_aux;

    } else if (N == 3)
    {
	//		w_aux = new double[nbins ];
	//		mu_aux = new double[nbins ];
	aux1 = new double[nbins];
	aux2 = new double[nbins];
	sigma2B = new double[nbins * nbins];

	double w0[nbins * nbins], w1[nbins * nbins], w2[nbins * nbins], mu0_aux[nbins], mu2_aux[nbins], mu0[nbins * nbins], mu2[nbins * nbins];

	vector_flip(aux1, P, nbins);
	cumsum(aux2, aux1, nbins);
	vector_flip(aux1, aux2, nbins);

	for (int i = 0; i < nbins; i++)
	    for (int j = 0; j < nbins; j++)
	    {
		w0[i * nbins + j] = w[i];
		w2[i * nbins + j] = aux1[j];
	    }

	double_vector_save_to_file("w0.txt", nbins*nbins, w0);
	double_vector_save_to_file("w2.txt", nbins*nbins, w2);

	divide_vectors(mu0_aux, mu, w, nbins);
	vector_flip(aux1, Pi, nbins);
	cumsum(aux2, aux1, nbins);

	double aux3[nbins];
	vector_flip(aux3, P, nbins);
	cumsum(aux1, aux3, nbins);

	divide_vectors(aux3, aux2, aux1, nbins);
	vector_flip(mu2_aux, aux3, nbins);

	for (int i = 0; i < nbins; i++)
	    for (int j = 0; j < nbins; j++)
	    {
		mu0[i * nbins + j] = mu0_aux[i];
		mu2[i * nbins + j] = mu2_aux[j];
	    }

	double_vector_save_to_file("mu0.txt", nbins*nbins, mu0);
	double_vector_save_to_file("mu2.txt", nbins*nbins, mu2);

	for (int i = 1; i < nbins * nbins; i++)
	{
	    w1[i] = 1 - w0[i] - w2[i];
	    //	    if (w1[i] <= 0)
	    //		w1[i] = NAN;
	}

	double aux4[nbins * nbins], aux5[nbins * nbins], t1[nbins * nbins], t2[nbins * nbins], t3[nbins * nbins], t4[nbins * nbins], t34[nbins * nbins], t34_aux[nbins * nbins];

	for (int i = 0; i < nbins * nbins; i++)
	    aux4[i] = mu0[i] - mu[nbins - 1];
	vector_pow(aux5, aux4, 2, nbins * nbins);
	multiplicate_vectors(t1, w0, aux5, nbins * nbins);
	multiplicate_vectors(t3, w0, aux4, nbins * nbins);

	for (int i = 0; i < nbins * nbins; i++)
	    aux4[i] = mu2[i] - mu[nbins - 1];
	vector_pow(aux5, aux4, 2, nbins * nbins);
	multiplicate_vectors(t2, w2, aux5, nbins * nbins);
	multiplicate_vectors(t4, w2, aux4, nbins * nbins);
	for (int i = 0; i < nbins * nbins; i++)
	    t34[i] = t3[i] + t4[i];
	vector_pow(t34_aux, t34, 2, nbins * nbins);
	divide_vectors(t34, t34_aux, w1, nbins * nbins);

	for (int i = 0; i < nbins * nbins; i++)
	{
	    sigma2B[i] = t1[i] + t2[i] + t34[i];
	    if (w1[i] <= 0)
		sigma2B[i] = 0;
	    //	    if (isnan(sigma2B[i]))
	    //		sigma2B[i] = 0;
	}

	double_vector_save_to_file("sigma2B.txt", nbins*nbins, sigma2B);

	vector_max(&sigma2Bmax, &ind_sigma2Bmax, sigma2B, nbins * nbins);
	std::cout << "Maximo: " << sigma2Bmax << std::endl << "Posicion del maximo: " << ind_sigma2Bmax << std::endl;

	int k1, k2;

	k1 = (int) ind_sigma2Bmax / nbins;
	k2 = ind_sigma2Bmax % nbins;

	std::cout << "k1 = " << k1 << ", k2 = " << k2 << std::endl;

	for (int i = 0; i < xsize * ysize; i++)
	    if (data[i] <= bins[k1])
		data_out[i] = 0.0;
	    else if (data[i] > bins[k1] & data[i] <= bins[k2])
		data_out[i] = 0.5;
	    else
		data_out[i] = 1.0;

	double* thr_aux;
	thr_aux = new double[2];
	thr_aux[0] = bins[k1];
	thr_aux[1] = bins[k2];
	*thr = thr_aux;
    }

    //    for (int i = 0; i < xsize * ysize; i++)
    //	if (isnan(data_out[i])){
    //	    data_out[i] = 0;
    //	    std::cout<<"La concha de la lora hay un nan"<<std::endl;
    //	}

    *Iseg = data_out;

    return 0;
}