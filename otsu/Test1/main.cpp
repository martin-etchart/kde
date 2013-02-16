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

int main(int argc, char** argv) {

    // Read image from text file
    int xsize;
    int ysize;
    double *data = NULL;
    const char * full_fname = "../trimodal2gaussian.txt";
    //const char * full_fname = "../peppers.txt";
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

    if (nbins < LEVELS) {
        counts = new int[unI_size];
        for (int i = 0; i < unI_size; i++)
            counts[i] = 0;
        bins = new double[unI_size];
        for (int i = 0; i < unI_size; i++)
            bins[i] = unI[i];
    } else if (nbins == LEVELS) {
        counts = new int[LEVELS];
        for (int i = 0; i < LEVELS; i++)
            counts[i] = 0;
        bins = new double[LEVELS];
        double step = (LEVELS - 1) / LEVELS / 2;
        bins[0] = step;
        for (int i = 1; i < LEVELS; i++)
            bins[i] = bins[i - 1] + step;
    }

    histogram(counts, xsize*ysize, nbins, data, bins);

    //  Display vector for matlab comparison
    std::cout << std::endl << "cBins = ";
    double_vector_print(nbins, bins);
    double_vector_save_to_file("cBins.txt", nbins, bins);
    std::cout << std::endl << "cHist = ";
    int_vector_print(nbins, counts);
    int_vector_save_to_file("cHist.txt", nbins, counts);

    return 0;
}