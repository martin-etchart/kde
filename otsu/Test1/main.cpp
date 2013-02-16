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

/*
 * 
 */
int main(int argc, char** argv) {

    // Read image from text file
    int xsize;
    int ysize;
    double *data = NULL;
    const char * full_fname = "../peppers.txt";
    file_read_into_array_doubles_mat(full_fname, &data, &xsize, &ysize);

    // Unique
    int unI_size;
    double* unI;
    unique(xsize*ysize, data, &unI_size, &unI);
    std::cout << "unI_c = "; double_vector_print(unI_size, unI);

    // Histogram
    double bins[11] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
    int counts[11] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    histogram(counts, xsize*ysize, data, bins);

    //  Display vector for matlab comparison
    std::cout << "cBins = ";
    double_vector_print(11, bins);
    std::cout << "cHist = ";
    int_vector_print(11, counts);

    return 0;
}