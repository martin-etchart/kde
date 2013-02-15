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

    // Histogram
    double bins[11] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
    int counts[11] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    
    histogram(counts, xsize*ysize, data, bins);
        
    //  Display vector for matlab comparison
    std::cout << "cBins = ";
    double_vector_print(11,bins);    
    std::cout << "cHist = ";
    int_vector_print(11,counts);
    
    // Unique
    double vector_prueba[] = {2, 4.5, 2.2, 512, 9.3, 4.5, 512, 2, 134};
    int size_out;
    double* vector_prueba_out = unique(sizeof(vector_prueba)/sizeof(vector_prueba[0]),vector_prueba, &size_out);

    std::cout<<"ASDF: "<<size_out<<std::endl;
    std::cout<<"ASDF: "<<vector_prueba_out<<std::endl;
    std::cout<<"ASDF: "<<vector_prueba_out[1]<<std::endl;
    
    double_vector_print(size_out,vector_prueba_out);

    
    return 0;
}