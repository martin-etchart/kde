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

int main(int argc, char** argv)
{
	int _verbose=1;
	// Read image from text file
	int xsize;
	int ysize;
	double *data = NULL;
	//        const char * full_fname = "../../pics/trimodal2gaussian.txt";
	const char * full_fname = "../../pics/peppers.txt";
	file_read_into_array_doubles_mat(full_fname, &data, &xsize, &ysize);
	int modes = 3;

	double* Iseg = new double[xsize * ysize];
	double* thr;
	otsu(Iseg, &thr, data, xsize, ysize, modes);

	if(_verbose)
	{
		double_vector_save_to_file("data.txt", xsize*ysize, data);
		double_vector_save_to_file("data_out.txt", xsize*ysize, Iseg);
		double_vector_save_to_file("thr.txt", modes-1, thr);
		std::cout << "thr = ";
		double_vector_print(modes-1, thr);
	}

	delete[] thr;
	free(data);
	delete[] Iseg;
}
