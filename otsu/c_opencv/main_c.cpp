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
	otsu(Iseg, &thr, data, xsize, ysize, modes);
	double_vector_save_to_file("data_out.txt", xsize*ysize, Iseg);

	if (modes == 3)
	{
		double_vector_save_to_file("thr.txt", 2, thr);
		std::cout << "thr = ";
		double_vector_print(2, thr);
	}
	else
	{
		double_vector_save_to_file("thr.txt", 1, thr);
		std::cout << "thr = ";
		double_vector_print(1, thr);
	}

	delete[] thr;
	free(data);
	delete[] Iseg;
}
