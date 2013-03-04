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
#include <cv.h>
#include "highgui.h"

using namespace std;



int main(int argc, char** argv)
{
	const char * full_fname = "../../pics/trimodal2gaussian.png";
	IplImage* img = cvLoadImage(full_fname);

	int xsize = img->width;
	int ysize = img->height;

	double data [xsize * ysize];

	for (int i = 0; i < ysize; i++)
		for (int j = 0; j < xsize; j++)
		{
			CvScalar c = cvGet2D(img, i, j);
			data[i * ysize + j] = c.val[0] / 255;
		}

	cvReleaseImage(&img);

	double_vector_save_to_file("data.txt", xsize*ysize, data);

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
	}
	else
	{
		double_vector_save_to_file("thr.txt", 1, thr);
		std::cout << "thr = ";
		double_vector_print(1, thr);
	}

}

