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

void otsuN(IplImage* img, IplImage* img_seg, int modes, double **thr);

int main(int argc, char** argv)
{
	const char * full_fname = "../../pics/trimodal2gaussian.png";
	IplImage* img = cvLoadImage(full_fname);

	int modes = 3;
	double* thr;

	IplImage* img_seg_cv=cvCreateImage(cvGetSize(img),img->depth,1);
	otsuN(img, img_seg_cv, modes, &thr);

	cvSaveImage("out.png",img_seg_cv);


	cvReleaseImage(&img_seg_cv);
	cvReleaseImage(&img);
}

void otsuN(IplImage* img, IplImage* img_seg, int modes, double **thr)
{
	int _verbose=0;

	int xsize = img->width;
	int ysize = img->height;

	double data [xsize * ysize];

	for (int i = 0; i < ysize; i++)
		for (int j = 0; j < xsize; j++)
		{
			CvScalar c = cvGet2D(img, i, j);
			data[i * ysize + j] = c.val[0] / 255;
		}


	double* Iseg;
	Iseg = new double[xsize * ysize];

	otsu(Iseg, thr, data, xsize, ysize, modes);

	for (int i = 0; i < ysize; i++)
		for (int j = 0; j < xsize; j++)
		{
		CvScalar c;
			double v=Iseg[i * ysize + j];
			c.val[0] =v*255;
			c.val[1] =v*255;
			c.val[2] =v*255;
			cvSetAt(img_seg,c, i, j);
		}

	if(_verbose)
	{
		double_vector_save_to_file("data.txt", xsize*ysize, data);
		double_vector_save_to_file("data_out.txt", xsize*ysize, Iseg);
		double_vector_save_to_file("thr.txt", modes-1, *thr);
		std::cout << "thr = ";
		double_vector_print(modes-1, *thr);
	}

	delete[] Iseg;
	//delete[] data;
}