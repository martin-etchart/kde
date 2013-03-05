/*
 * File:   main.cpp
 * Author: mat
 *
 * Created on February 13, 2013, 8:26 PM
 */

#include <cstdlib>
#include <iostream>
#include "otsu.h"
#include <math.h>
#include <vector>

using namespace std;

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
	delete[] thr;
}

