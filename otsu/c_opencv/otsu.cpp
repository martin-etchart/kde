#include "otsu.h"
#include <stdio.h>
#include <iostream>
#include <algorithm>
#include <vector>
#include <math.h>



#ifdef __APPLE__
#include <malloc/malloc.h>
#else
#include <malloc.h>
#endif

int file_read_into_array_doubles_mat(const char *filename, double **out_data, int *xsize, int *ysize)
{
	FILE *in_file;
	in_file = fopen(filename, "r");
	double *data = NULL;

	if (in_file == NULL)
	{
		return -1;
	}
	else
	{
		fscanf(in_file, "%d", xsize);
		fscanf(in_file, "%d", ysize);
		printf("size: %d, %d\n", *xsize, *ysize);
		int total = (*xsize)*(*ysize);
		data = (double*) malloc(total * sizeof (*data));

		for (int j = 0; j < total; j++)
		{
			fscanf(in_file, "%lg", data + j);
		}
		fclose(in_file);
		*out_data = data;
	}
	return 0;
}

int double_vector_print(int l, double* v)
{
	std::cout << "[ ";
	for (int i = 0; i < l - 1; i++)
		std::cout << v[i] << ", ";
	std::cout << v[l - 1] << " ];" << std::endl;
}

int int_vector_print(int l, int* v)
{
	std::cout << "[ ";
	for (int i = 0; i < l - 1; i++)
		std::cout << v[i] << ", ";
	std::cout << v[l - 1] << " ];" << std::endl;
}

int double_vector_save_to_file(char* filename, int l, double* v)
{
	FILE *file;
	file = fopen(filename, "w");
	for (int i = 0; i < l; i++)
		fprintf(file, "%lg\n", v[i]);
	fclose(file);
}

int int_vector_save_to_file(char* filename, int l, int* v)
{
	FILE *file;
	file = fopen(filename, "w");
	for (int i = 0; i < l; i++)
		fprintf(file, "%d\n", v[i]);
	fclose(file);
}

int unique(int l, double* v_in, int* l_out, double** v_out)
{
	int _verbose = 0;

	double v[l];
	for (int i = 0; i < l; i++)
		v[i] = v_in[i];
	std::sort(v, v + l); // Sort

	if (_verbose)
	{
		std::cout << "size: " << l << std::endl;
		std::cout << "Al principio: ";
		double_vector_print(l, v_in);
		std::cout << "Despues del sort: ";
		double_vector_print(l, v);
	}

	std::vector<double> vec; // (vector_prueba, vector_prueba+(sizeof(vector_prueba)/sizeof(vector_prueba[0])));
	for (int j = 0; j < l - 1; j++)
		if (v[j] != v[j + 1])
			vec.push_back(v[j]);
	vec.push_back(v[l - 1]);

	long unsigned int len = vec.size();
	*l_out = (int) len;

	double *v_aux = NULL;
	v_aux = (double*) malloc(*l_out * sizeof (*v_aux));

	//    std::cout << "Salida del Unique: [ ";
	int i = 0;
	for (std::vector<double>::iterator it = vec.begin(); it != (vec.end() - 1); ++it)
	{
		//        std::cout << *it << ", ";
		v_aux[i] = *it;
		i++;
	}
	//    std::cout << *(vec.end()-1) << ']' << std::endl;
	v_aux[i] = *(vec.end() - 1);

	*v_out = v_aux;

	return 0;
}

int histogram(int* counts, int len, int nbins, double* data, double* bins)
{

	double step;
	int ind;
	for (int i = 0; i < len; i++)
	{
		ind = 0;
		step = bins[1] - bins[0];
		while ((data[i] >= bins[ind]+(step / 2)) & (ind + 1 < nbins))
		{
			ind++;
			step = bins[ind + 1] - bins[ind];
		}
		counts[ind]++;
	}

}

int cumsum(double* b, double* a, int N)
{
	/* b = cumsum(a) */
	b[0] = a[0];
	for (int i = 1; i < N; i++)
	{
		b[i] = b[i - 1] + a[i];
	}
}

int vector_pow(double* b, double* a, int power, int N)
{
	/* b = a^power */
	for (int i = 0; i < N; i++)
		b[i] = pow(a[i], power);
}

int divide_vectors(double* c, double* a, double* b, int N)
{
	/* c = a./b */
	for (int i = 0; i < N; i++)
		c[i] = a[i] / b[i];
}

int multiplicate_vectors(double* c, double* a, double* b, int N)
{
	/* c = a.*b */
	for (int i = 0; i < N; i++)
		c[i] = a[i] * b[i];
}

int vector_max(double* m, int* index, double* v, int N)
{
	*m = v[0];
	*index = 0;
	for (int i = 1; i < N; i++)
		if (v[i] >= *m)
		{
			*m = v[i];
			*index = i;
		}
}

int vector_flip(double* b, double* a, int N)
{
	/*
	 * b = flipud(a)
	 * b = fliplr(a)
	 */

	for (int i = 0; i < N; i++)
		b[i] = a[N - 1 - i];
}

int otsu(double* data_out, double** thr, double* data, int xsize, int ysize, int N)
{
	int _verbose=0;

	int LEVELS = 256;
	//double data_out[xsize * ysize];

	// Unique
	int unI_size;
	double* unI;
	unique(xsize*ysize, data, &unI_size, &unI);

	if(_verbose>=1)
	{
	std::cout << std::endl << "unI_c = ";
	double_vector_print(unI_size, unI);
	}

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
	}
	else if (nbins == LEVELS)
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

	free(unI);

	histogram(counts, xsize*ysize, nbins, data, bins);

	//  Display vector for matlab comparison
	if(_verbose)
	{
	std::cout << std::endl << "cBins = ";
	double_vector_print(nbins, bins);
	double_vector_save_to_file("cBins.txt", nbins, bins);
	std::cout << std::endl << "cHist = ";
	int_vector_print(nbins, counts);
	int_vector_save_to_file("cHist.txt", nbins, counts);
	}

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

	if(_verbose)
	{
	double_vector_save_to_file("P.txt", nbins, P);
	double_vector_save_to_file("w.txt", nbins, w);
	double_vector_save_to_file("mu.txt", nbins, mu);
	}

	double* w_aux;
	double* mu_aux = NULL;
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

		vector_max(&sigma2Bmax, &ind_sigma2Bmax, sigma2B, nbins - 2);
		if(_verbose)
		{
		double_vector_save_to_file("sigma2B.txt", nbins - 2, sigma2B);
		std::cout << "Maximo: " << sigma2Bmax << std::endl << "Posicion del maximo: " << ind_sigma2Bmax << std::endl;
		std::cout << "pixval: " << bins[ind_sigma2Bmax + 1] << std::endl;
		}



		for (int i = 0; i < xsize * ysize; i++)
			if (data[i] <= bins[ind_sigma2Bmax + 1])
				data_out[i] = 0.0;
			else
				data_out[i] = 1.0;

		double* th_aux;
		th_aux = new double[1];
		th_aux[0] = bins[ind_sigma2Bmax + 1];
		*thr = th_aux;

	}
	else if (N == 3)
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

	//*Iseg = data_out;
	delete[] counts;
	delete[] aux1;
	delete[] aux2;
	delete[] bins;
	delete[] w_aux;
	if (mu_aux) delete[] mu_aux;
	delete[] sigma2B;
	return 0;
}

#ifdef OTSU_WITH_OPENCV
void otsuN(IplImage* img, IplImage* img_seg, int modes, double **thr)
{
	int _verbose=1;

	int xsize = img->width;
	int ysize = img->height;

	double data [xsize * ysize];


	for (int i = 0; i < ysize; i++)
		for (int j = 0; j < xsize; j++)
		{
			CvScalar c = cvGet2D(img, i, j);
			data[i * xsize + j] = c.val[0] / 255.0;
		}

	if(_verbose)
		printf("img: %d %d %d %d\n",xsize,ysize,img->nChannels,img->depth);

	double* Iseg=new double[xsize * ysize];

	otsu(Iseg, thr, data, xsize, ysize, modes);

	for (int i = 0; i < ysize; i++)
		for (int j = 0; j < xsize; j++)
		{
			CvScalar c;
			double v = Iseg[i * xsize + j];
			c.val[0] = v * 255;
			c.val[1] = v * 255;
			c.val[2] = v * 255;
			cvSetAt(img_seg, c, i, j);
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

#endif