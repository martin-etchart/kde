#include "utils.h"
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
    } else
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
    //    std::cout << "Al principio: ";
    //    double_vector_print(l, v);
    double v[l];
    for (int i = 0; i < l; i++)
	v[i] = v_in[i];
    std::sort(v, v + l); // Sort
    //    std::cout << "Despues del sort: ";
    //    double_vector_print(l, v);

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