#include "utils.h"
#include <stdio.h>
#include <iostream>
#include <algorithm>
#include <vector>

#ifdef __APPLE__
#include <malloc/malloc.h>
#else
#include <malloc.h>
#endif

int file_read_into_array_doubles_mat(const char *filename, double **out_data, int *xsize, int *ysize) {
    FILE *in_file;
    in_file = fopen(filename, "r");
    double *data = NULL;

    if (in_file == NULL) {
        return -1;
    } else {
        fscanf(in_file, "%d", xsize);
        fscanf(in_file, "%d", ysize);
        printf("size: %d, %d\n", *xsize, *ysize);
        int total = (*xsize)*(*ysize);
        data = (double*) malloc(total * sizeof (*data));

        for (int j = 0; j < total; j++) {
            fscanf(in_file, "%lg", data + j);
        }
        fclose(in_file);
        *out_data = data;
    }
    return 0;
}

int double_vector_print(int l, double* v) {
    std::cout << "[ ";
    for (int i = 0; i < l - 1; i++)
        std::cout << v[i] << ", ";
    std::cout << v[l - 1] << " ]" << std::endl;
}

int int_vector_print(int l, int* v) {
    std::cout << "[ ";
    for (int i = 0; i < l - 1; i++)
        std::cout << v[i] << ", ";
    std::cout << v[l - 1] << " ]" << std::endl;
}

double* unique(int l, double* v, int* l_out) {
    std::cout << "Al principio: ";
    double_vector_print(l, v);

    std::sort(v, v + l); // Sort
    std::cout << "Despues del sort: ";
    double_vector_print(l, v);

    std::vector<double> vec; // (vector_prueba, vector_prueba+(sizeof(vector_prueba)/sizeof(vector_prueba[0])));
    for (int j = 0; j < l - 1; j++)
        if (v[j] != v[j + 1])
            vec.push_back(v[j]);
    vec.push_back(v[l - 1]);
    
    long unsigned int len = vec.size();
    *l_out=(int)len;
    double v_out[len];
    
    std::cout << "Salida del Unique: [ ";
    int i=0;
    for (std::vector<double>::iterator it = vec.begin(); it != (vec.end()-1); ++it) {
        std::cout << *it << ", ";
        v_out[i]=*it;
        i++;
    }
    std::cout << *(vec.end()-1) << ']' << std::endl;
    v_out[i]=*(vec.end()-1);
    
    return 0;
}

int histogram(int* counts, int len, double* data, double* bins) {

    double step = bins[1] - bins[0];
    int ind;
    for (int i = 0; i < len; i++) {
        ind = 0;
        while (data[i] >= bins[ind]+(step / 2))
            ind++;
        counts[ind]++;
    }
}