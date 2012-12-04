#include <stdio.h>
#include  "gsl/gsl_randist.h"
#include  "gsl/gsl_sf.h"
#include  "gsl/gsl_roots.h"
#include "kde_util.h"

typedef struct {
  int N;
  double *It;
  double *a2;
  int n;
}my_f_params;

double my_f(double x, void *params)
{
    my_f_params *p=(my_f_params *)params;
    return fixed_point(x,p->N,p->It,p->a2,p->n);
  }


int fzero(double *x_star, double N, double* It, double* a2, int n)
{

  double x=0.5,a=1;
  double b=2;
  int status;
  int max_iter=100;
  double r;
  double x0=0, x1=0.7;
  double r_expected=0;


    const gsl_root_fsolver_type *T = gsl_root_fsolver_bisection;
   gsl_root_fsolver *s = gsl_root_fsolver_alloc(T);

  gsl_function F;

  F.function=&my_f;

  my_f_params params;
  params.N=N;
  params.It=It;
  params.a2=a2;
  params.n=n;
  F.params=&params;

  int iter=0;

  gsl_root_fsolver_set (s, &F, x0,x1);

  printf ("using %s method\n",
	  gsl_root_fsolver_name (s));

  printf ("%-5s %10s %10s %10s\n",
	  "iter", "root", "err", "err(est)");

  do
    {
      iter++;
      status = gsl_root_fsolver_iterate (s);
      x0=x;
      x = gsl_root_fsolver_root (s);
      status = gsl_root_test_delta (x, x0, 0, 1e-3);

      if (status == GSL_SUCCESS)
	printf ("Converged:\n");
      
      printf ("%5d %10.7f %+10.7f %10.7f\n",
                   iter, x, x - r_expected, x - x0);
    }
  while (status == GSL_CONTINUE && iter < max_iter);
  
  *x_star=x;
  
  return status;


}
