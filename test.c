#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <gsl/gsl_sf_coupling.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf.h>

#include "recurse.h"

#define __SUCCESS 0
#define __FAILURE 1

/* void */
/* gough_test (void) */
/* { */
/*   int tj; */
/*   int N = 1e9; */
/*   gsl_set_error_handler_off(); */

/*   for (tj = 0; tj < N; tj = (tj +1)*100-1) */
/*     { */
/*       gsl_sf_result rc, rw; */
/*       int ec = gsl_sf_coupling_3j_e (0, tj, tj, 0, tj, -tj, &rc); */
/*       int ew = gsl_sf_wigner_3j_e (0, tj, tj, 0, tj, -tj, &rw); */
/*       double exact = 1.0 / sqrt(tj + 1.0); */

/*       printf("j=%d/2 c3j(0 j j| 0 j -j)= % .18e %.5e (status=%d)\n", tj, rc.val, rc.err, ec); */
/*       printf("j=%d/2 w3j(0 j j| 0 j -j)= % .18e %.5e (status=%d)\n", tj, rw.val, rw.err, ew); */
/*       printf("j=%d/2 exact             = % .18e\n",tj, exact); */
/*       printf("j=%d/2 wigner-exact      = % .18e\n",tj, rw.val - exact); */
/*     } */
/* } */

int
check_3j_family_j_exact_1 (const int two_jmax)
{
  int two_j;

  for (two_j = 0; two_j <= two_jmax; two_j = (two_j + 1) * 100 - 1)
    {
      double exact = 1.0 / sqrt (two_j + 1.0);
      double gsl = gsl_sf_coupling_3j (two_j, two_j, 0, two_j, -two_j, 0);
      double *a;
      int tjmax, tjmin;
      size_t dim;

      dim = wigner3j_family_j_dim (two_j, two_j, two_j, -two_j);
      a = malloc (dim * sizeof (double));
      if (a == NULL)
	{
	  fprintf (stderr, "%s:%d memory allocation error\n", __FILE__,
		   __LINE__);
	  return __FAILURE;
	}

      wigner3j_family_j (two_j, two_j, two_j, -two_j, a, &tjmax, &tjmin);

      printf ("two_j: %d\trecursive: % .18e\t exact: % .18e\tdiff: % .18e\n",
	      two_j, a[0], exact, a[0] - exact);

      printf ("two_j: %d\t gsl: % .18e\t exact: % .18e\tdiff: % .18e\n",
	      two_j, gsl, exact, gsl - exact);


      free (a);
    }

  return __SUCCESS;
}

int
check_3j_family_j_exact_2 (const int two_jmax)
{
  int two_j;

  printf ("two_j\trecursive\tgsl\texact\trec-ex\tgsl-ex\n");

  for (two_j = 0; two_j <= two_jmax; two_j = (two_j + 1) * 100)
    {
      double exact = 1.0 / sqrt (two_j + 1.0);
      double gsl = gsl_sf_coupling_3j (two_j, two_j, 0, 0, 0, 0);
      double *a;
      int tjmax, tjmin;
      size_t dim;

      dim = wigner3j_family_j_dim (two_j, two_j, 0, 0);
      a = malloc (dim * sizeof (double));
      if (a == NULL)
	{
	  fprintf (stderr, "%s:%d memory allocation error\n", __FILE__,
		   __LINE__);
	  return __FAILURE;
	}

      wigner3j_family_j (two_j, two_j, 0, 0, a, &tjmax, &tjmin);

      printf ("%d\t%g\t%g\t%g\t%g\t%g\n", two_j, a[0], gsl, exact,
	      a[0] - exact, gsl - exact);

      free (a);
    }

  return __SUCCESS;
}

int
hammer_3j_j (const int two_jmax)
{
  double *a;
  int two_j2, two_j3, two_m2, two_m3, tjmax, tjmin;

  for (two_j2 = 0; two_j2 <= two_jmax; two_j2++)
    for (two_m2 = -two_j2; two_m2 <= two_j2; two_m2 += 2)
      for (two_j3 = 0; two_j3 <= two_jmax; two_j3++)
	for (two_m3 = -two_j3; two_m3 <= two_j3; two_m3 += 2)
	  {
	    size_t dim;

	    dim = wigner3j_family_j_dim (two_j2, two_j3, two_m2, two_m3);
	    a = malloc (dim * sizeof (double));
	    if (a == NULL)
	      {
		fprintf (stderr, "%s:%d memory allocation error\n", __FILE__,
			 __LINE__);
		return __FAILURE;
	      }
	    printf ("%d %d %d %d\n", two_j2, two_j3, two_m2, two_m3);
	    wigner3j_family_j (two_j2, two_j3, two_m2, two_m3, a, &tjmax,
			       &tjmin);
	    free (a);
	  }

  return __SUCCESS;
}

int
check_3j_family_j (const int two_j1, const int two_j2,
		   const int two_m1, const int two_m2)
{
  double *a;
  int two_jmin, two_jmax, imax, i;
  size_t dim;

  dim = wigner3j_family_j_dim (two_j1, two_j2, two_m1, two_m2);
  a = malloc (dim * sizeof (double));
  if (a == NULL)
    {
      fprintf (stderr, "%s:%d memory allocation error\n", __FILE__, __LINE__);
      return __FAILURE;
    }

  wigner3j_family_j (two_j1, two_j2, two_m1, two_m2, a, &two_jmin,
		     &two_jmax);

  imax = (two_jmax - two_jmin) / 2;

  printf ("two_j1: %d two_j2: %d two_m1: %d two_m2: %d\n",
	  two_j1, two_j2, two_m1, two_m2);

  for (i = 0; i <= imax; i++)
    {
      int j = two_jmin + i;
      printf ("\ttwo_j3: %d\t %g\n", j, a[i]);
    }

  free (a);

  return __SUCCESS;
}

int
check_3j_family_j_gsl (const int two_j1, const int two_j2,
		       const int two_m1, const int two_m2)
{
  double *a;
  int two_jmin, two_jmax, imax, i;
  int two_m3 = -(two_m1 + two_m2);
  int ret;
  size_t dim;

  dim = wigner3j_family_j_dim (two_j1, two_j2, two_m1, two_m2);
  a = malloc (dim * sizeof (double));
  if (a == NULL)
    {
      fprintf (stderr, "%s:%d memory allocation error\n", __FILE__, __LINE__);
      return __FAILURE;
    }

  ret = wigner3j_family_j (two_j1, two_j2, two_m1, two_m2,
			   a, &two_jmin, &two_jmax);

  if (ret)
    {
      printf ("Invalid angular momenta\n");
      free (a);
      return __FAILURE;
    }

  imax = (two_jmax - two_jmin) / 2;

  printf ("two_j1: %d two_j2: %d two_m1: %d two_m2: %d\n",
	  two_j1, two_j2, two_m1, two_m2);

  for (i = 0; i <= imax; i++)
    {
      int two_j3 = two_jmin + i * 2;

      double gsl = gsl_sf_coupling_3j (two_j1, two_j2, two_j3,
				       two_m1, two_m2, two_m3);
      printf ("\ttwo_j3: %d\t\twigner: %g\tgsl: %g\tdiff: %g\n",
	      two_j3, a[i], gsl, a[i] - gsl);
    }

  free (a);

  return __SUCCESS;
}

int
check_3j_family_m_gsl (const int two_j1, const int two_j2,
		       const int two_j3, const int two_m1)
{
  double *a;
  int two_mmin, two_mmax, imax, i;
  size_t dim;

  dim = wigner3j_family_m_dim (two_j1, two_j2, two_j3, two_m1);
  a = malloc (dim * sizeof (double));
  if (a == NULL)
    {
      fprintf (stderr, "%s:%d memory allocation error\n", __FILE__, __LINE__);
      return __FAILURE;
    }

  wigner3j_family_m (two_j1, two_j2, two_j3, two_m1, a, &two_mmin,
		     &two_mmax);

  imax = (two_mmax - two_mmin) / 2;

  printf ("two_j1: %d two_j2: %d two_j3: %d two_m1: %d\n",
	  two_j1, two_j2, two_j3, two_m1);

  for (i = 0; i <= imax; i++)
    {
      int two_m = two_mmin + i * 2;
      double gsl = gsl_sf_coupling_3j (two_j1, two_j2, two_j3,
				       two_m1, two_m, -two_m1 - two_m);
      printf ("\ttwo_m: %d\t\twigner: %g\tgsl: %g\tdiff: %g\n",
	      two_m, a[i], gsl, a[i] - gsl);
    }

  free (a);

  return __SUCCESS;
}


int
main ()
{
/*   double *a; */
/*   int i, two_jmin, two_jmax, imax; */

  //hammer_3j_j (2);

  check_3j_family_j_gsl (1, 1, -1, 1);
  check_3j_family_j_gsl (9, 1, -1, 1);
  check_3j_family_j_gsl (9, 8, -1, 0);
  check_3j_family_j_gsl (90, 80, -30, 10);

  printf ("These should all be NAN, since j is half integer, m is 0\n");
  check_3j_family_j_gsl (99, 99, 0, 0);

  check_3j_family_j_gsl (99, 99, 1, 1);

  check_3j_family_j_gsl (100, 100, 0, 0);

  check_3j_family_j_exact_1 (1e4);
  //check_3j_family_j_exact_2(10000);

/*   free (a); */
/*   wigner3j_family_j (2, 2, -2, 2, a, &two_jmin, &two_jmax); */
/*   free (a); */
/*   wigner3j_family_j (1, 2, -1, 0, a, &two_jmin, &two_jmax); */

  check_3j_family_m_gsl (1, 1, 2, -1);
  check_3j_family_m_gsl (9, 8, 9, 3);
  check_3j_family_m_gsl (90, 80, 90, 10);
  check_3j_family_m_gsl (180, 180, 90, 10);

  return 0;
}
