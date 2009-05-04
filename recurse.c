/* Functions for calculating angular momentum coupling coefficients using
   recurrsion relations. The algorithms used are those of Schulten and Gordon[1]
   augmented with the relations of Luscombe and Luben[2]. Commentary throughout
   refers to the paper of Luscombe and Luben (LL98).

   [1] K. Schulten and R. G. Gordon, J. Math. Phys. 16, 1961 (1975).
   [2] J. H. Luscombe and M. Luben, Phys. Rev. E 57, 7274 (1998).
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "recurse.h"

#define ODD(n) ((n) & 1)
#define SMALL 1.0e-10

#define SUCCESS 0
#define FAIL 1

int
LL98 (double **psi, const int two_nmin, const int two_nmax, void *params,
      double (*X) (const double, const void *),
      double (*Y) (const double, const void *),
      double (*Z) (const double, const void *),
      void (*normalize) (double *, const double, const int, const void *),
      double (*single_val) (const void *))
/* This is the generic LL98 recurssion strategy common to the 3j and 6j
   calculations. */
{
  double *_psi;
  double y, *rs;
  double nmin = two_nmin / 2.0, nmax = two_nmax / 2.0;
  int nmax_idx = (two_nmax - two_nmin) / 2;
  int ndim = nmax_idx + 1, nminus_idx=0, nplus_idx=nmax_idx, i;
  int iter_up = 1, iter_down = 1;

  *psi = malloc (ndim * sizeof (double));
  if (psi == NULL)
    {
      fprintf (stderr, "LL98: Memory allocation error (1)\n");
      return FAIL;
    }
  _psi = *psi;

  if (ndim == 1) /* Only a single value is possible, requires special handling.*/
    {
      _psi[0] = single_val (params);
      return SUCCESS;
    }

  rs = malloc (ndim * sizeof (double));
  if (rs == NULL)
    {
      fprintf (stderr, "LL98: Memory allocation error (2)\n");
      free(_psi);
      return FAIL;
    }

  /* Iterate LL98 Eq. 3 from nmin upwards unless the first term is undefined. */
  y = Y (nmin, params);

  if (fabs (y) > SMALL)
    {
      rs[0] = -X (nmin, params) / y;

      for (i = 1; i <= nmax_idx; i++)
	{
	  double n, denom;

	  if (rs[i - 1] > 1.0)
	    {
	      nminus_idx = i - 1;
	      break;
	    }

	  n = nmin + i;
	  denom = Y (n, params) + Z (n, params) * rs[i - 1];

	  if (fabs (denom) > SMALL)
	    rs[i] = -X (n, params) / denom;
	  else
	    {
	      nminus_idx = i - 1;
	      break;
	    }
	}

      /* Generate psi(n_minus-k)/psi(n_minus) == Psi_minus(n) using LL98
	 Eq. 5'. */
      if (nminus_idx > 0)
	{
	  _psi[nminus_idx-1] = rs[nminus_idx-1];

	  for (i=nminus_idx-2; i>=0; i--)
	    _psi[i] = _psi[i+1]*rs[i];
	}
    }
  else 
    {
      /* If Y is zero there are two possibilities: 
       
	 a) X != 0. In this case, first term s(nmin) is infinity because
	 psi(nmin + 1) = 0. However, psi(nmin) is not nescessarily 0 in this
	 case though. This implies we're actually in the classically allowed
	 region at nmin, and so we can later use the 3 term recursion to iterate
	 up from nmin.
       
	 b) X = 0. In this case the first term is undefined, and we're unable to
	 iterate upwards from nmin using either the 2 or 3 term recursions. */
      nminus_idx = 0;

      if (fabs(X (nmin, params)) < SMALL) 
	iter_up = 0;
    }

  /* Iterate LL98 Eq. 2 from nmax downwards, unless the first term is undefined. */
  y = Y (nmax, params);

  if (fabs (y) > SMALL)
    {
      rs[nmax_idx] = - Z (nmax, params) / y;

      for (i = nmax_idx - 1; i >= 0; i--) 
	/* This could also be (i = nmax_idx - 1; i > nminus_idx; i--), but it makes
	   no difference. */
	{
	  double n, denom;

	  if (rs[i + 1] > 1.0)
	    {
	      nplus_idx = i + 1;
	      break;
	    }

	  n = nmin + i;

	  denom = Y (n, params) + X (n, params) * rs[i + 1];

	  if (fabs (denom) > SMALL)
	    rs[i] = -Z (n, params) / denom;
	  else
	    {
	      nplus_idx = i + 1;
	      break;
	    }
	}

      /* Generate psi(n_plus+k)/psi(n_plus) == Psi_plus(n) using LL98 Eq. 4'. */
      if (nplus_idx < nmax_idx)
	{
	  _psi[nplus_idx+1] = rs[nplus_idx+1];

	  for (i=nplus_idx+2; i<=nmax_idx; i++)
	    _psi[i] = _psi[i-1]*rs[i];
	}
    }
  else
    {
      /* If Y is zero there are two possibilities: 

	 a) Z != 0. In this case, first term r(nmax) is infinity because
	 psi(nmax - 1) = 0. However, psi(nmax) is not nescessarily 0 in this
	 case though. This implies we're actually in the classically allowed
	 region at nmax, and so we can later use the 3 term recursion to iterate
	 up from nmin.

	 b) Z = 0. In this case the first term is undefined, and we're unable to
	 iterate upwards from nmin using either the 2 or 3 term recursions. */
      nplus_idx = nmax_idx;

      if (fabs (Z (nmax, params)) < SMALL) 
	iter_down = 0;
    }

  free (rs);

  /* Iterate in the classical region using three term recursion LL98 Eq. 1.  */
  if (iter_up) /* Iterate upwards from nminus, chosing nc = nplus. */
    {
      double a;
      int iter_up_start_idx;

      /* Note that this initialization stuff can't be done inside the logic of
	 iterating LL98 Eq. 3 above, since it can potentially be clobbered during
	 the subsequent iteration of LL98 Eq. 4 if that section was also to
	 contain initialization logic for iterating downwards in the classical
	 region below. Really, tempting though it is, don't move this earlier. */
      if (nminus_idx < 2)
	{
	  _psi[0] = 1.0;
	  _psi[1] = -Y(nmin, params) / X(nmin, params); /* Since psi(nmin - 1) = 0 */
	  iter_up_start_idx = 2;
	}
      else
	{
	  _psi[nminus_idx] = 1.0;
	  iter_up_start_idx = nminus_idx + 1;
	}
      
      for (i = iter_up_start_idx; i <= nplus_idx; i++)
	{
	  double nn = nmin - 1.0 + i;	/* n - 1 */
	  _psi[i] = -(Y (nn, params) * _psi[i - 1] +
			Z (nn, params) * _psi[i - 2]) / X (nn, params);
	}

      /* Since we choose nc=nplus, Psi_plus(nc)=1, and we multiply
	 Psi_minus(nmin...nplus) by Psi_plus(nc)/Psi_minus(nc) ==
	 1/Psi_minus(n_plus) to give us Psi_plus(nmin...nplus). */
      a = 1.0 / _psi[nplus_idx];
      
      for (i = 0; i <= nplus_idx; i++)
	_psi[i] *= a;
      
      normalize (_psi, nmin, nmax_idx, params);
      return SUCCESS;
    }

  if (iter_down) /* Iterate downwards from nplus, chosing nc = nminus. */
    {
      double a;
      int iter_down_start_idx;

      /* Note that this initialization stuff could be done inside the logic of
	 iterating LL98 Eq. 2 above. However following that design leads to some
	 rather obscure corner cases and errors, so it's cleaner to do it
	 here. Really, don't move it. */
      if (nplus_idx > nmax_idx - 2)
	{
	  _psi[nplus_idx] = 1.0;
	  _psi[nplus_idx - 1] = -Y(nmax, params) / Z(nmax, params);
	  iter_down_start_idx = nplus_idx - 2;
	}
      else
	{
	  _psi[nplus_idx] = 1.0;
	  iter_down_start_idx = nplus_idx - 1;
	}

      for (i = iter_down_start_idx; i >= nminus_idx; i--)
	{
	  double nn = nmin + 1.0 + i;	/* n + 1 */
	  _psi[i] = -(X (nn, params) * _psi[i + 2] +
			Y (nn, params) * _psi[i + 1]) / Z (nn, params);
	}
	  
      /* Since we choose nc=nminus, Psi_minus(nc)=1, and we multiply
	 Psi_plus(nminus...nmax) by Psi_minus(nc)/Psi_plus(nc) ==
	 1/Psi_plus(n_plus) to give us Psi_minus(nminus...nmax). */
      a = 1.0 / _psi[nminus_idx];
      
      for (i = nmax_idx; i >= nminus_idx; i--)
	_psi[i] *= a;
      
      normalize (_psi, nmin, nmax_idx, params);
      return SUCCESS;
    }

  fprintf (stderr, "LL98: Could not iterate in either direction\n");
  return FAIL;

}

/* Specific functions for calculation of 3j coefficients using j recurrsion -
   first column of Table 1 of LL98. */
typedef struct params_3j_j
{
  int two_j2, two_j3, two_m2, two_m3, two_jmin, two_jmax;
  double j2, j3, m2, m3;
} params_3j_j;

double
A (const double j, const double j2, const double j3,
   const double m2, const double m3)
{
  double a = j * j;
  double b = j2 - j3;
  double c = j2 + j3 + 1.0;
  double d = m2 + m3;

  return sqrt ((a - b * b) * (c * c - a) * (a - d * d));
}

double
B (const double j, const double j2, const double j3,
   const double m2, const double m3)
{
  double a = m2 + m3;
  double b = j2 * (j2 + 1.0) - j3 * (j3 + 1.0);
  double c = (m2 - m3) * j * (j + 1.0);

  return (2.0 * j + 1.0) * (a * b - c);
}

double
X_3j_j (const double j, const void *params)
{
  params_3j_j *p = (params_3j_j *) params;
  return j * A (j + 1.0, p->j2, p->j3, p->m2, p->m3);
}

double
Y_3j_j (const double j, const void *params)
{
  params_3j_j *p = (params_3j_j *) params;
  return B (j, p->j2, p->j3, p->m2, p->m3);
}

double
Z_3j_j (const double j, const void *params)
{
  params_3j_j *p = (params_3j_j *) params;
  return (j + 1.0) * A (j, p->j2, p->j3, p->m2, p->m3);
}

void
normalize_3j_j (double *f, const double jmin, const int jmax_idx,
		const void *params)
{
  params_3j_j *p = (params_3j_j *) params;
  double a = 0.0, phase;
  int i;

  for (i = 0; i <= jmax_idx; i++)
    {
      double j = jmin + i;
      double ff = f[i];
      a += (2.0 * j + 1.0) * ff * ff;
    }

  a = 1.0 / sqrt (a);

  if (ODD ((p->two_j2 - p->two_j3 + p->two_m2 + p->two_m3) / 2))
    phase = -1.0;
  else
    phase = 1.0;

  if ((f[jmax_idx] / phase) < 0)
    a = -a;

  for (i = 0; i <= jmax_idx; i++)
    f[i] *= a;
}

double
single_val_3j_j (const void *params)
{
  params_3j_j *p = (params_3j_j *) params;
  double a = 1.0 / sqrt(p->two_jmin + 1.0);
  
  if (ODD ((p->two_jmin + p->two_m2 + p->two_m3) / 2))
    return -a;
  else
    return a;
}

void
wigner3j_family_j (const int two_j2, const int two_j3,
		   const int two_m2, const int two_m3,
		   double **family, int *two_jmin, int *two_jmax)
{
  // TODO: Add checking for vald inputs!
  params_3j_j p;
  int a = abs (two_j2 - two_j3);
  int b = abs (two_m2 + two_m3);
  
  *two_jmin = a > b ? a : b;
  *two_jmax = two_j2 + two_j3;

  p.two_j2 = two_j2;
  p.two_j3 = two_j3;
  p.two_m2 = two_m2;
  p.two_m3 = two_m3;

  p.j2 = two_j2 / 2.0;
  p.j3 = two_j3 / 2.0;
  p.m2 = two_m2 / 2.0;
  p.m3 = two_m3 / 2.0;

  p.two_jmin=*two_jmin;
  p.two_jmax=*two_jmax;
  
  LL98 (family, *two_jmin, *two_jmax, &p, X_3j_j, Y_3j_j, Z_3j_j,
	normalize_3j_j, single_val_3j_j);
}

/* End of specifics for 3j calculation by j recurrsion. */


#undef ODD
#undef SMALL
