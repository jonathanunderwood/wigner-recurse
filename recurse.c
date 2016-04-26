/* This file is part of wigner-recurse, code for calculating Wigner 3j
 * and 6j symbols using recursion.
 *
 * Copyright (C) 2008-2016 Jonathan * G. Underwood.
 *
 * wigner-recurse is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * wigner-recurse is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with wigner-recurse.  If not, see
 * <http://www.gnu.org/licenses/>.
 */


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
#include <stdbool.h>

#define INLINE
#include "recurse.h"

#define __ODD(n) ((n) & 1)
#define __EVEN(n) (!(__ODD(n)))
#define __MAX(a,b) ((a) > (b) ? (a) : (b))
#define __MIN(a,b) ((a) < (b) ? (a) : (b))
#define __SMALL 1.0e-15

#define __SUCCESS 0
#define __FAILURE 1

static inline int
LL98 (double psi[], const int two_nmin, const int two_nmax, void *params,
      double (*X) (const double, const void *),
      double (*Y) (const double, const void *),
      double (*Z) (const double, const void *),
      void (*normalize) (double *, const double, const int, const void *),
      double (*single_val) (const void *))
/* This is the generic LL98 recurssion strategy common to the 3j and 6j
   calculations. */
{
  int nmax_idx = (two_nmax - two_nmin) / 2;
  int ndim = nmax_idx + 1, nminus_idx = 0, nplus_idx = nmax_idx, i;
  bool iter_up = true, iter_down = true;
  double y;
  double nmin = two_nmin / 2.0, nmax = two_nmax / 2.0;
  double rs[ndim];

  if (ndim == 1)		/* Only a single value is possible, requires special handling. */
    {
      psi[0] = single_val (params);
      return __SUCCESS;
    }

  /* Iterate LL98 Eq. 3 from nmin upwards unless the first term is undefined. */
  y = Y (nmin, params);

  if (fabs (y) > __SMALL)
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

	  if (fabs (denom) > __SMALL)
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
	  psi[nminus_idx - 1] = rs[nminus_idx - 1];

	  for (i = nminus_idx - 2; i >= 0; i--)
	    psi[i] = psi[i + 1] * rs[i];
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

      if (fabs (X (nmin, params)) < __SMALL)
	iter_up = false;
    }

  /* Iterate LL98 Eq. 2 from nmax downwards, unless the first term is undefined. */
  y = Y (nmax, params);

  if (fabs (y) > __SMALL)
    {
      rs[nmax_idx] = -Z (nmax, params) / y;

      for (i = nmax_idx - 1; i > nminus_idx; i--)
	{
	  double n, denom;

	  if (rs[i + 1] > 1.0)
	    {
	      nplus_idx = i + 1;
	      break;
	    }

	  n = nmin + i;

	  denom = Y (n, params) + X (n, params) * rs[i + 1];

	  if (fabs (denom) > __SMALL)
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
	  psi[nplus_idx + 1] = rs[nplus_idx + 1];

	  for (i = nplus_idx + 2; i <= nmax_idx; i++)
	    psi[i] = psi[i - 1] * rs[i];
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

      if (fabs (Z (nmax, params)) < __SMALL)
	iter_down = false;
    }

  /* Iterate in the classical region using three term recursion LL98 Eq. 1.  */
  if (iter_up)			/* Iterate upwards from nminus, chosing nc = nplus. */
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
	  psi[0] = 1.0;
	  psi[1] = -Y (nmin, params) / X (nmin, params);	/* Since psi(nmin - 1) = 0 */
	  iter_up_start_idx = 2;
	}
      else
	{
	  psi[nminus_idx] = 1.0;
	  iter_up_start_idx = nminus_idx + 1;
	}

      for (i = iter_up_start_idx; i <= nplus_idx; i++)
	{
	  double nn = nmin - 1.0 + i;	/* n - 1 */
	  psi[i] = -(Y (nn, params) * psi[i - 1] +
		      Z (nn, params) * psi[i - 2]) / X (nn, params);
	}

      /* Since we choose nc=nplus, Psi_plus(nc)=1, and we multiply
         Psi_minus(nmin...nplus) by Psi_plus(nc)/Psi_minus(nc) ==
         1/Psi_minus(n_plus) to give us Psi_plus(nmin...nplus). */
      a = 1.0 / psi[nplus_idx];

      for (i = 0; i <= nplus_idx; i++)
	psi[i] *= a;

      normalize (psi, nmin, nmax_idx, params);
      return __SUCCESS;
    }

  if (iter_down)		/* Iterate downwards from nplus, chosing nc = nminus. */
    {
      double a;
      int iter_down_start_idx;

      /* Note that this initialization stuff could be done inside the logic of
         iterating LL98 Eq. 2 above. However following that design leads to some
         rather obscure corner cases and errors, so it's cleaner to do it
         here. Really, don't move it. */
      if (nplus_idx > nmax_idx - 2)
	{
	  psi[nplus_idx] = 1.0;
	  psi[nplus_idx - 1] = -Y (nmax, params) / Z (nmax, params);
	  iter_down_start_idx = nplus_idx - 2;
	}
      else
	{
	  psi[nplus_idx] = 1.0;
	  iter_down_start_idx = nplus_idx - 1;
	}

      for (i = iter_down_start_idx; i >= nminus_idx; i--)
	{
	  double nn = nmin + 1.0 + i;	/* n + 1 */
	  psi[i] = -(X (nn, params) * psi[i + 2] +
		      Y (nn, params) * psi[i + 1]) / Z (nn, params);
	}

      /* Since we choose nc=nminus, Psi_minus(nc)=1, and we multiply
         Psi_plus(nminus...nmax) by Psi_minus(nc)/Psi_plus(nc) ==
         1/Psi_plus(n_plus) to give us Psi_minus(nminus...nmax). */
      a = 1.0 / psi[nminus_idx];

      for (i = nmax_idx; i >= nminus_idx; i--)
	psi[i] *= a;

      normalize (psi, nmin, nmax_idx, params);
      return __SUCCESS;
    }

  fprintf (stderr, "LL98: Could not iterate in either direction\n");
  return __FAILURE;

}

/* Specific functions for calculation of 3j coefficients using j recurrsion -
   first column of Table 1 of LL98. */
typedef struct params_3j_j
{
  int two_j2, two_j3, two_m2, two_m3, two_jmin, two_jmax;
  double j2, j3, m2, m3;
} params_3j_j;

static inline double
A (const double j, const double j2, const double j3,
   const double m2, const double m3)
{
  double a = j * j;
  double b = j2 - j3;
  double c = j2 + j3 + 1.0;
  double d = m2 + m3;

  return sqrt ((a - b * b) * (c * c - a) * (a - d * d));
}

static inline double
B (const double j, const double j2, const double j3,
   const double m2, const double m3)
{
  double a = m2 + m3;
  double b = j2 * (j2 + 1.0) - j3 * (j3 + 1.0);
  double c = (m2 - m3) * j * (j + 1.0);

  return (2.0 * j + 1.0) * (a * b - c);
}

static inline double
X_3j_j (const double j, const void *params)
{
  params_3j_j *p = (params_3j_j *) params;
  return j * A (j + 1.0, p->j2, p->j3, p->m2, p->m3);
}

static inline double
Y_3j_j (const double j, const void *params)
{
  params_3j_j *p = (params_3j_j *) params;
  return B (j, p->j2, p->j3, p->m2, p->m3);
}

static inline double
Z_3j_j (const double j, const void *params)
{
  params_3j_j *p = (params_3j_j *) params;
  return (j + 1.0) * A (j, p->j2, p->j3, p->m2, p->m3);
}

static void
normalize_3j_j (double *f, const double jmin, const int jmax_idx,
		const void *params)
{
  params_3j_j *p = (params_3j_j *) params;
  double a = 0.0;
  int i;

  for (i = 0; i <= jmax_idx; i++)
    {
      double j = jmin + i;
      double ff = f[i];
      a += (2.0 * j + 1.0) * ff * ff;
    }

  a = 1.0 / sqrt (a); /* Normalization factor */

  /* Change sign of a to give correct sign of f(j_max) - see last row
     of table 1 in LL98 */
  if (__ODD ((p->two_j2 - p->two_j3 + p->two_m2 + p->two_m3) / 2))
    {
      if (f[jmax_idx] > 0)
        a = -a;
    }
  else
    {
      if (f[jmax_idx] < 0)
        a = -a;
    }

  for (i = 0; i <= jmax_idx; i++)
    f[i] *= a;
}

static inline double
single_val_3j_j (const void *params)
{
  params_3j_j *p = (params_3j_j *) params;
  double a = 1.0 / sqrt (p->two_jmin + 1.0);

  if (__ODD ((p->two_j2 - p->two_j3 + p->two_m2 + p->two_m3) / 2))
    return -a;
  else
    return a;
}

int
wigner3j_family_j (const int two_j2, const int two_j3,
		   const int two_m2, const int two_m3,
		   double family[], int *two_jmin, int *two_jmax)
{
  params_3j_j p;
  int a, b;

  if ((abs (two_m2) > two_j2) || (abs (two_m3) > two_j3))
    return __FAILURE;

  /* Check (j,m) pairs are both integer or both half integer */
  if (__ODD (two_m2 + two_j2) || __ODD (two_m3 + two_j3))
    return __FAILURE;

  a = abs (two_j2 - two_j3);
  b = abs (two_m2 + two_m3);

  *two_jmin = __MAX(a, b);
  *two_jmax = two_j2 + two_j3;

  p.two_j2 = two_j2;
  p.two_j3 = two_j3;
  p.two_m2 = two_m2;
  p.two_m3 = two_m3;

  p.j2 = two_j2 / 2.0;
  p.j3 = two_j3 / 2.0;
  p.m2 = two_m2 / 2.0;
  p.m3 = two_m3 / 2.0;

  p.two_jmin = *two_jmin;
  p.two_jmax = *two_jmax;

  LL98 (family, *two_jmin, *two_jmax, &p, X_3j_j, Y_3j_j, Z_3j_j,
	normalize_3j_j, single_val_3j_j);

  return __SUCCESS;
}

static inline bool
is_triangle (const int two_ja, const int two_jb, const int two_jc)
/* Returns 0 if the triangle condition is met, and the sum of twice
   the angular momenta is even. Arguments are twice the value of the
   angular momenta. */
{
  if ((two_jc <= two_ja + two_jb) && (two_jc >= abs (two_ja - two_jb))
      && (__EVEN (two_ja + two_jb + two_jc)))
    return true;
  else
    return false;
}

double
wigner3j (const int two_j1, const int two_j2, const int two_j3,
          const int two_m1, const int two_m2, const int two_m3)
/* Returns a single Wigner 3J symbol. Computation is done by
   calculating the family of 3J symbols by J recursion with two_j2 and
   two_j3 fixed, and returning the symbol corresponding to the
   specified value of two_j1. In most cases this will require less
   computation than using M recursion. */
{
  int a = abs (two_j2 - two_j3);
  int b = abs (two_m2 + two_m3);
  int two_jmin = __MAX(a, b);
  int two_jmax = two_j2 + two_j3;
  int dim = 1 + (two_jmax - two_jmin) / 2;
  double buff[__MAX (dim, 1)]; /* Don't want 0 or negative size here. */

  /* Check that |m| <= j */
  if ((abs (two_m1) > two_j1) ||
      (abs (two_m2) > two_j2) ||
      (abs (two_m3) > two_j3))
    return NAN;

  /* Check (j,m) pairs are both integer or both half integer */
  if (__ODD (two_m2 + two_j2) ||
      __ODD (two_m2 + two_j2) ||
      __ODD (two_m3 + two_j3))
    return NAN;

  /* Check for positive j */
  if (two_j1 < 0 || two_j2 < 0 || two_j3 < 0)
    return NAN;

  /* And, finally check triangle condition is met */
  if (!is_triangle (two_j1, two_j2, two_j3))
    return 0.0;

  wigner3j_family_j (two_j2, two_j3, two_m2, two_m3,
                     buff, &two_jmin, &two_jmax);

  return buff[(two_j1 - two_jmin) / 2];
}

/* End of specifics for 3j calculation by j recurrsion. */


/* Specific functions for calculation of 3j coefficients using m recurrsion -
   second column of Table 1 of LL98. */
typedef struct params_3j_m
{
  int two_j1, two_j2, two_j3, two_m1, two_mmin, two_mmax;
  double j1, j2, j3, m1;
} params_3j_m;

static inline double
C (const double m, const double j1, const double j2, const double j3,
   const double m1)
{
  return sqrt ((j2 - m + 1) * (j2 + m) * (j3 - m - m1 + 1.0) * (j3 + m + m1));
}

static inline double
D (const double m, const double j1, const double j2, const double j3,
   const double m1)
{

  return j2 * (j2 + 1.0) + j3 * (j3 + 1.0) - j1 * (j1 + 1.0) - 2.0 * m * (m +
									  m1);
}

static inline double
X_3j_m (const double m, const void *params)
{
  params_3j_m *p = (params_3j_m *) params;
  return C (m + 1.0, p->j1, p->j2, p->j3, p->m1);
}

static inline double
Y_3j_m (const double m, const void *params)
{
  params_3j_m *p = (params_3j_m *) params;
  return D (m, p->j1, p->j2, p->j3, p->m1);
}

static inline double
Z_3j_m (const double m, const void *params)
{
  params_3j_m *p = (params_3j_m *) params;
  return C (m, p->j1, p->j2, p->j3, p->m1);
}

static inline void
normalize_3j_m (double *f, const double mmin, const int mmax_idx,
		const void *params)
{
  params_3j_m *p = (params_3j_m *) params;
  double a = 0.0, phase;
  int i;

  for (i = 0; i <= mmax_idx; i++)
    {
      double ff = f[i];
      a += ff * ff;
    }

  a *= 2.0 * p->j1 + 1.0;
  a = 1.0 / sqrt (a);

  if (__ODD ((p->two_j2 - p->two_j3 - p->two_m1) / 2))
    phase = -1.0;
  else
    phase = 1.0;

  if ((f[mmax_idx] / phase) < 0)
    a = -a;

  for (i = 0; i <= mmax_idx; i++)
    f[i] *= a;
}

static inline double
single_val_3j_m (const void *params)
{
  params_3j_m *p = (params_3j_m *) params;
  double a = 1.0 / sqrt (p->two_j1 + 1.0);

  if (__ODD ((p->two_j1 - p->two_m1) / 2))
    return -a;
  else
    return a;
}

int
wigner3j_family_m (const int two_j1, const int two_j2, const int two_j3,
		   const int two_m1, double family[], int *two_mmin,
		   int *two_mmax)
{
  params_3j_m p;

  if (!is_triangle (two_j1, two_j2, two_j3))
    return __FAILURE;

  if (abs (two_m1) > two_j1)
    return __FAILURE;

  /* Check j1 and m1 are both integer or both half integer */
  if (__ODD (two_m1 + two_j1))
    return __FAILURE;

  *two_mmin = __MAX(-two_j2, -two_j3 - two_m1);
  *two_mmax = __MIN(two_j2, two_j3 - two_m1);

  p.two_j1 = two_j1;
  p.two_j2 = two_j2;
  p.two_j3 = two_j3;
  p.two_m1 = two_m1;

  p.j1 = two_j1 / 2.0;
  p.j2 = two_j2 / 2.0;
  p.j3 = two_j3 / 2.0;
  p.m1 = two_m1 / 2.0;

  p.two_mmin = *two_mmin;
  p.two_mmax = *two_mmax;

  LL98 (family, *two_mmin, *two_mmax, &p, X_3j_m, Y_3j_m, Z_3j_m,
	normalize_3j_m, single_val_3j_m);

  return __SUCCESS;
}


/* End of specifics for 3j calculation by m recurrsion. */



#undef __ODD
#undef __EVEN
#undef __MAX
#undef __MIN
#undef __SMALL
#undef __SUCCESS
#undef __FAILURE
