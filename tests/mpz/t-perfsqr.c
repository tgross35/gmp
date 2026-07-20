/* Test mpz_perfect_square_p and mpz_perfect_square_root.

Copyright 2000-2002, 2026 Free Software Foundation, Inc.

This file is part of the GNU MP Library test suite.

The GNU MP Library test suite is free software; you can redistribute it
and/or modify it under the terms of the GNU General Public License as
published by the Free Software Foundation; either version 3 of the License,
or (at your option) any later version.

The GNU MP Library test suite is distributed in the hope that it will be
useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
Public License for more details.

You should have received a copy of the GNU General Public License along with
the GNU MP Library test suite.  If not, see https://www.gnu.org/licenses/.  */

#include <stdio.h>
#include <stdlib.h>

#include "gmp-impl.h"
#include "tests.h"

#include "mpn/perfsqr.h"


/* check_modulo() exercises mpz_perfect_square_p/mpz_perfect_square_root on
   squares which cover each possible quadratic residue to each divisor used
   within mpn_perfect_square_p, ensuring those residues aren't incorrectly
   claimed to be non-residues.

   Each divisor is taken separately.  It's arranged that n is congruent to 0
   modulo the other divisors, 0 of course being a quadratic residue to any
   modulus.

   The values "(j*others)^2" cover all quadratic residues mod divisor[i],
   but in no particular order.  j is run from 1<=j<=divisor[i] so that zero
   is excluded.  A literal n==0 doesn't reach the residue tests.  */

void
check_modulo (void)
{
  static const unsigned long  divisor[] = PERFSQR_DIVISORS;
  unsigned long  i, j;

  mpz_t  alldiv, others, n, root, rop;

  mpz_init (alldiv);
  mpz_init (others);
  mpz_init (root);
  mpz_init (n);
  mpz_init (rop);

  /* product of all divisors */
  mpz_set_ui (alldiv, 1L);
  for (i = 0; i < numberof (divisor); i++)
    mpz_mul_ui (alldiv, alldiv, divisor[i]);

  for (i = 0; i < numberof (divisor); i++)
    {
      /* product of all divisors except i */
      mpz_set_ui (others, 1L);
      for (j = 0; j < numberof (divisor); j++)
        if (i != j)
          mpz_mul_ui (others, others, divisor[j]);

      for (j = 1; j <= divisor[i]; j++)
        {
          /* square */
          mpz_mul_ui (root, others, j);
          mpz_mul (n, root, root);
          if (! mpz_perfect_square_p (n))
            {
              printf ("mpz_perfect_square_p got 0, want 1\n");
              mpz_trace ("  n", n);
              abort ();
            }
          if (! mpz_perfect_square_root (rop, n))
            {
              printf ("mpz_perfect_square_root got 0, want 1\n");
              mpz_trace ("  n", n);
              abort ();
            }
          if (! mpz_cmp (rop, n))
            {
              gmp_printf ("mpz_perfect_square_root rop %Zd, want %Zd\n", rop, n);
              abort ();
            }
        }
    }

  mpz_clear (alldiv);
  mpz_clear (others);
  mpz_clear (root);
  mpz_clear (n);
  mpz_clear (rop);
}

/* Check negative, 0, and 1 */
void
check_edge_cases (void)
{
  mpz_t  n, root;
  mpz_init (n);
  mpz_init (root);

  mpz_set_si(n, -4);
  if (mpz_perfect_square_p (n))
    {
      printf ("mpz_perfect_square_p -4 not a square\n");
      abort ();
    }
  if (mpz_perfect_square_root (root, n))
    {
      printf ("mpz_perfect_square_root -4 not a square\n");
      abort ();
    }
  MPZ_CHECK_FORMAT (root);

  // 0 and 1 are both perfect squares with respective root 0 and 1.
  for (unsigned int m = 0; m < 2; m++)
    {
      mpz_set_ui(n, m);
      if (! mpz_perfect_square_p (n))
        {
          printf ("mpz_perfect_square_p %u should be a square\n", m);
          abort ();
        }
      if (! mpz_perfect_square_root (root, n) || mpz_cmp_ui (root, m) != 0)
        {
          printf ("mpz_perfect_square_root %u should have root = %u\n", m, m);
          abort ();
        }
    }

  mpz_clear (root);
  mpz_clear (n);
}

/* Exercise mpz_perfect_square_p compared to what mpz_sqrt says. */
void
check_sqrt (int reps)
{
  mpz_t x2, x2t, x, rop;
  mp_size_t x2n;
  int res;
  int i;
  int want;
  int cnt = 0;
  gmp_randstate_ptr rands = RANDS;
  mpz_t bs;

  mpz_init (bs);

  mpz_init (x2);
  mpz_init (x);
  mpz_init (x2t);
  mpz_init (rop);

  for (i = 0; i < reps; i++)
    {
      mpz_urandomb (bs, rands, 9);
      x2n = mpz_get_ui (bs);
      mpz_rrandomb (x2, rands, x2n);

      res = mpz_perfect_square_p (x2);
      mpz_sqrt (x, x2);
      mpz_mul (x2t, x, x);
      want = mpz_cmp(x2, x2t) == 0;

      if (res != want)
        {
          printf    ("mpz_perfect_square_p and mpz_sqrt differ\n");
          mpz_trace ("   x  ", x);
          mpz_trace ("   x2 ", x2);
          mpz_trace ("   x2t", x2t);
          printf    ("   mpz_perfect_square_p %d\n", res);
          printf    ("   mpz_sqrt             %d\n", want);
          abort ();
        }

      res = mpz_perfect_square_root (rop, x2);
      if (res != want)
        {
          printf    ("mpz_perfect_square_root and mpz_sqrt differ\n");
          mpz_trace ("   x  ", x);
          mpz_trace ("   x2 ", x2);
          mpz_trace ("   x2t", x2t);
          printf    ("   mpz_perfect_square_root %d\n", res);
          printf    ("   mpz_sqrt                %d\n", want);
          abort ();
        }
      MPZ_CHECK_FORMAT (rop);
      if (res && mpz_cmp(x, rop) != 0)
        {
          printf    ("mpz_perfect_square_root and mpz_sqrt differ\n");
          mpz_trace ("   x  ", x);
          mpz_trace ("   x2 ", x2);
          mpz_trace ("   x2t", x2t);
          mpz_trace ("   mpz_perfect_square_root rop", rop);
          mpz_trace ("   mpz_sqrt", x);
          abort ();
        }
      /* Check that same variable as input and output works */
      mpz_set(rop, x2);
      res = mpz_perfect_square_root (rop, rop);
      MPZ_CHECK_FORMAT (rop);
      if (res != want || (res && mpz_cmp(x, rop) != 0))
        {
          printf    ("mpz_perfect_square_root differ when output == input\n");
          mpz_trace ("   x2 ", x2);
          printf    ("   mpz_perfect_square_root %d\n", res);
          mpz_trace ("   mpz_perfect_square_root rop ", rop);
          abort ();
        }

      cnt += res != 0;
    }

  if (reps > 1000 && cnt == 0) {
    printf("No perfect squares found in %d reps\n", reps);
  }
  /* printf ("%d/%d perfect squares\n", cnt, reps); */

  mpz_clear (bs);
  mpz_clear (x2);
  mpz_clear (x);
  mpz_clear (x2t);
  mpz_clear (rop);
}

/* Exercise mpz_perfect_square_p on large squares. */
void
check_sqares (int reps)
{
  mpz_t x2, x, rop;
  mp_size_t xn;
  int res;
  int i;
  int want;
  mp_bitcnt_t l, h;
  gmp_randstate_ptr rands = RANDS;
  mpz_t bs;

  mpz_init (bs);

  mpz_init (x2);
  mpz_init (x);
  mpz_init (rop);

  for (i = 0; i < reps; i++)
    {
      mpz_urandomb (bs, rands, 13);
      xn = mpz_get_ui (bs) + 100;
      mpz_rrandomb (x, rands, xn);

      mpz_mul (x2, x, x);
      res = mpz_perfect_square_root (rop, x2);
      want = 1;

      if (res != want || mpz_cmp (rop, x) !=0)
        {
          printf    ("mpz_perfect_square_p did not detect a square [%i]\n", i);
          mpz_trace ("  x2  ", x2);
          mpz_trace ("  x   ", x);
          mpz_trace ("  rop ", rop);
	  mpz_sub (rop, rop, x);
          mpz_trace ("rop-x ", rop);
          printf    ("  mpz_perfect_square_p %d\n", res);
          printf    ("  want                 %d\n", want);
          abort ();
        }

      l = mpz_scan1 (x2, 0);
      h = mpz_sizeinbase (x2, 2);

      h = MAX (h, l + 5);
      mpz_combit (x2, l + 3 + gmp_urandomm_ui (rands, h - l - 4));

      res = mpz_perfect_square_root (rop, x2);
      MPZ_CHECK_FORMAT (rop);
      want = 0;

      if (res != want)
        {
          printf    ("mpz_perfect_square_p did not detect a non square [%i]\n", i);
          mpz_trace ("  x2  ", x2);
          mpz_trace ("  x   ", x);
          mpz_trace ("  rop ", rop);
	  mpz_sub (rop, rop, x);
          mpz_trace ("rop-x ", rop);
          printf    ("  mpz_perfect_square_p %d\n", res);
          printf    ("  want                 %d\n", want);
          abort ();
        }

    }

  mpz_clear (bs);
  mpz_clear (x2);
  mpz_clear (x);
  mpz_clear (rop);
}


int
main (int argc, char **argv)
{
  int reps = 200000;

  tests_start ();
  mp_trace_base = -16;

  if (argc == 2)
     reps = atoi (argv[1]);

  check_edge_cases ();
  check_modulo ();
  check_sqrt (reps);
  check_sqares (reps >> 2);

  tests_end ();
  exit (0);
}
