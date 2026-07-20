/* Copyright 2012, 2015, 2026 Free Software Foundation, Inc.

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


#include <stdlib.h>		/* for abort */
#include <stdio.h>		/* for printf */

#include "gmp-impl.h"
#include "tests/tests.h"

#define MAX_LIMBS 150
#define COUNT 500

int
main (int argc, char **argv)
{
  gmp_randstate_ptr rands;

  mp_ptr ap, rp, pp, sp;
  mp_limb_t before_rp, before_sp;
  int count = COUNT;
  unsigned i;
  TMP_DECL;

  TMP_MARK;

  TESTS_REPS (count, argv, argc);

  tests_start ();
  rands = RANDS;

  ap = TMP_ALLOC_LIMBS (MAX_LIMBS);
  rp = TMP_ALLOC_LIMBS (MAX_LIMBS + 2) + 1;
  pp = TMP_ALLOC_LIMBS (MAX_LIMBS);
  sp = TMP_ALLOC_LIMBS (3*MAX_LIMBS + 2) + 1;

  before_rp = rp [-1] = gmp_urandomm_ui (rands, GMP_NUMB_MAX);
  before_sp = sp [-1] = gmp_urandomm_ui (rands, GMP_NUMB_MAX);
  for (i = 0; i < count; i++)
    {
      mp_size_t n;
      mp_bitcnt_t bn;
      mp_limb_t after_rp, after_sp;
      int res;

      n = 1 + gmp_urandomm_ui (rands, MAX_LIMBS);

      if (i & 1)
	mpn_random2 (ap, n);
      else
	mpn_random (ap, n);

      if (i & 0xf)
	ap[0] = (ap[0] | 7) ^ 6;

      bn = 1 + gmp_urandomm_ui (rands, GMP_NUMB_BITS - 2*(n == 1));

      after_rp = rp [n - (bn >= GMP_NUMB_BITS - 1)] = gmp_urandomm_ui (rands, GMP_NUMB_MAX);
      after_sp = sp [3 * n] = gmp_urandomm_ui (rands, GMP_NUMB_MAX);
      res = mpn_bsqrt (rp, ap, n * GMP_NUMB_BITS - bn, sp);
      if (rp [n - (bn >= GMP_NUMB_BITS - 1)] != after_rp || rp [-1] != before_rp ||
	  sp [3 * n] != after_sp || sp [-1] != before_sp)
	{
	  gmp_fprintf (stderr,
		       "mpn_bsqrt memoty bounds violated: %u limbs, - %u bits, res %i [%i]\n",
		       (unsigned) n, (unsigned) bn, res, i);
	  gmp_fprintf (stderr, "before_rp: %Mx <> %Mx, \t", before_rp, rp[-1]);
	  gmp_fprintf (stderr, " after_rp: %Mx <> %Mx\n", after_rp, rp[n - (bn >= GMP_NUMB_BITS - 1)]);
	  gmp_fprintf (stderr, "before_sp: %Mx <> %Mx, \t", before_sp, sp[-1]);
	  gmp_fprintf (stderr, " after_sp: %Mx <> %Mx\n", after_sp, sp[3 * n]);
	  gmp_fprintf (stderr, "a     = %Nx\n", ap, n);
	  abort ();
	}

      if (!res && ((*ap & (7 >> ((n == 1) && (bn == GMP_NUMB_BITS - 1)))) != 1))
	continue;

      mpn_sqrlo (pp, rp, n);

      if (!res || ((n!=1) && (mpn_cmp (pp, ap, n - 1) != 0)) ||
	  ((bn != GMP_NUMB_BITS) && ((pp[n - 1] ^ ap[n - 1]) & GMP_NUMB_MAX >> bn != 0)))
	{
	  gmp_fprintf (stderr,
		       "mpn_bsqrt returned bad result: %u limbs, - %u bits, res %i [%i]\n",
		       (unsigned) n, (unsigned) bn, res, i);
	  gmp_fprintf (stderr, "a     = %Nx\n", ap, n);
	  gmp_fprintf (stderr, "r     = %Nx\n", rp, n);
	  gmp_fprintf (stderr, "r^2   = %Nx\n", pp, n);
	  abort ();
	}
    }
  TMP_FREE;
  tests_end ();
  return 0;
}
