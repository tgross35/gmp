/* mpn_bsqrt, a^{1/2} (mod 2^n), for odd a.

Copyright 2009, 2010, 2012, 2015, 2026 Free Software Foundation, Inc.

This file is part of the GNU MP Library.

The GNU MP Library is free software; you can redistribute it and/or modify
it under the terms of either:

  * the GNU Lesser General Public License as published by the Free
    Software Foundation; either version 3 of the License, or (at your
    option) any later version.

or

  * the GNU General Public License as published by the Free Software
    Foundation; either version 2 of the License, or (at your option) any
    later version.

or both in parallel, as here.

The GNU MP Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received copies of the GNU General Public License and the
GNU Lesser General Public License along with the GNU MP Library.  If not,
see https://www.gnu.org/licenses/.  */

#include "gmp-impl.h"

/* tp needs (1 + bnb / GMP_NUMB_BITS) limbs of space + the scratch
   used by mpn_bsqrtinv, i.e 3*(1 + bnb / GMP_NUMB_BITS) */
int
mpn_bsqrt (mp_ptr rp, mp_srcptr ap, mp_bitcnt_t nb, mp_ptr tp)
{
  mp_ptr sp;
  mp_size_t n;

  ASSERT (nb > 0);

  n = 1 + nb / GMP_NUMB_BITS;
  sp = tp + n;

  MPN_FILL (tp, n, CNST_LIMB(0));
  if (! mpn_bsqrtinv (tp, ap, nb, sp))
    return 0;

  if (n == 1)
    *rp = *tp * *ap;
  else
    mpn_mullo_n (rp, tp, ap, n);

  return 1;
}
