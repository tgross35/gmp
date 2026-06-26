/* mpz_init_setbit(integer, val) -- Initialize and assign INTEGER with 2^VAL.

Copyright 2024 Free Software Foundation, Inc.

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

void
mpz_init_setbit (mpz_ptr dest, mp_bitcnt_t bit_idx)
{
  mp_ptr dp;
  mp_size_t limb_idx = bit_idx / GMP_NUMB_BITS;
  mp_limb_t mask = CNST_LIMB(1) << (bit_idx % GMP_NUMB_BITS);

  ALLOC (dest) = SIZ (dest) = limb_idx + 1;
  PTR (dest) = dp = __GMP_ALLOCATE_FUNC_LIMBS (limb_idx + 1);
  MPN_ZERO (dp, limb_idx);
  dp[limb_idx] = mask;
}
