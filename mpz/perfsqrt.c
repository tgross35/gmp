/* mpz_perfect_square_root(rop, op) -- Return non-zero if op is a perfect
   square and sets rop to root, otherwise returns zero.

Copyright 2026 Free Software Foundation, Inc.

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

int
mpz_perfect_square_root (mpz_ptr rop, mpz_srcptr op)
{
  mp_size_t u_size = SIZ (op);

  /* Handle <= 0 */
  if (__GMP_UNLIKELY ( u_size <= 0))
    {
      /* Zero is a square with rop = 0 */
      if ( u_size == 0)
      {
        SIZ(rop) = 0;
        return 1;
      }
      /* Negatives are never perfect squares. */
      return 0;
    }

  if (! mpn_probab_perfect_square_p (PTR (op), u_size))
    return 0;

  /* This is the precise size of sqrt(op) because leading limb is non-zero. */
  mp_size_t rop_size = (SIZ (op) + 1) / 2;

  /* mpn_sqrtrem doesn't allow sp == np */
  if (rop == op)
    {
      TMP_DECL;
      TMP_MARK;

      mp_ptr rop_ptr = TMP_ALLOC_LIMBS (rop_size);
      int res = ! mpn_sqrtrem (rop_ptr, NULL, PTR (op), u_size);

      /* It might be better to always copy so the value in rop doesn't
         depend on if the arguments are aliased or not. */
      if (res)
        {
          MPN_COPY (PTR(rop), rop_ptr, rop_size);
          SIZ(rop) = rop_size;
        }

      TMP_FREE;
      return res;
    }

  /* Make sure rop has exact space for the square root. */
  SIZ(rop) = rop_size;
  return ! mpn_sqrtrem (MPZ_REALLOC (rop, rop_size), NULL, PTR (op), u_size);
}
