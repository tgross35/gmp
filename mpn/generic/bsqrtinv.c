/* mpn_bsqrtinv, compute r such that r^2 * y = 1 (mod 2^{b+1}).

   Contributed to the GNU project by Martin Boij (as part of perfpow.c),
   optimized by Marco Bodrato.

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
#include "longlong.h"

/* Compute r such that r^2 * y = 1 (mod 2^{b+1}).
   Return non-zero if such an integer r exists.

   Iterates
     t  <-- (r^2 y - 1) / 2
     r' <-- r - r t , (or its negation, sometimes)
   using Hensel lifting.  Since we divide by two, the Hensel lifting is
   somewhat degenerates.  Therefore, we lift from 2^b to 2^{b+1}-1.

   Somewhere, the Halley recursion is used, lifting from n to 3n-1.
     t  <-- (r^2 y - 1) / 2
     r' <-- r - r t + r 3t^2 / 2

   FIXME:
     (1) Simplify to do precision book-keeping in limbs rather than bits.

     (2) Take advantage of zero low part of r^2 y - 1.

     (3) Use wrap-around trick.
*/

#ifndef BSQRTINV_DONT_USE_TABLE
/* Generated with GP-Pari:
   b=8;v=Vecsmall(-binary(2^2^b-1));
   forstep(i=1,2^(b+1),2,v[lift(Mod(i,2^(b+3))^-2)>>3+1]=i\2);
   for(i=1,2^b,print1(v[i],",");if(i%16==0,print(),print1(" ")))
 */
static const unsigned char binvsqrttab[256] = /* The least significant 1 was removed */
  {  0, 170, 172, 102, 184, 253, 219, 129, 240, 218, 227,  22, 168,  50, 148,  46,
    31, 245, 115,  57, 152, 157,   4, 222, 208, 197,   3, 137, 136, 146, 139, 113,
    63, 149, 108, 217, 120,  61, 228,  62, 176, 101, 220, 214, 104, 242,  84, 238,
    95,  53, 179, 134,  88,  34,  59,  97, 144,   5,  67,  54,  72, 173, 203,  78,
   127,  42,  44,  25,  56, 130, 164, 254, 112,  90, 156, 105,  40,  77,  20,  81,
   159, 138, 243, 185,  24, 226, 123,  94,  80, 186, 131, 246,   8,  18, 244, 241,
   191, 234,  19, 166,   7, 189, 100,  65,  48, 229,  92,  86,  23, 114,  43, 110,
   223, 181, 204,   6,  39,  93, 187, 225,  16, 133, 195,  73,  55, 210, 180,  49,
   255,  85,  83, 153,  71,   2,  36, 126,  15,  37,  28, 233,  87, 205, 107, 209,
   224,  10, 140, 198, 103,  98, 251,  33,  47,  58, 252, 118, 119, 109, 116, 142,
   192, 106, 147,  38, 135, 194,  27, 193,  79, 154,  35,  41, 151,  13, 171,  17,
   160, 202,  76, 121, 167, 221, 196, 158, 111, 250, 188, 201, 183,  82,  52, 177,
   128, 213, 211, 230, 199, 125,  91,   1, 143, 165,  99, 150, 215, 178, 235, 174,
    96, 117,  12,  70, 231,  29, 132, 161, 175,  69, 124,   9, 247, 237,  11,  14,
    64,  21, 236,  89, 248,  66, 155, 190, 207,  26, 163, 169, 232, 141, 212, 145,
    32,  74,  51, 249, 216, 162,  68,  30, 239, 122,  60, 182, 200,  45,  75, 206};
#endif

/* tp needs 2*(1 + bnb / GMP_NUMB_BITS) limbs of space */
int
mpn_bsqrtinv (mp_ptr rp, mp_srcptr yp, mp_bitcnt_t bnb, mp_ptr tp)
{
  mp_limb_t y0 = *yp;
  ASSERT (bnb > 0);

#ifndef BSQRTINV_RP_NOT_ZEROED
  ASSERT (mpn_zero_p (rp, 1 + bnb / GMP_NUMB_BITS));
#endif
  if (UNLIKELY (bnb == 1))
    {
      *rp = y0;
      return (y0 & 3) == 1;
    }
  else
    {
      mp_ptr tp2 = tp + 1 + bnb / GMP_NUMB_BITS;
      mp_size_t bn, order[GMP_LIMB_BITS + 1];
      mp_limb_t t0, r0;
      int i;

      if ((y0 & 7) != 1)
	return 0;

      /* Rekte : 4 */
      /* 4 -4x-> 11 -4x-> 32 (ne 33, 8x) */
      /* 4 -4x-> 11 -3x-> 21 -3x-> preter 33 (10x) [4,33]*/

      /* 4 -4x-> 11 -4x-> 32 -3x-> 63 (ne 64, 11x) [3,4]*/
      /* 4 -4x-> 11 -4x-> 32 -4x-> preter 65 (12x) */

      /* 4 -3x-> 7 -3x-> 13 -3x-> 25 -4x-> preter 65 (13x) */
      /* 4 -3x-> 7 -3x-> 13 -3x-> 25 -3x-> preter 33 (12x) */
      /* 4 -3x-> 7 -3x-> 13 -4x-> preter 33 (10x) */

      /* Tabulo: 10 */
      /* 10 -4x-> 29 -4x-> preter 65, 8x [4,4]*/
      /* 10 -3x-> 19 -3x-> preter 33, 6x [3,3]*/

#ifdef BSQRTINV_DONT_USE_TABLE
      /* 16-bits computations are enough */
      unsigned ru = y0 + ((y0 & 8) >> 2) + ((y0 & 16) >> 1);

      unsigned tu = ru * ru * (unsigned) y0 >> 1;
      ASSERT ((tu & (GMP_NUMB_MAX >> (GMP_NUMB_BITS - 4))) == 0);
      ru += ru * tu * ((tu >> 1) + tu - 1); /* Halley -> 11 */
      /* Better sequences are possible from size 11, but this
	 code is currently not used. */
#else /* ! defined(BSQRTINV_DONT_USE_TABLE) */
      unsigned ru = binvsqrttab[(y0 >> 3) & 0xff];
      ru = (ru << 1) + 1;
#endif
      r0 = ru;

#if GMP_NUMB_BITS < 10 * 2 - 2
      const mp_bitcnt_t precomputed_bits = 10;
#else /* GMP_NUMB_BITS >= 10 * 2 - 2 */
      t0 = r0 * r0 * y0 >> 1;
      ASSERT ((t0 & (GMP_NUMB_MAX >> (GMP_NUMB_BITS - 10))) == 0);
#if GMP_NUMB_BITS < 19 * 2 - 2
      r0 -= r0 * t0;
      const mp_bitcnt_t precomputed_bits = 19;
#else /* GMP_NUMB_BITS >= 19 * 2 - 2 */
      r0 += r0 * t0 * ((t0 >> 1) + t0 - 1); /* Halley -> 29*/
#if GMP_NUMB_BITS < 29 * 3 - 2
      const mp_bitcnt_t precomputed_bits = 29;
#else /* GMP_NUMB_BITS >= 29 * 3 - 2 */
#error Not implemented, yet.
#endif
#endif
#endif

      i = 0;
      for (; bnb > GMP_NUMB_BITS + 1; bnb = (bnb + 2) >> 1)
	order[i++] = bnb;
      if (bnb > precomputed_bits) {
	if (bnb >= GMP_NUMB_BITS) {
	  mp_limb_t r0h = r0 >> 1;
	  /* We could gain the third bit with (r0h|1)*((r0h+1)>>1) */
	  mp_limb_t r0sqm1 = r0h * (r0h + 1); /* r0*r0 >> 2 */
	  mp_limb_t yh = (y0 >> 2) + (yp[1] << (GMP_NUMB_BITS - 2));
	  mp_limb_t yrrm1d4 = y0 * r0sqm1 + yh; /* (r0*r0*y0-1) >> 2 */
	  ASSERT ((yrrm1d4 & (GMP_NUMB_MAX >> (GMP_NUMB_BITS - precomputed_bits + 1))) == 0);
#if GMP_NUMB_BITS < 19 * 2 - 2
	  mp_limb_t rt = r0h - r0 * yrrm1d4;  /* (r - r*(r0*r0*y0-1)/2) >>1 */
#else /* GMP_NUMB_BITS >= 19 * 2 - 2 */
	  mp_limb_t rt = r0h + r0 * yrrm1d4 * (yrrm1d4 * 3 - 1); /* Halley, x3 - 1 */
#endif
	  r0 = (rt << 1) ^ ((rt & GMP_LIMB_HIGHBIT) ? GMP_NUMB_MAX : CNST_LIMB(1));
#ifdef BSQRTINV_RP_NOT_ZEROED
	  rp[1] = 0;
#endif
	} else {
	  t0 = r0 * r0 * y0 >> 1;
#if GMP_NUMB_BITS < 19 * 2 - 2
	  r0 -= r0 * t0;
#else /* GMP_NUMB_BITS >= 19 * 2 - 2 */
	  r0 += r0 * t0 * ((t0 >> 1) + t0 - 1); /* Halley, x3 - 1 */
#endif
	  ASSERT ((t0 & (GMP_NUMB_MAX >> (GMP_NUMB_BITS - precomputed_bits))) == 0);
	}
      }

      if (i) {
	mp_limb_t t4, t3, t2, t1, r1;

	umul_ppmm (t1, t0, r0, r0); /* [t1,t0] <- r^2 */
	if (bnb <= GMP_NUMB_BITS) {
	  umul_ppmm (t3, t2, y0, t0);
	  t3 += y0 * t1 + yp[1] * t0;
	  t2 = ((t2 >> 1) | (t3 << (GMP_NUMB_BITS - 1))) & GMP_NUMB_MAX;
	  t3 = t3 >> 1 ; /* [t3,t2] <- (r^2 y - 1) / 2 */

	  /* [r1,t4] <- r (r^2 y - 1) / 2 */
	  umul_ppmm (r1, t4, r0, t2);
	  r1 += r0 * t3;

	  /* r (r^2 y - 1) / 2 - r */
	  sub_ddmmss(rp[1], rp[0], r1, t4, 0, r0);
	} else {
	  t0 = (t0 >> 2) | (t1 << (GMP_NUMB_BITS -2));
	  t1 = (t1 >> 2); /* [t1,t0] = r0*r0 >> 2 */
	  umul_ppmm (t3, t2, y0, t0);
	  t3 += y0 * t1 + yp[1] * t0;
	  t3 += (yp[1] >> 2) + (yp[2] << (GMP_NUMB_BITS - 2)) + (t2 != 0);
	  ASSERT (t2 + (y0 >> 2) + (yp[1] << (GMP_NUMB_BITS - 2)) == 0);
	  /* [t3,0] <- (r0^2 y - 1) / 2 / 2 */

	  t4 = t3 * r0 - 1;
	  /* [2*t4+1,-r0] <- r0*(r0^2 y-1)/2 - r0 */
	  if (t4 & GMP_LIMB_HIGHBIT) {
	    rp[0] = r0 & GMP_NUMB_MAX;
	    rp[1] = ~t4 << 1;
	  } else {
	    rp[0] = -r0 & GMP_NUMB_MAX;
	    rp[1] = (t4 << 1) ^ 1;
	  }
#ifdef BSQRTINV_RP_NOT_ZEROED
	  rp[2] = 0;
#endif
	}
	--i;

	for (bn = 2 + (bnb > GMP_NUMB_BITS); --i >= 0;)
	  {
	    mp_size_t pbn = bn;
	    /* sqr may partially overwrite tp2, but that part is unused here */
	    mpn_sqr (tp, rp, bn); /* tp <- r^2 */

	    bnb = order[i];
	    bn = 1 + bnb / GMP_LIMB_BITS;

#ifdef  BSQRTINV_USE_MULMID
	    if (pbn > 20) /* Should be tuned, if it makes sense */
	      {
		int offset = 2; /* The larger, the safer */
		mpn_mullo_n (tp2 + pbn - offset, yp, tp + pbn - offset, bn - pbn + offset);
		/* mulmid partially overwrites tp2, but that part is unused here */
		mpn_mulmid (tp + pbn - offset, yp, bn, tp, pbn - offset);
		ASSERT (tp [bn + 2] < GMP_NUMB_MAX); /* MPN_INCR will be stopped here */
		MPN_INCR_U (tp + pbn, bn - pbn + 3, ! mpn_zero_p (tp + pbn - offset, offset) |
			    ! mpn_zero_p (tp2 + pbn - offset, offset - 1));
		mpn_add_nc (tp2 + pbn - 1, tp2 + pbn - 1, tp + pbn,
			    bn - pbn + 1, (tp2[pbn - 1] ^ tp[pbn]) & 1);
	      }
	    else
#endif
	      {
		mpn_mullo_n (tp2, yp, tp, bn); /* tp2 <- rp^2 y */
		ASSERT (tp2[0] == CNST_LIMB (1));
		ASSERT (pbn == 2 || mpn_zero_p (tp2 + 1, pbn - 2));
		/* tp2 <- (rp^2 y - 1) / 2 (skip the lowest limbs) */
	      }
	    ASSERT_NOCARRY (mpn_rshift (tp2 + pbn - 1, tp2 + pbn - 1, bn - pbn + 1, 1));

	    /* tp <- r (r^2 y - 1) / 2 (only the relevant limbs) */
#ifdef BSQRTINV_RP_NOT_ZEROED
	    rp [pbn] = 0;
#endif
	    ASSERT (pbn >= bn - pbn + 1 || (pbn == bn - pbn && rp [pbn] == 0));
	    mpn_mullo_n (tp, rp, tp2 + pbn - 1, bn - pbn + 1);

	    /* mpn_rsub_1 (rp + pbn - 1, tp, bn - pbn + 1, rp[pbn - 1]) */
	    int borrow;
	    SUBC_LIMB (borrow, rp[pbn - 1], rp[pbn - 1], tp[0]);

	    if (borrow)
	      mpn_com (rp + pbn, tp + 1, bn - pbn);
	    else
	      mpn_neg (rp + pbn, tp + 1, bn - pbn);
	    /* rp <- r - r (r^2 y - 1) / 2 */
	  }
      } else {
	*rp = r0 & GMP_NUMB_MAX;
      }
    }
  return 1;
}
