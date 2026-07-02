/* Alternate implementations of binvert_limb to compare speeds. */

/*
Copyright 2000, 2002, 2022 Free Software Foundation, Inc.

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
#include "speed.h"

#define binvert_limb_uintfast(inv,n)					\
  do {									\
    mp_limb_t  __n = (n);						\
    mp_limb_t  __inv;							\
    UHWtype  __inv8;							\
    ASSERT ((__n & 1) == 1);						\
									\
    __inv8 = binvert_limb_table[(__n & 0xff)/2]; /*  8 */		\
    if (GMP_NUMB_BITS > 16) {						\
      UHWtype __inv16 = 2 * __inv8 - __inv8 * __inv8 * __n;		\
    if (GMP_NUMB_BITS > 32) {						\
      UHWtype __inv32 = 2 * __inv16 - __inv16 * __inv16 * __n;		\
      __inv = 2*(mp_limb_t)__inv32 - (mp_limb_t)__inv32*__inv32*__n;	\
    } else {\
      __inv = 2*(mp_limb_t)__inv16 - (mp_limb_t)__inv16*__inv16*__n;	\
    };									\
    } else {\
      __inv = 2 * (mp_limb_t)__inv8 - (mp_limb_t)__inv8 * __inv8 * __n;	\
    };									\
									\
    if (GMP_NUMB_BITS > 64)						\
      {									\
	int  __invbits = 64;						\
	do {								\
	  __inv = 2 * __inv - __inv * __inv * __n;			\
	  __invbits *= 2;						\
	} while (__invbits < GMP_NUMB_BITS);				\
      }									\
									\
    ASSERT ((__inv * __n & GMP_NUMB_MASK) == 1);			\
    (inv) = __inv & GMP_NUMB_MASK;					\
  } while (0)

/* A (not exact) copy of the binvert function contained in
   mpn/generic/sec_powm.c .
 */

static mp_limb_t
sec_binvert_limb (mp_limb_t n)
{
  mp_limb_t inv, t;
  ASSERT ((n & 1) == 1);
  /* 3 + 2 -> 5 */
  UHWtype th, invh = n + (((n + 1) << 1) & 0x18);
  /* UHWtype th, invh = (n * 3) ^ 2; */

  th = n * invh;
#if GMP_NUMB_BITS <= 10
  /* 5 x 2 -> 10 */
  inv = 2 * invh - invh * th;
#else /* GMP_NUMB_BITS > 10 */
  /* 5 x 2 + 2 -> 12 */
  invh = 2 * invh - invh * th + ((invh<<10)&-(th&(1<<5)));
#endif /* GMP_NUMB_BITS <= 10 */

  if (GMP_NUMB_BITS > 12)
    {
      t = n * invh - 1;
      if (GMP_NUMB_BITS <= 36)
	{
	  /* 12 x 3 -> 36 */
	  inv = invh + invh * t * (t - 1);
	}
      else /* GMP_NUMB_BITS > 36 */
	{
	  mp_limb_t t2 = t * t;
#if GMP_NUMB_BITS <= 60
	  /* 12 x 5 -> 60 */
	  inv = invh + invh * (t2 + 1) * (t2 - t);
#else /* GMP_NUMB_BITS > 60 */
	  /* 12 x 5 + 4 -> 64 */
	  inv = invh * ((t2 + 1) * (t2 - t) + 1 - ((t<<48)&-(t&(1<<12))));

	  /* 64 -> 128 -> 256 -> ... */
	  for (int todo = (GMP_NUMB_BITS - 1) >> 6; todo != 0; todo >>= 1)
	    inv = 2 * inv - inv * inv * n;
#endif /* GMP_NUMB_BITS <= 60 */
	}
    }

  ASSERT ((inv * n & GMP_NUMB_MASK) == 1);
  return inv & GMP_NUMB_MASK;
}

#define binvert_limb_sec(inv,n) inv = sec_binvert_limb (n)

/* Like the standard version in gmp-impl.h, but with a different path
   for bit sizes larger than 32, with concurrent multiplications.  */

static mp_limb_t
binvert_limb_pipe_f (mp_limb_t n)
{
  mp_limb_t  __inv = 3 * n ^ 2; /* 5 */
  mp_limb_t  __y0 = CNST_LIMB (1) - __inv * n;
  mp_limb_t  __y1 = __y0 * __y0;
  mp_limb_t  __y2 = __y1 * __y1;
  __inv *= (1 + __y0) * (1 + __y1) * (1 + __y2);

  if (GMP_NUMB_BITS > 40)
    {
      int  __invbits = 40;
      do {
	__y2 *= __y2;
	__inv *= __y2 + 1;
	__invbits *= 2;
      } while (__invbits < GMP_NUMB_BITS);
    }

  ASSERT ((__inv * __n & GMP_NUMB_MASK) == 1);
  return __inv & GMP_NUMB_MASK;
}
#define binvert_limb_pipe(inv,n) inv = binvert_limb_pipe_f (n)

/* Like the standard version in gmp-impl.h, but with the expressions using a
   "1-" form.  This has the same number of steps, but "1-" is on the
   dependent chain, whereas the "2*" in the standard version isn't.
   Depending on the CPU this should be the same or a touch slower.  */

#if GMP_LIMB_BITS <= 32
#define binvert_limb_mul1(inv,n)                                \
  do {                                                          \
    mp_limb_t  __n = (n);                                       \
    mp_limb_t  __inv;                                           \
    ASSERT ((__n & 1) == 1);                                    \
    __inv = binvert_limb_table[(__n&0xFF)/2]; /*  8 */          \
    __inv = (1 - __n * __inv) * __inv + __inv;  /* 16 */        \
    __inv = (1 - __n * __inv) * __inv + __inv;  /* 32 */        \
    ASSERT (__inv * __n == 1);                                  \
    (inv) = __inv;                                              \
  } while (0)
#endif

#if GMP_LIMB_BITS > 32 && GMP_LIMB_BITS <= 64
#define binvert_limb_mul1(inv,n)                                \
  do {                                                          \
    mp_limb_t  __n = (n);                                       \
    mp_limb_t  __inv;                                           \
    ASSERT ((__n & 1) == 1);                                    \
    __inv = binvert_limb_table[(__n&0xFF)/2]; /*  8 */          \
    __inv = (1 - __n * __inv) * __inv + __inv;  /* 16 */        \
    __inv = (1 - __n * __inv) * __inv + __inv;  /* 32 */        \
    __inv = (1 - __n * __inv) * __inv + __inv;  /* 64 */        \
    ASSERT (__inv * __n == 1);                                  \
    (inv) = __inv;                                              \
  } while (0)
#endif


/* The loop based version used in GMP 3.0 and earlier.  Usually slower than
   multiplying, due to the number of steps that must be performed.  Much
   slower when the processor has a good multiply.  */

#define binvert_limb_loop(inv,n)                \
  do {                                          \
    mp_limb_t  __v = (n);                       \
    mp_limb_t  __v_orig = __v;                  \
    mp_limb_t  __make_zero = 1;                 \
    mp_limb_t  __two_i = 1;                     \
    mp_limb_t  __v_inv = 0;                     \
                                                \
    ASSERT ((__v & 1) == 1);                    \
                                                \
    do                                          \
      {                                         \
        while ((__two_i & __make_zero) == 0)    \
          __two_i <<= 1, __v <<= 1;             \
        __v_inv += __two_i;                     \
        __make_zero -= __v;                     \
      }                                         \
    while (__make_zero);                        \
                                                \
    ASSERT (__v_orig * __v_inv == 1);           \
    (inv) = __v_inv;                            \
  } while (0)


/* Another loop based version with conditionals, but doing a fixed number of
   steps. */

#define binvert_limb_cond(inv,n)                \
  do {                                          \
    mp_limb_t  __n = (n);                       \
    mp_limb_t  __rem = (1 - __n) >> 1;          \
    mp_limb_t  __inv = GMP_LIMB_HIGHBIT;        \
    int        __count;                         \
                                                \
    ASSERT ((__n & 1) == 1);                    \
                                                \
    __count = GMP_LIMB_BITS-1;               \
    do                                          \
      {                                         \
        __inv >>= 1;                            \
        if (__rem & 1)                          \
          {                                     \
            __inv |= GMP_LIMB_HIGHBIT;          \
            __rem -= __n;                       \
          }                                     \
        __rem >>= 1;                            \
      }                                         \
    while (-- __count);                         \
                                                \
    ASSERT (__inv * __n == 1);                  \
    (inv) = __inv;                              \
  } while (0)


/* Another loop based bitwise version, but purely arithmetic, no
   conditionals. */

#define binvert_limb_arith(inv,n)                                       \
  do {                                                                  \
    mp_limb_t  __n = (n);                                               \
    mp_limb_t  __rem = (1 - __n) >> 1;                                  \
    mp_limb_t  __inv = GMP_LIMB_HIGHBIT;                                \
    mp_limb_t  __lowbit;                                                \
    int        __count;                                                 \
                                                                        \
    ASSERT ((__n & 1) == 1);                                            \
                                                                        \
    __count = GMP_LIMB_BITS-1;                                       \
    do                                                                  \
      {                                                                 \
        __lowbit = __rem & 1;                                           \
        __inv = (__inv >> 1) | (__lowbit << (GMP_LIMB_BITS-1));      \
        __rem = (__rem - (__n & -__lowbit)) >> 1;                       \
      }                                                                 \
    while (-- __count);                                                 \
                                                                        \
    ASSERT (__inv * __n == 1);                                          \
    (inv) = __inv;                                                      \
  } while (0)


double
speed_binvert_limb_mul1 (struct speed_params *s)
{
  SPEED_ROUTINE_MODLIMB_INVERT (binvert_limb_mul1);
}
double
speed_binvert_limb_loop (struct speed_params *s)
{
  SPEED_ROUTINE_MODLIMB_INVERT (binvert_limb_loop);
}
double
speed_binvert_limb_cond (struct speed_params *s)
{
  SPEED_ROUTINE_MODLIMB_INVERT (binvert_limb_cond);
}
double
speed_binvert_limb_arith (struct speed_params *s)
{
  SPEED_ROUTINE_MODLIMB_INVERT (binvert_limb_arith);
}
double
speed_binvert_limb_sec (struct speed_params *s)
{
  SPEED_ROUTINE_MODLIMB_INVERT (binvert_limb_sec);
}
double
speed_binvert_limb_pipe (struct speed_params *s)
{
  SPEED_ROUTINE_MODLIMB_INVERT (binvert_limb_pipe);
}
double
speed_binvert_limb_uintfast (struct speed_params *s)
{
  SPEED_ROUTINE_MODLIMB_INVERT (binvert_limb_uintfast);
}
