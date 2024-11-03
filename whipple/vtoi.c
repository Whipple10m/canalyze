/* VAX to IEEE floating point representation conversion. Both floating
   (single precision) and double (double prcision).

   Written by:
    R.W. Lessard, April 1994
    University College Dublin
    Department of Experimental Physics
    lessard@ferdia.ucd.ie

   References:
    Stevenson, D., "A Proposed Standard for Binary Floating-Point Arithmetic,
     Draft 8.0 of IEEE Task P754", IEEE Computer, 14(3), pg 51, March 1981.
    Brink, J., Spillman, R., Computer Architecture and VAX Assembly Language,
     The Benjamin/Cummings Publishing Company, 1987.
    Levy, H.M., Eckhouse, R.H. Jr., Computer Programming and Architecture -
     The VAX-11, Digital Press, 1980.
     
VAX ->	value = (-1)^s * 0.f * 2^(exponent - 128)
Bias = 128
Special cases:
 1) exponent = 0, sign = 0, value = zero (exactly)
 2) exponent = 0, sign = 1, value = error
Little endian bit pattern:
 FLOAT:
 high_mantissa(7):exponent(8):sign(1):low_mantissa(16)
 DOUBLE:
 a_mantissa(7):exponent(8):sign(1):b_mantissa(16):c_mantissa(16):d_mantissa(16)
Note:
 f is always normalised, i.e. the first bit of the mantissa must be one
Smallest float/double: 0.100..00 * 2^(-127) = 2.939e(-39)
Largest  float/double: 0.111..11 * 2^(127)  = 1.7e(39)

IEEE ->
Special cases:
 FLOAT
  if e = 255 and f <> 0 value = NAN (Not a number)
     e = 255 and f = 0  value = (-1)^s * Infinity
     0 < e < 255        value = (-1)^s * 1.f * 2^(exponent - 127)
     e = 0 and f <> 0   value = (-1)^s * 0.f * 2^(-126)
     e = 0 and f = 0    value = (-1)^s * 0 (exactly)
  Bias = 127
  Little endian bit pattern:
   mantissa(23):exponent(8):sign(1)
 DOUBLE
  if e = 2047 and f <> 0 value = NAN (Not a number)
     e = 2047 and f = 0  value = (-1)^s * Infinity
     0 < e < 2047        value = (-1)^s * 1.f * 2^(exponent - 1023)
     e = 0 and f <> 0    value = (-1)^s * 0.f * 2^(-1022)
     e = 0 and f = 0     value = (-1)^s * 0 (exactly)
  Bias = 1023
  Little endian bit pattern:
   low_mantissa(32):high_mantissa(20):exponent(11):sign(1)
*/
/** Includes **/
#include <stdio.h>
#include "vaxtoieee.h"

void vaxtoieee_float(float *x)
{
  union ieee_ffloat v;
  union vax_ffloat a;
  unsigned unbiased_exponent;

  a.f = *x;
  v.b.sign = a.b.sign;
  switch(a.b.exponent){
  case 0:
    if(a.b.sign == 0){/*VAX defined zero*/
      v.b.exponent = 0;
      v.b.mantissa = 0;
    }else{
      fprintf(stderr,"VAX FLOATING POINT ERROR, INCORRECT VAX REPRESENTATION");
      v.b.exponent = 0;
      v.b.mantissa = 0;
    }
    break;
  case 1:
    unbiased_exponent = a.b.exponent - VAX_F_BIAS;
    v.b.exponent = unbiased_exponent + IEEE_F_BIAS;
    v.b.mantissa = a.b.high_mantissa << 16;
    v.b.mantissa = v.b.mantissa | a.b.low_mantissa;
    /* Put back the VAX hidden bit */
    v.b.mantissa >>= 1;
    v.b.mantissa |= 0x400000;
    break;
  default:
    unbiased_exponent = a.b.exponent - VAX_F_BIAS;
    v.b.exponent = unbiased_exponent + IEEE_F_BIAS - 1;
    v.b.mantissa = a.b.high_mantissa << 16;
    v.b.mantissa = v.b.mantissa | a.b.low_mantissa;
  break;
  }
  *x = v.f;

  return;
}

void vaxtoieee_double(double *x)
{
  union ieee_dfloat v;
  union vax_dfloat a;
  unsigned unbiased_exponent;

  a.f = *x;
  v.b.sign = a.b.sign;
  unbiased_exponent = a.b.exponent - VAX_F_BIAS;
  v.b.exponent = unbiased_exponent + IEEE_D_BIAS - 1;
  v.b.high_mantissa = (a.b.a_man << 13) | (a.b.b_man >> 3);
  v.b.low_mantissa = (a.b.b_man << 29) | (a.b.c_man << 13) |
    (a.b.d_man >> 3);

  *x = v.f;

  return;
}
