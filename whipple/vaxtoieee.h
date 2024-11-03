/** Defines **/
#define IEEE_F_BIAS 127
#define IEEE_D_BIAS 1023
#define VAX_F_BIAS  128

/** Global Variables **/
struct vax_fbits{
    unsigned high_mantissa:  7;
    unsigned exponent     :  8;
    unsigned sign         :  1;
    unsigned low_mantissa : 16;
};
struct vax_dbits{
    unsigned a_man    :  7;
    unsigned exponent :  8;
    unsigned sign     :  1;
    unsigned b_man    : 16;
    unsigned c_man    : 16;
    unsigned d_man    : 16;
};
struct ieee_fbits{
    unsigned mantissa : 23;
    unsigned exponent :  8;
    unsigned sign     :  1;
};
struct ieee_dbits{
    unsigned low_mantissa  : 32;
    unsigned high_mantissa : 20;
    unsigned exponent      : 11;
    unsigned sign          :  1;
};
union vax_ffloat{
  float f;
  struct vax_fbits b;
};
union vax_dfloat{
  double f;
  struct vax_dbits b;
};
union ieee_ffloat{
  float f;
  struct ieee_fbits b;
};
union ieee_dfloat{
  double f;
  struct ieee_dbits b;
};

/** Function Declarations **/
void ieeetovax_double();
void vaxtoieee_double();
void vaxtoieee_float();
void ieeetovax_float();
