#ifndef __MATRIX_H
#define __MATRIX_H
#include <math.h>

const double ident[9] = {1,0,0,0,1,0,0,0,1};

// **** Matrix manipulation, all done in 9-v notation

// DET: Returns the det. of an integer matrix
inline int det (const int x[9]) 
{
  return (x[0]*(x[4]*x[8] - x[5]*x[7])
    + x[1]*(x[6]*x[5] - x[3]*x[8])
    + x[2]*(x[3]*x[7] - x[6]*x[4]));
}

inline double det (const double x[9]) 
{
  return (x[0]*(x[4]*x[8] - x[5]*x[7])
    + x[1]*(x[6]*x[5] - x[3]*x[8])
    + x[2]*(x[3]*x[7] - x[6]*x[4]));
}

// MULT: evaluates the product of a * b, returns in c:
inline void mult (const int a[9], const int b[9], int c[9]) 
{
  c[0] = a[0]*b[0] + a[1]*b[3] + a[2]*b[6];
  c[1] = a[0]*b[1] + a[1]*b[4] + a[2]*b[7];
  c[2] = a[0]*b[2] + a[1]*b[5] + a[2]*b[8];
  c[3] = a[3]*b[0] + a[4]*b[3] + a[5]*b[6];
  c[4] = a[3]*b[1] + a[4]*b[4] + a[5]*b[7];
  c[5] = a[3]*b[2] + a[4]*b[5] + a[5]*b[8];
  c[6] = a[6]*b[0] + a[7]*b[3] + a[8]*b[6];
  c[7] = a[6]*b[1] + a[7]*b[4] + a[8]*b[7];
  c[8] = a[6]*b[2] + a[7]*b[5] + a[8]*b[8];
}

inline void mult (const double a[9], const int b[9], double c[9]) 
{
  c[0] = a[0]*b[0] + a[1]*b[3] + a[2]*b[6];
  c[1] = a[0]*b[1] + a[1]*b[4] + a[2]*b[7];
  c[2] = a[0]*b[2] + a[1]*b[5] + a[2]*b[8];
  c[3] = a[3]*b[0] + a[4]*b[3] + a[5]*b[6];
  c[4] = a[3]*b[1] + a[4]*b[4] + a[5]*b[7];
  c[5] = a[3]*b[2] + a[4]*b[5] + a[5]*b[8];
  c[6] = a[6]*b[0] + a[7]*b[3] + a[8]*b[6];
  c[7] = a[6]*b[1] + a[7]*b[4] + a[8]*b[7];
  c[8] = a[6]*b[2] + a[7]*b[5] + a[8]*b[8];
}

inline void mult (const int a[9], const double b[9], double c[9]) 
{
  c[0] = a[0]*b[0] + a[1]*b[3] + a[2]*b[6];
  c[1] = a[0]*b[1] + a[1]*b[4] + a[2]*b[7];
  c[2] = a[0]*b[2] + a[1]*b[5] + a[2]*b[8];
  c[3] = a[3]*b[0] + a[4]*b[3] + a[5]*b[6];
  c[4] = a[3]*b[1] + a[4]*b[4] + a[5]*b[7];
  c[5] = a[3]*b[2] + a[4]*b[5] + a[5]*b[8];
  c[6] = a[6]*b[0] + a[7]*b[3] + a[8]*b[6];
  c[7] = a[6]*b[1] + a[7]*b[4] + a[8]*b[7];
  c[8] = a[6]*b[2] + a[7]*b[5] + a[8]*b[8];
}

inline void mult (const double a[9], const double b[9], double c[9]) 
{
  c[0] = a[0]*b[0] + a[1]*b[3] + a[2]*b[6];
  c[1] = a[0]*b[1] + a[1]*b[4] + a[2]*b[7];
  c[2] = a[0]*b[2] + a[1]*b[5] + a[2]*b[8];
  c[3] = a[3]*b[0] + a[4]*b[3] + a[5]*b[6];
  c[4] = a[3]*b[1] + a[4]*b[4] + a[5]*b[7];
  c[5] = a[3]*b[2] + a[4]*b[5] + a[5]*b[8];
  c[6] = a[6]*b[0] + a[7]*b[3] + a[8]*b[6];
  c[7] = a[6]*b[1] + a[7]*b[4] + a[8]*b[7];
  c[8] = a[6]*b[2] + a[7]*b[5] + a[8]*b[8];
}

// MULT3 routines: D = ABC
// only dii is a little different... the rest are all the same
inline void mult (const int a[9], const int b[9], const int c[9],
		  int d[9]) 
{ int temp[9];  mult(a,b,temp);  mult(temp,c,d);}

inline void mult (const double a[9], const int b[9], const int c[9],
		  double d[9]) 
{ int temp[9];  mult(b,c,temp);  mult(a, temp,d);}

inline void mult (const int a[9], const double b[9], const int c[9],
		  double d[9]) 
{ double temp[9];  mult(a,b,temp);  mult(temp,c,d);}

inline void mult (const double a[9], const double b[9], const int c[9],
		  double d[9]) 
{ double temp[9];  mult(a,b,temp);  mult(temp,c,d);}

inline void mult (const int a[9], const int b[9], const double c[9],
		  double d[9]) 
{ int temp[9];  mult(a,b,temp);  mult(temp,c,d);}

inline void mult (const double a[9], const int b[9], const double c[9],
		  double d[9]) 
{ double temp[9];  mult(a,b,temp);  mult(temp,c,d);}

inline void mult (const int a[9], const double b[9], const double c[9],
		  double d[9]) 
{ double temp[9];  mult(a,b,temp);  mult(temp,c,d);}

inline void mult (const double a[9], const double b[9], const double c[9],
		  double d[9]) 
{ double temp[9];  mult(a,b,temp);  mult(temp,c,d);}

// MULT_SCALAR: Multiplies matrix a by scalar b, result in c
inline void mult (const int a[9], const int b, int c[9]) 
{ for (int i=0; i<9; ++i) c[i] = b*a[i]; }

inline void mult (const int b, const int a[9], int c[9]) 
{ for (int i=0; i<9; ++i) c[i] = b*a[i]; }

inline void mult (const int a[9], const double b, double c[9]) 
{ for (int i=0; i<9; ++i) c[i] = b*a[i]; }

inline void mult (const double b, const int a[9], double c[9]) 
{ for (int i=0; i<9; ++i) c[i] = b*a[i]; }

inline void mult (const double a[9], const int b, double c[9]) 
{ for (int i=0; i<9; ++i) c[i] = b*a[i]; }

inline void mult (const int b, const double a[9], double c[9]) 
{ for (int i=0; i<9; ++i) c[i] = b*a[i]; }

inline void mult (const double a[9], const double b, double c[9]) 
{ for (int i=0; i<9; ++i) c[i] = b*a[i]; }

inline void mult (const double b, const double a[9], double c[9]) 
{ for (int i=0; i<9; ++i) c[i] = b*a[i]; }


// SQUARE: Calculate xT * x
//   Note: the product is symmetric, so we get to do a smaller number
//   of multiplications.

inline void square (const int x[9], int s[9]) 
{
  s[0] = x[0]*x[0] + x[3]*x[3] + x[6]*x[6];
  s[1] = x[0]*x[1] + x[3]*x[4] + x[6]*x[7];
  s[2] = x[0]*x[2] + x[3]*x[5] + x[6]*x[8];
  s[3] = s[1];
  s[4] = x[1]*x[1] + x[4]*x[4] + x[7]*x[7];
  s[5] = x[1]*x[2] + x[4]*x[5] + x[7]*x[8];
  s[6] = s[2];
  s[7] = s[5];
  s[8] = x[2]*x[2] + x[5]*x[5] + x[8]*x[8];
}

inline void square (const int x[9], double s[9]) 
{
  s[0] = x[0]*x[0] + x[3]*x[3] + x[6]*x[6];
  s[1] = x[0]*x[1] + x[3]*x[4] + x[6]*x[7];
  s[2] = x[0]*x[2] + x[3]*x[5] + x[6]*x[8];
  s[3] = s[1];
  s[4] = x[1]*x[1] + x[4]*x[4] + x[7]*x[7];
  s[5] = x[1]*x[2] + x[4]*x[5] + x[7]*x[8];
  s[6] = s[2];
  s[7] = s[5];
  s[8] = x[2]*x[2] + x[5]*x[5] + x[8]*x[8];
}

inline void square (const double x[9], double s[9]) 
{
  s[0] = x[0]*x[0] + x[3]*x[3] + x[6]*x[6];
  s[1] = x[0]*x[1] + x[3]*x[4] + x[6]*x[7];
  s[2] = x[0]*x[2] + x[3]*x[5] + x[6]*x[8];
  s[3] = s[1];
  s[4] = x[1]*x[1] + x[4]*x[4] + x[7]*x[7];
  s[5] = x[1]*x[2] + x[4]*x[5] + x[7]*x[8];
  s[6] = s[2];
  s[7] = s[5];
  s[8] = x[2]*x[2] + x[5]*x[5] + x[8]*x[8];
}

// INVERSE: Evaluates inverse of x; returns it in inv / inverse (return value
// is det x).  Need to divide by det to get correct inverse.
inline int inverse (const int x[9], int inv[9])
{
  inv[0] = x[4]*x[8] - x[5]*x[7];
  inv[1] = x[2]*x[7] - x[1]*x[8];
  inv[2] = x[1]*x[5] - x[2]*x[4];

  inv[3] = x[5]*x[6] - x[3]*x[8];
  inv[4] = x[0]*x[8] - x[2]*x[6];
  inv[5] = x[2]*x[3] - x[0]*x[5];

  inv[6] = x[3]*x[7] - x[4]*x[6];
  inv[7] = x[1]*x[6] - x[0]*x[7];
  inv[8] = x[0]*x[4] - x[1]*x[3];
  
  return det(x);
}

inline double inverse (const double x[9], double inv[9])
{
  inv[0] = x[4]*x[8] - x[5]*x[7];
  inv[1] = x[2]*x[7] - x[1]*x[8];
  inv[2] = x[1]*x[5] - x[2]*x[4];

  inv[3] = x[5]*x[6] - x[3]*x[8];
  inv[4] = x[0]*x[8] - x[2]*x[6];
  inv[5] = x[2]*x[3] - x[0]*x[5];

  inv[6] = x[3]*x[7] - x[4]*x[6];
  inv[7] = x[1]*x[6] - x[0]*x[7];
  inv[8] = x[0]*x[4] - x[1]*x[3];
  
  return det(x);
}

// ********************************* magnsq ****************************
// Returns the squared magnitude of a vector, where metric[9] is the
// SYMMETRIC metric, and u is the vector in direct coordinates
inline double magnsq (const double metric[9], const double u[3]) 
{
  return metric[0]*u[0]*u[0] + metric[4]*u[1]*u[1] + metric[8]*u[2]*u[2]
    + 2*metric[1]*u[0]*u[1] + 2*metric[2]*u[0]*u[2] + 2*metric[5]*u[1]*u[2];
}

inline double magnsq (const int metric[9], const double u[3]) 
{
  return metric[0]*u[0]*u[0] + metric[4]*u[1]*u[1] + metric[8]*u[2]*u[2]
    + 2*metric[1]*u[0]*u[1] + 2*metric[2]*u[0]*u[2] + 2*metric[5]*u[1]*u[2];
}

inline double magnsq (const double metric[9], const int u[3]) 
{
  return metric[0]*u[0]*u[0] + metric[4]*u[1]*u[1] + metric[8]*u[2]*u[2]
    + 2*metric[1]*u[0]*u[1] + 2*metric[2]*u[0]*u[2] + 2*metric[5]*u[1]*u[2];
}

inline int magnsq (const int metric[9], const int u[3]) 
{
  return metric[0]*u[0]*u[0] + metric[4]*u[1]*u[1] + metric[8]*u[2]*u[2]
    + 2*metric[1]*u[0]*u[1] + 2*metric[2]*u[0]*u[2] + 2*metric[5]*u[1]*u[2];
}

// ZERO: determine if matrix is zero
inline int zero (const int a[9]) 
{
  return (a[0]==0) && (a[1]==0) &&(a[2]==0) &&
    (a[3]==0) && (a[4]==0) && (a[5]==0) &&
    (a[6]==0) && (a[7]==0) && (a[8]==0);
}

inline int zero (const double a[9]) 
{
  return zero(a[0]) && zero(a[1]) && zero(a[2]) &&
    zero(a[3]) && zero(a[4]) && zero(a[5]) &&
    zero(a[6]) && zero(a[7]) && zero(a[8]);
}

// ZERO_VECT: same thing for vectors
inline int zero_vect(const int a[3]) {
  return (a[0]==0) && (a[1]==0) && (a[2]==0);
}

inline int zero_vect(const double a[3]) {
  return zero(a[0]) && zero(a[1]) && zero(a[2]);
}

// CAREFUL_INVERSE: inverts a matrix, symmetric or not, but polishes
// the final result to within a set tolerance.
// NOTE: unlike all of our other inverses, we go ahead and scale by
// the determinant.
const double CAREFUL_INV_TOL = 1e-15;
const int MAX_ITER = 20; // maximum number of iterations

void careful_inverse (const double a[9], double inv[9], double &err) 
{
  double b[9], ba[9], bab[9];
  double deta;
  int d;
  
  deta = inverse(a, b); // our first crack at it...
  if (deta == 0) {
    err = 1; // a pretty obvious flag...
    return;
  }
  mult(b, 1./deta, b); // now b should be scaled accordingly.
  
  mult(b, a, ba);
  for (err=0, d=0; d<9; ++d) err += fabs(ba[d]-ident[d]);
  // polish the result, until err is below our desired tolerance
  // but don't do more than n steps.
  for (int n=0; (n<MAX_ITER) && (err > CAREFUL_INV_TOL); ++n) {
    mult(ba, b, bab);
    // polish in place:
    for (d=0; d<9; ++d) b[d] = 2.*b[d] - bab[d];
    mult(b, a, ba);
    for (err=0, d=0; d<9; ++d) err += fabs(ba[d]-ident[d]);
  }
  // final result:
  for (d=0; d<9; ++d) inv[d] = b[d];
}

inline void careful_inverse (const double a[9], double inv[9]) 
{
    double err;
    careful_inverse(a, inv, err);
}

// Do the transpose "in place," with a loop "unrolled" version:

#define __ISWAP__(x,y) ({temp = x; x = y; y = temp;})
#define __DSWAP__(x,y) ({temp = x; x = y; y = temp;})

inline void self_transpose (int a[9]) 
{
    int temp;
    __ISWAP__(a[1], a[3]);
    __ISWAP__(a[2], a[6]);
    __ISWAP__(a[5], a[7]);
}

inline void self_transpose (double a[9]) 
{
    double temp;
    __DSWAP__(a[1], a[3]);
    __DSWAP__(a[2], a[6]);
    __DSWAP__(a[5], a[7]);
}

#undef __ISWAP__
#undef __DSWAP__

#endif