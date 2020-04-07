// dcomplex.cc
// Includes complicated function definitions for complex class

#include "dcomplex.h"

// Complex division
complex operator / (const complex& oper1, const complex& oper2) 
{
  // denominator
  double den = fabs(oper2.real()) + fabs(oper2.imaginary());

  // real and imaginary parts of the oper1 factor
  double oper1_real_den = oper1.real()/den;
  double oper1_imag_den = oper1.imaginary()/den;

  // real and imaginary parts of the oper2 factor
  double oper2_real_den = oper2.real()/den;
  double oper2_imag_den = oper2.imaginary()/den;

  // normalization
  double normalization = (oper2_real_den*oper2_real_den +
			  oper2_imag_den*oper2_imag_den);

  return complex((oper1_real_den*oper2_real_den + 
		  oper1_imag_den*oper2_imag_den)/normalization,
		 (oper1_imag_den*oper2_real_den -
		  oper1_real_den*oper2_imag_den)/normalization);
}

// Complex divide by
complex& complex::operator /= (const complex& oper2)
{
  // denominator
  double den = fabs(oper2.real()) + fabs(oper2.imaginary());

  // real and imaginary parts of the oper1 factor
  double oper1_real_den = real_part/den;
  double oper1_imag_den = imaginary_part/den;

  // real and imaginary parts of the oper2 factor
  double oper2_real_den = oper2.real()/den;
  double oper2_imag_den = oper2.imaginary()/den;

  // normalization
  double normalization = (oper2_real_den*oper2_real_den +
			  oper2_imag_den*oper2_imag_den);

  real_part = (oper1_real_den*oper2_real_den + 
	       oper1_imag_den*oper2_imag_den)/normalization;
  imaginary_part = (oper1_imag_den*oper2_real_den -
		    oper1_real_den*oper2_imag_den)/normalization;
  return (*this);
}

/*
istream &operator >> (istream &in_file, complex &number)
{
  double real,imaginary;
  char ch;      // Random character used to verify input

  number.set(0.0,0.0);

  //  in_file.ipfx(1);     // Tell the I/O system we are reading formatted
  in_file >> ws;       // Skip white space

  if (in_file.bad()) return (in_file);
  in_file >> ch;       // Get character after white space
  if (ch != '{') {
    in_file.setf(ios::failbit);     // We have an error
    return (in_file);
  }

  in_file >> real;

  if (in_file.bad()) return (in_file);
  if (ch != ',') {
    in_file.setf(ios::failbit);     // We have an error
    return (in_file);
  }

  in_file >> imaginary;

  in_file >> ws >> ch;
  if (ch != '}') {
    in_file.setf(ios::failbit);     // We have an error
    return (in_file);
  }

  number.set(real,imaginary);
  return (in_file);
}
*/

// Returns absolute value of complex number
double Cabs(const complex oper1)
{
  double x,y,ans,temp;
  
  x = fabs(oper1.real());
  y = fabs(oper1.imaginary());

  if (x == 0.0) {
    ans = y;
  }
  else if (y == 0.0) {
    ans = x;
  }
  else if (x > y) {
    temp = y/x;
    ans = x*sqrt(1.0+temp*temp);
  } else {
    temp = x/y;
    ans = y*sqrt(1.0+temp*temp);
  }
  return ans;
}
