// dcomplex.h
// Contains definitions and simple functions of class complex
// Taken from Practical C++ Programming

#ifndef __complex_h__
#define __complex_h__

#include <iostream>
#include <math.h>

class complex
{
 private:
  double real_part;
  double imaginary_part;

 public:
  // Default constructor, zero everything
  complex(void) {
    real_part = 0.0;
    imaginary_part = 0.0;
  }

  // Copy constructor
  complex(const complex& other_complex) {
    real_part = other_complex.real_part;
    imaginary_part = other_complex.imaginary_part;
  }

  // Construct a complex out of two real numbers
  complex(double init_real, double init_imaginary = 0.0) {
    real_part = init_real;
    imaginary_part = init_imaginary;
  }

  // Destructor does nothing
  ~complex() {}

  // Functions to return pieces of the number
  double real(void) const {
    return (real_part);
  }

  double imaginary(void) const {
    return (imaginary_part);
  }

  // Functions to set parts of number
  void set(double real, double imaginary) {
    real_part = real;
    imaginary_part = imaginary;
  }

  void set_real(double real) {
    real_part = real;
  }

  void set_imaginary(double imaginary) {
    imaginary_part = imaginary;
  }

  complex operator = (const complex& oper2) {
    set(oper2.real_part,oper2.imaginary_part);
    return (*this);
  }

  complex& operator +=(const complex& oper2) {
    real_part += oper2.real();
    imaginary_part += oper2.imaginary();
    return (*this);
  }

  complex& operator += (double oper2) {
    real_part += oper2;
    return (*this);
  }

  complex& operator -= (const complex& oper2) {
    real_part -= oper2.real();
    imaginary_part -= oper2.imaginary();
    return (*this);
  }

  complex& operator -= (double oper2) {
    real_part -= oper2;
    return (*this);
  }

  complex& operator *= (const complex& oper2) {
    // Place to hold the real part of the result while we compute imag part
    double real_result = (real_part*oper2.real() - 
			  imaginary_part*oper2.imaginary());
    imaginary_part = (real_part*oper2.imaginary() +
		      imaginary_part*oper2.real());
    real_part = real_result;
    return (*this);
  }

  complex& operator *= (double oper2) {
    real_part *= oper2;
    imaginary_part *= oper2;
    return (*this);
  }

  complex& operator /= (const complex& oper2);

  complex& operator /= (double oper2) {
    real_part /= oper2;
    imaginary_part /= oper2;
    return (*this);
  }
};

inline complex operator + (const complex& oper1, const complex& oper2)
{
  return complex(oper1.real() + oper2.real(),
		 oper1.imaginary() + oper2.imaginary());
}

inline complex operator + (const complex& oper1, double oper2)
{
  return complex(oper1.real() + oper2,
		 oper1.imaginary());
}

inline complex operator + (double oper1, const complex& oper2)
{
  return complex(oper1 + oper2.real(),
		 oper2.imaginary());
}

inline complex operator - (const complex& oper1, const complex& oper2)
{
  return complex(oper1.real() - oper2.real(),
		 oper1.imaginary() - oper2.imaginary());
}

inline complex operator - (const complex& oper1, double oper2)
{
  return complex(oper1.real() - oper2,
		 oper1.imaginary());
}

inline complex operator - (double oper1, const complex& oper2)
{
  return complex(oper1 - oper2.real(),
		 -oper2.imaginary());
}

inline complex operator * (const complex& oper1, const complex& oper2)
{
  return complex(oper1.real()*oper2.real() - oper1.imaginary()*
		 oper2.imaginary(),
		 oper1.real()*oper2.imaginary() + oper1.imaginary()*
		 oper2.real());
}

inline complex operator * (const complex& oper1, const double oper2)
{
  return complex(oper1.real()*oper2,
		 oper1.imaginary()*oper2);
}

inline complex operator * (const double oper1, const complex& oper2)
{
  return complex(oper1*oper2.real(),
		 oper1*oper2.imaginary());
}

extern complex operator / (const complex &oper1, const complex &oper2);

inline complex operator /  (const double &oper1, const complex &oper2) {
  return (complex(oper1,0.0)/oper2);
}

inline complex operator / (const complex &oper1, const double &oper2) {
  return (oper1/complex(oper2,0.0));
}

inline int operator == (const complex& oper1, const complex& oper2) {
  return ((oper1.real() == oper2.real()) && 
	  (oper1.imaginary() == oper2.imaginary()));
}

inline int operator != (const  complex& oper1, const complex& oper2) {
  return (!(oper1 == oper2));
}

inline complex operator - (const complex& oper1)
{
  return complex(-oper1.real(), -oper1.imaginary());
}

inline complex operator + (const complex& oper1)
{
  return complex(+oper1.real(), +oper1.imaginary());
}

inline complex conjg(const complex oper1)
{
  return complex(oper1.real(), -oper1.imaginary());
}

extern double Cabs(const complex oper1);

/*
inline ostream &operator << (ostream &out_file, const complex &number)
{
  out_file << '{' << number.real() << ',' << number.imaginary()
	   << '}';
  return (out_file);
}
*/
/* extern istream &operator >> (istream &in_file, complex &number); */

#endif /* __complex_h__ */
