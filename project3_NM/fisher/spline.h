// spline.h
//
// Class for creating splines of a table of data
// different class instances for different data sets saving 
// redundancy of code when you need to spline multiple quantities

#ifndef SPLINE_H
#define SPLINE_H

#include "dnumrecipes.h"
#include <sstream>
using namespace std;

class Spline {

 public:
  Spline();
  ~Spline();

  //member function
  int getN();
  double getXMax();
  double getXMin();
  double getXelement(int i);
  int checkInit();

  //
  void setSplineSP(int n1, double x[], double y[]);

  void splineSP(double x[], double y[], int n, double yp1, double ypn,
		double y2[]);
  void splintSP(double xa[], double ya[], double y2a[], int n, double x,
	      double *y);
  double returnValue(double x);
  void hunt(double xx[], int n, double x, int *jlo);

  void loadFileSpline(string file, int ncol, int xcol, int ycol);
  void cleanSpline();

 protected:
  double *xas, *yas, *y2as;
  double xmin,xmax;
  int n;
  int init_flag;

};


#endif
