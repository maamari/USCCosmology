//inflation.h
//
//  Class for inflation using slow roll parameters
//

#ifndef INFLATION_H
#define INFLATION_H

const double EULER_CONSTANT(0.5772156649);  //Euler-Mascheroni constant

#include "spline.h"
using namespace std;

// Inflation class declaration

class Inflation
{

 public:
Inflation();
~Inflation();

//spectral parameters
double getNscalar(double epsilon, double eta, double xi, double lambda4);
double getTS(double epsilon, double eta, double xi, double lambda4);
double getRunning(double epsilon, double eta, double xi, double lambda4);
double getNtensor(double epsilon, double eta, double xi, double lambda4);

 protected:

};

#endif