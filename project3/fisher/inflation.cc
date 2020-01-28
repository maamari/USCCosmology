/* inflation.cc
 * Contains functions for calculating inflationary parameters
  */

#include <math.h>
#include <iostream>
#include <fstream>
#include "dcosmology2.h"
#include "dnumrecipes.h"
#include "spline.h"
#include <sstream>
#include "inflation.h"

using namespace std;

/*********************************************************************
 **************** Constructor/Destructor *****************************
 ********************************************************************/

 Inflation::Inflation()
{
	//cout<<"Inflation constructor called"<<endl;
}

 Inflation::~Inflation()
{
	//cout<<"Inflation destructor called"<<endl;

}

/////////////////////////////////////////////////////////////////////
// Convert slow roll parameters to spectral parameters
/////////////////////////////////////////////////////////////////////
double Inflation::getNscalar(double epsilon, double eta, double xi, double lambda4)
{
	double C(4.0*(EULER_CONSTANT+log(2.0))-5.0);
	double nscal;

	nscal=1.0+2.0*eta-4.0*epsilon;
	nscal-=2.0*(1.0+C)*epsilon*epsilon;
	nscal-=0.5*(3.0-5.0*C)*epsilon*eta;
	nscal+=0.5*(3.0-C)*xi;

	return nscal;
}

double Inflation::getTS(double epsilon, double eta, double xi, double lambda4)
{
	double C(-2.0+EULER_CONSTANT+log(2.0));
	double ts;

	ts=16.0*epsilon*(1.0+2.0*C*(epsilon-eta));

	return ts;
}

double Inflation::getRunning(double epsilon, double eta, double xi, double lambda4)
{
	double C(4.0*(EULER_CONSTANT+log(2.0))-5.0);
	double alpha;
	double detadN, depsilondN, dxidN;

	depsilondN=2.0*epsilon*(eta-epsilon);
	detadN= -epsilon*eta+xi;
	dxidN=xi*(eta-2.0*epsilon)+lambda4;

	alpha=2.0*detadN-4.0*depsilondN-4.0*(1.0+C)*epsilon*detadN;
 	alpha-=0.5*(3.0-5.0*C)*(epsilon*detadN+eta*depsilondN);
	alpha+=0.5*(3.0-C)*dxidN;
	alpha/= -(1.0-epsilon);

	return alpha;
}

double Inflation::getNtensor(double epsilon, double eta, double xi, double lambda4)
{
	double C(4.0*(EULER_CONSTANT+log(2.0))-5.0);
	double nt;

	nt= -2.0*epsilon-(3.0+C)*epsilon*epsilon;
	nt+=(1.0+C)*epsilon*eta;

	return nt;
}