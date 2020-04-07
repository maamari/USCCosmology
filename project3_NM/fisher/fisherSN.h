//fisherSN.h
// information for calculating the fisher matrix for a galaxy survey
// FisherSN is a subclass of Fisher
//
// Calculates Fisher matrix for observations at a single specified redshift
// 
// Extra parameters: bias  - one per redshift : cannot connect between 
//						redshift bins
//

#ifndef FISHERSN_H
#define FISHERSN_H

#include "spline.h"
#include "fisher.h"
using namespace std;


// FisherSN class declaration

class FisherSN : public Fisher
{

 public:
FisherSN();
~FisherSN();

//initialisation routines
void initSurvey(string name1, double nbin, double zmin1, double zmax1, double sigmaz1, double sigmaD1);

// member functions
string getSurveyName();
double getSurveyZ();
double getZMin();
double getZMax();
double getNbin();
void setCosmology(CosmParam fiducial);

//magnitude calculation
double magnitudeSNfull(double z, CosmParam *c, double magoff, double magL, double magQ, double zbias);
double magnitudeSN(double z, CosmParam *c);

//calculate fisher matrix
void fisherMatrixSN();

//add redshift bins
void addRedshiftBin(double zmin1, double zmax1, double nbin, double sigmaz);

//dark energy
void projectDarkEnergy(int de_flag);
double deDerivatives(string tagn, string tagi);

protected:

  //survey parameters
  double survey_z;
  double volume;
  double nbin;
  double zmin;
  double zmax;
  double sigmaz;  //redshift errors
  double sigmaD;  //dispersion of SN

};

#endif
