//fisherBAO.h
// information for calculating the fisher matrix for a galaxy survey
// FisherBAO is a subclass of Fisher
//
// Calculates Fisher matrix for observations at a single specified redshift
// 
// Extra parameters: bias  - one per redshift : cannot connect between 
//						redshift bins
//

#ifndef FISHERBAO_H
#define FISHERBAO_H

#include "spline.h"
#include "fisher.h"
using namespace std;


// FisherBAO class declaration

class FisherBAO : public Fisher
{

 public:
FisherBAO();
~FisherBAO();

//initialisation routines
void initSurvey(string name1, double obsarea, double zmin1, double zmax1, double sigmaz1);

// member functions
string getSurveyName();
double getSurveyZ();
double getZMin();
double getZMax();
double getVolume();
void setCosmology(CosmParam fiducial);

//calculate fisher matrix
void fisherMatrixBAO();
double getFnl(double z);

//systematic errors
void addSystematicErrors(double deltaz);

//add redshift bins
void addRedshiftBin(double zmin1, double zmax1);
void addDistanceMeasurement(double zbin, double sigmaD, double sigmaH);

protected:

  //survey parameters
  double survey_z;
  double volume;
  double obsarea;
  double zmin;
  double zmax;
  double sigmaz;  //redshift errors
  double sigmar;  //distance errors
  double ncal;    //number of photometric z calibration sources

};

#endif
