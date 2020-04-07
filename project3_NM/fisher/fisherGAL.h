//fisherGAL.h
// information for calculating the fisher matrix for a galaxy survey
// FisherGAL is a subclass of Fisher
//
// Calculates Fisher matrix for observations at a single specified redshift
// 
// Extra parameters: bias  - one per redshift : cannot connect between 
//						redshift bins
//

#ifndef FISHERGAL_H
#define FISHERGAL_H


#include "spline.h"
#include "fisher.h"
using namespace std;


//const int NPRECISION(14); //sigfig for file io
const int nparamGALX(3);
const int nl_flag(1);

// FisherGAL class declaration

class FisherGAL : public Fisher
{

 public:
FisherGAL();
~FisherGAL();

//initialisation routines
void initSurvey(string name1, double z1, double volume1, double density1, double kmin1, double kmax1, double sigma8gal1, double sigmaz1);

// member functions
string getSurveyName();
double getSurveyZ();
double getKMin();
double getKMax();
double getKNL(double z);
double getKMaxUse();
double getVolume();
double getDensity();
double getSigperp2();
double getSigpara2();

//cosmology supplemental routines
double getBeta(double z, CosmParam *c);
double getBias(CosmParam *c);

//set galaxy files
double setGalaxyFiles();
double setGalaxyPair(int iflag);

//calculate fisher matrix
void fisherMatrixGAL();

// calculate experimental uncertainty
double volumeEffective(double k, double mu);

//calculate derivatives of power spectrum
double derivePowerGAL(double k, double mu, CosmParam c_p, CosmParam c_m, string tag, int iflag);
void setDerivSplines(string tag1, string tag2, int iflag);

//power spectrum code
double getPowerSpectrumObs(double kref, double muref,CosmParam *c, Spline *spline);
double powerSpectrumError(double k);
void outPowerErrors();
double averagedEffectiveVolume(double k);
void derivGAL(string tag,string file);



protected:

  //survey parameters
  double survey_z;
  double volume;
  double density;
  double kmin;
  double kmax;
  double sigma8gal;
  double sigmaz;  //redshift errors
  double sigmar;  //distance errors
  double sigpara2;
  double sigperp2;

  //redshift parameters for galaxy surveys
 // double lda;
 // double lh;
 // double lg;
 // double lbeta;
 // double pshot;
 // double beta;

  //growth parameters to speed calculation
//  double growthZ;
//  double growth0;
//  double daratio;
//  double hratio;

  //splines
  Spline splineFid;
  Spline splinePow1p;
  Spline splinePow1m;
  Spline splinePow2p;
  Spline splinePow2m;

  //flags
 // int nl_flag;

};

//fisher calculation
double getFisherKernelK(double lk);
double getFisherKernel(double mu);
double setFisherKernel(double mu, string file_fid, FisherGAL *survey1,string file_dlPi1, string file_dlPj1,int iflag);

double setKErrorKernelGAL(double k1, double mu, FisherGAL *survey1, int iflag);
double getKErrorKernelGAL(double theta);

double setKVolumeKernelGAL(double k1, double mu, FisherGAL *survey1, int iflag);
double getKVolumeKernelGAL(double theta);

#endif
