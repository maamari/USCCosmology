
#ifndef FISHERGAL_CLUSTER_FNL_H
#define FISHERGAL_CLUSTER_FNL_H

#include "fnl_massf.h"
#include "fisherGAL_CLUSTER.h"


using namespace std;

// FisherGAL_CLUSTER_fnl class declaration

class FisherGAL_CLUSTER_FNL : public Fisher
{
 public:
  double zmaxx,intzz;
  double SCCC,fnl;
  int fisher_flag,  dFlag, survey;
  string CLL_name,nameSave,nameSavePS, CLL_namePS, BZfr_name;

  FisherGAL_CLUSTER_FNL();
 ~FisherGAL_CLUSTER_FNL();

 // Cosmology function
 double Omegamz(double z, double wm);
 double OmegaLambdaz(double z, double wl, double w, double wa);
 static double wzKernel(double z, void *params);
 double Eofz(double wm, double wl, double wr, double w, double wa, double z);
 double angdiamdist (double ommhh, double wl, double wr, double w, double wa, double z);
 double comovingdistance(double z, Cosmoparams cosmo);
 double ComovingVolume(double z, double ommhh, double wr, double wl, double w, double wa);
 double getVolume(double z, double ommhh, double wl, double wr, double w, double wa);
 static double k (double z, void * params);
 static double setVol (double z, void * params);
 double growthFac(double z,double wm, double wl, double w);

 // Number Count Fisher Matrix
 double cal_mlim(double ommhh, double wl, double wr, double w, double wa, double z);
 double fff(double x);
 double M200toM500(double M200, double overdensity, double z);
 double MhtoM200(double Mh, double overdensity, double z);
 double NofzFR(double fnl,double wa, double w, double wmmmhh, double wlll, double wrrr, double wbbbhh, double sigmaa8, double nspecc, double zzz, double Mmm1, double Mmm2,double Bm0, double sigma0,double alpha, double beta);
 double    NFR(double fnl,double wa, double w, double wmmmhh, double wlll, double wrrr, double wbbbhh, double sigmaa8, double nspecc, double zzz1,double zzz2, double Mmm1, double Bm0, double sigma0,double alpha, double beta);
 void fisherMatrix();
 void derivatives();

 // Power Spectrum
 double biasPS(double z, double mass,My_params_N my_params);
 double PowSFR(double wa, double w,double wmmm, double wlll, double wbbb,double wrrr,  double sigmaa8, double nspecc, double Bm0, double sigma0, double alpha, double beta,double zzz, double KK);
double ComBiasFR(double kk, double fnl, double mm1, double mm2, double wa, double w, double wmmmhh, double wlll, double wbbbhh,double wrrr,  double sigmaa8, double nspecc, double Bm0,double sigma0,double alpha, double beta, double zzz);
double effPk(double kk , double fnl, double mm1,double MM1, double mm2, double MM2, double wa, double w,double wmmmhh, double wlll, double wbbbhh,double wrrr,  double sigmaa8, double nspecc, double Bm0,double sigma0,double alpha, double beta, double zzz1, double zzz2);

 // Power Spectrum Fisher Matrix
 void derivatives_PS();
 void fisherMatrix_PS();


 void initExperiment(int sur, string label, double SkyCov, double zmax, double intzzz, double FNL, int der_flag);

 protected:

};

#endif
