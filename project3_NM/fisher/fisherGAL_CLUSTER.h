
#ifndef FISHERGAL_CLUSTER_H
#define FISHERGAL_CLUSTER_H

#include "fisher.h"
#include "matrix.h"
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <math.h>
#include <iostream>
#include <valarray>
#include <fstream>
#include <sstream>
#include "dcosmology2.h"
#include <gsl/gsl_sf_coupling.h>
#include <gsl/gsl_sf_result.h>
#include "dnumrecipes.h"
#include "spline.h"


using namespace std;

// Constants
#ifdef TCMB
#define TCMB 2.725
#endif
#ifndef Qe
#define Qe 1.60217646e-19         // electron charge in Coulombs
#endif
#ifndef cl
#define cl 299792458.0           // speed of light in m/s
#endif
#ifndef me
#define me 9.19838188e-31        // electron mass in kg
#endif
#ifndef me2
#define me2 0.511                // electron mass in MeV/c2
#endif
#ifndef sigT
#define sigT 6.6524e-29          // Thompson scattering cross section in m^2
#endif
#ifndef kB
#define kB 1.3806503e-23         // Boltzmann's constant in m^2 kg s^-2 K^-1
#endif
#ifndef hp
#define hp 6.626068e-34          // Planck's constant in J/Hz
#endif
#ifndef G
#define G 6.674e-11              // Gravtational constant in N (m/kg)^2
#endif
#ifndef G2
#define G2 4.515e-48             // Gravitational constant in Mpc^3 Msun^-1 s^-2
#endif
#ifndef PI
#define PI 4.0*atan(1.0)          // Pi
#endif
#ifndef pi
#define pi PI
#endif
#ifndef d2r
#define d2r (pi/180.0)
#endif
#ifndef convd
#define convd 3.2477649e-23      // m to Mpc
#endif
#ifndef Mpc
#define Mpc 3.085678e22          // Mpc to m
#endif
#ifndef Msun
#define Msun 1.99e30             // Mass of Sun in kg
#endif
#ifndef fsdegree
#define fsdegree 41252.96125       // full sky sq degree
#endif


// Cosmological parameters
//#ifndef sigma8
//#define sigma8 0.796             // mass variance at R=8h^-1 Mpc
//#endif
//#ifndef h
//#define h 0.72                  // Hubble scaling
//#endif

class Cosmo;


class My_params {
  public:
    double z, M, r, k, moment, fr0;
    int l;
    // Cosmo cosmo;
    double wm, wl, wr, wb, sigma8, ns,Bm0, sigma0,alpha,beta,w,wa,h;
    double m1,m2;
    bool frflag;
    string filename, mdyfile;
    int survey;
    CosmParam c;
    //FisherGAL_CLUSTER *cluster;
};


// FisherGAL_CLUSTER class declaration

class FisherGAL_CLUSTER : public Fisher
{
 public:
  double zmaxx,intzz;
  double SCCC,fr_0;
  int fisher_flag,  dFlag, survey;
  string Mdyn_name,CLL_name,nameSave,nameSavePS,FR_name,nameSaveFR,BZfr_name, CLL_namePS, FR_bias, FR_Pk;

  FisherGAL_CLUSTER();
 ~FisherGAL_CLUSTER();

 double getParamFromTag_cluster(string tag, CosmParam *c);
 void specifyParams_cluster(int *choiceVec, CosmParam fiducial);
 void invertM(double **fisher, int nparam, double **ifisher);

 // class to contain cosmological parameters
class Cosmo {
  public:
    double wm,      // omega m
           wl,      // omega lambda
           wr,      // omega radiation
           w, wa,h;
};

 // Cosmology function
 double Omegamz(double z, double wm);
 double OmegaLambdaz(double z, double wl, double w, double wa);
 static double wzKernel(double z, void *params);
 double Eofz(double wm, double wl, double wr, double w, double wa, double z);
 double angdiamdist (double ommhh, double wl, double wr, double w, double wa, double z);
 double comovingdistance(double z, Cosmo cosmo);
 double ComovingVolume(double z, double ommhh, double wr, double wl, double w, double wa);
 double getVolume(double z, double ommhh, double wl, double wr, double w, double wa);
 static double k (double z, void * params);
 static double setVol (double z, void * params);
 //double omegaZ(double zCurrent, double wm, double wl);
 double growthFac(double z,double wm, double wl, double w);
 // Number Counts
 double dndlMJenkins(double z, double tM, My_params my_params);
 double cal_mlim(double ommhh, double wl, double wr, double w, double wa, double z);
 double fff(double x);
 double M200toM500(double M200, double overdensity, double z);
 double MhtoM200(double Mh, double overdensity, double z);
 double findlambdaC(double fr0);
 double findMyn( double z, double m, double wm);
 double findMyn_z( double z, double m, double wm);
 double findNnn(double z, double m,string filename, double wm);
 double findNnn_z(double z, double m,string filename, double wm);
 static double  SetdndM(double m, void * params);
 static double  SetN(double z, void *params);
 //double  Nofz(double wmmm, double wlll, double wrrr,   double wb, double sigmaa8, double nspecc, double zzz, double Mminn, double Mmaxx);
 //double  Nofz2(double wmmm, double wlll, double wrrr, double wbbb,  double sigmaa8, double nspecc, double zzz, double Mmm1, double Mmm2,double Bm0, double sigma0,double alpha, double beta);
 double  NofzFR(bool frflag, double wa, double w, double wmmmhh, double wlll, double wrrr, double wbbbhh, double sigmaa8, double nspecc, double zzz, double Mmm1, double Mmm2,double Bm0, double sigma0,double alpha, double beta);
 double  NFR   (bool frflag, double h, double wa, double w, double wmmm, double wlll, double wrrr, double wbbb,  double sigmaa8, double nspecc, double zzz1,double zzz2, double Mmm1, double Mmm2,double Bm0, double sigma0,double alpha, double beta);

 // Number Count Fisher Matrix
 void fisherMatrix();
 //void integrateFR();
 void prova();
 void derivatives();
 void combMatrix(string outfile);
 void reduceMatrix(string label, double FRR);

 // Power Spectrum
 double findPk( double z, double k);
 double biasPS(double z, double mass,My_params my_params);
 double PowSFR(bool frflag,double wa, double w,double wmmm, double wlll, double wbbb,double wrrr,  double sigmaa8, double nspecc, double Bm0, double sigma0, double alpha, double beta,double zzz, double KK);
double ComBiasFR(bool frflag, double mm1, double mm2, double wa, double w, double wmmmhh, double wlll, double wbbbhh,double wrrr,  double sigmaa8, double nspecc, double Bm0,double sigma0,double alpha, double beta, double zzz);
 double findbias(double z, double m,string filename, double wm);
 double findbias_z(double z, double m,string filename, double wm);

 // Power Spectrum Fisher Matrix
 void derivatives_PS();
 void fisherMatrix_PS();


 void initExperiment(int sur,string label,  string Mdynname,string FRname,string FRbias, string FRPk, double SkyCov, double zmax, double intzzz, double FRR, int der_flag);

 protected:

};

double getNFR(double z);
double setNFR(double z, My_params *c1, int iflag);
double getNofzFR(double m);
double setNofzFR(double m, My_params *c1, int iflag);
//double omegaZ(double zCurrent, double wm, double wl);

#endif
