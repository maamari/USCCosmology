
#ifndef FISHERVPAIR_H
#define FISHERVPAIR_H

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

 // class to contain cosmological parameters
class Cosmo {
  public:
    double wm,      // omega m
           wl,      // omega lambda
           wr,      // omega radiation
           w, wa,h,
          sigma8,nspec,omnu,wb;
};



class My_params {
  public:
    double z, M, r, k, moment, fr0;
    int l;
    Cosmo cosmo;
};

class My_params_cov {
  public:
    double z1,       // redshift 1 to the cluster
           z2,       // redshift 2 to the cluster
           r1,       // separation1
           r2,       // separation2
           deltar,   // r bin size 
           Mmin,       // Mass of the cluster, in units of 1e15 Msun
           fr0;     // f(R) parameter
           //wm,      // Omega matter
           //wl;      // lambda
    int l;      // multipole
    string sfr0;   // fr0 filename suffix
    Cosmo cosmo;
};

// FisherVPAIR class declaration

class FisherVPAIR : public Fisher
{
 public:
  double zmaxx,intzz, Rmin, Rmax, deltar;
  double SCCC,fr_0;
  int fisher_flag,  dFlag, survey;
  string CLL_name, nameSave,iCOV_name, COV_name;

  FisherVPAIR();
 ~FisherVPAIR();

 double getParamFromTag_cluster(string tag, CosmParam *c);
 void specifyParams_cluster(int *choiceVec, CosmParam fiducial);
 void invertM(double **fisher, int nparam, double **ifisher);


//Cosmology
static double Growth(double z, void *params);
double Dv(double z, Cosmo cosmo);
double Eofz(double z, Cosmo cosmo);


//vpair theory
double theory_vpair_ferreira(double r, double z, double fr0, double minmass,double w1,double w0,double wm,double wl,double wb,double sigma8, double nscal);
double xi_halo_bar(My_params my_params);
static double xi_halo_bar_kernelr(double r, void *params);
double xi_halo(My_params my_params);
static double xi_halo_kernel(double k, void *params);
double effective_b(My_params my_params, double Mmin);
//covariance matrix
void calccov();
void cosmic_variance(Cosmo cosmo,  double z, double Mmin);
double Wx(double x);
static double cv_kernel(double k, void *params);



 //  Fisher Matrix
 void fisherMatrix();
 void derivatives();
 void combMatrix(string outfile);
 void reduceMatrix(string label, double FRR);
 

 void initExperiment(int sur,string label, double SkyCov, double zmax, double intzzz, double rmin, double rmax, double intrrr, int der_flag);

 protected:

};

#endif
