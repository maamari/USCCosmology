//fisher.h
//
// Contains basic structure for accessing CAMB and manioulating fisher
// matrices

#ifndef FISHER_H
#define FISHER_H

#include "spline.h"
#include <string>
#include <vector>
#include "dcosmology2.h"
#include "dnumrecipes.h"
#include <cstring>
using namespace std;

const int NPRECISION(14); //sigfig for file io
const double SIGMATHOMSON(6.65246e-25); //Thomson scattering x-section (cgs)
//const int transfer_spline_flag(0);

struct fisherData {
	double exp_z;
	int nparam;
	string fisher_file;
	vector<string> param_tags;
	vector<double> param_values;
};

// Cosmological parameter structure
struct CosmParam {
  // Cosmological parameters
  double h;
  double omb;
  double omm;
  double omk;
  double oml;
  double omnu;
  double tau;
  double ts;
  double nscal;      // scalar spectral index
  double alpha;      // running of spectral index
  double ascal2;     // scalar amplitude
  double yHe;        // primordial helium fraction
  double sigma8;
  double bias;       // multiplicative factor for power spectrum

  int physical_flag;

  double ombhh;
  double ommhh;
  double omnuhh;

  double w0;         //DE parameters
  double w1;

  double exp_z;  // redshift that experiment observes

  string name;
  string file;
  string cl_file;
  string clustern_file;
  string filebase;

  //derivative step
  double p_step;

  //reference cosmology factors
  double daratio;
  double hratio;

  //distances
  double lda;    //log of the angular diameter distance
  double lh;	 //log of H(z)
  double lg;	 //log of G(z) the growth factor

  //neutrinos
  int hierarchy_flag;
  double numass1;
  double numass2;
  double numass3;

  //inflation
  double ntens;

  //reionization
      double zri;
      double delta_zri;

  // Cluster mass observable relations
  double alpha_cluster;
  double beta;
  double Bm0;
  double sigma0;

  // f(R) field amplitude
  double fR0;
  double currentZ;
  int frflag;
  double m1, m2;
};

// Fisher class declaration

const int nparamCORE(6);

class Fisher {

 public:
Fisher();
~Fisher();

//member functions
void setFisherMatrix(double **fisher, int nparam, CosmParam *fiducial);
CosmParam setCosmParam(double h, double omm, double omb, double oml, double omnu, double tau, double ascal2, double nscal, double alpha, double ts, double w0, double w1, double sigma8,string name);

void setName(string name);
void setExpZ(double z1);
int setHierarchyFlag(int hflag);
int setSlowrollFlag(int srflag);
int setLensingFlag(int lflag);
int setAltNeutrinoParam();
void assignSpline(Spline *spline, string name, string tag);
void assignSpline2(Spline *spline, string name, string tag, int flag);
void setDirbase(string name);

double getParamFromTag(string tag, CosmParam *c);
double getErrorFromTag(string tag);
double getFisherFromTag(string tag);


int getFisherSize();
int getFisherMatrix(double **fisher);
int getInverseFisherMatrix(double **ifisher);
vector<string> returnParam_Tags();
vector<double> returnParam_Values();
fisherData returnFisherData();
string returnDirbase();
string numToString(double z);
string numToStringLong(double z,int np);

void copyMatrix(double **A, int nmat, double **B);

//cosmology functions
double getZCMB(CosmParam *fiducial);
double angDiamDistance(double zCall, CosmParam *c);
double coAngDiamDistance(double zCall, CosmParam *c);
double lumDistance(double zCall, CosmParam *c);
double coordDistanceNum(double zCall, CosmParam *c);
double hubbleZ(double z, CosmParam *c);
double omegaMZ(double z, CosmParam *c);
double growthZ(double z, CosmParam *c);
double getDVolumeDOmega(double z, CosmParam *c);
double volumeFromArea(double zmin, double zmax, double area);
double getGrowth(double ause,CosmParam *c);

//reionization for CAMB
double getXE(double z, CosmParam *c);
double getXH(double z, CosmParam *c);
double getTauCMB(CosmParam *c);
void setZRI(double zri, double delta_zri);

//File I/O
void saveFisher(string tag);
void loadFisher(string name);
void loadFisherFile(string name, double **fisher, int nmat);
void loadParamValues(string file);
void specifyParams(int *choiceVec, CosmParam fiducial);
void addParam(string tag);
void addParamWithValue(string tag, double value);
void removeParam(string tag);
void setCosmology(CosmParam fiducial);

//translation
double omnuhhFromMnu(double Mnu);
double mnuFromOmnuhh(double omnuhh);
void CosmParamFromCosmology(CosmParam cosm, Cosmology c, string name);
double daratio(CosmParam *c);
double hratio(CosmParam *c);

//cosmological parameter initialisation routines
double saveCosmParam(CosmParam *c, string file);
CosmParam loadCosmParam(string file);
vector<double> callCAMB(CosmParam *c);
double powerSpectrum(double k, CosmParam *c, Spline *spline);
double powerSpectrumFid(double k, double z);

//control calls to CAMB and setting cosmologies
double setCosmFiles();
double setParamPair(CosmParam *fiducial, int iflag, int clflag);
void derivCMB(int lmax, double p_step,CosmParam *cm, CosmParam *cp,string file);
void derivGAL(string tag,string file);

//derivatives w.r.t fixed d_A
void derivCMB_OMK(int lmax, CosmParam *c, string name);
void derivCMB_W0(int lmax, CosmParam *c, string name);
double fixDA(CosmParam *c, double thetas);
double fixThetas(CosmParam *c, double thetas);

//useful scales
double angularScale(CosmParam *c);
double soundHorizon(CosmParam *c);

//Error ellipses
int fisherContour();
int fisherContourNotMarg();

//reduce size of fisher matrix
void reduceFisherMatrix(int ichoice[]);
void combineFisherMatrices(fisherData fisherData1, fisherData fisherData2, CosmParam *fiducial);
void copyFisherMatrix(fisherData fisherData1, CosmParam *fiducial);

//out put errors
void printErrors();

//latex tables
void latexTable();
string getLatexFromTag(string tag);

// cross-correlation
void normaliseFisher();

//dark energy
void projectDarkEnergy(int de_flag);
double deDerivatives(string tagn, string tagi);
double figureOfMerit(vector<string> params);
double marginalisedErrors(vector<string> params);

//projection onto wp-wa
void projectDarkEnergyWP();
double deDerivativesWP(string tagn1, string tagi, double ap1);

//slow roll parameters
void projectInflation(int de_flag);
double infDerivatives(string tagn, string tagi);

//neutrino parameters
void projectNeutrinoMasses(int de_flag);
double mnuDerivatives(string tagn, string tagi);

//omegac
void projectOmegaC();
double omegacDerivatives(string tagn1, string tagi);

//omegadm
void projectOmegaDM();
double omegadmDerivatives(string tagn1, string tagi);

//f_nu
void projectFNu();
double fnuDerivatives(string tagn1, string tagi);

//systematic errors
void addSystematicErrors(double sigmaD, double sigmaH);
void addSystematicErrorFromTag(string param, double sigma);
void addPriorFromTag(string param, double sigma);


 protected:

  //fiducial cosmology for fisher matrix
  CosmParam fiducialFC;

  //fisher matrix details
  string exp_name;
  string z_tag;
  double exp_z;
  int fisher_flag;
  int nparamFC;
  double **fisherFC;
  double **ifisherFC;

  string dirbase;

  //store parameter names
  int nparamFC_file;
  string fisher_fileFC;  //name of file where fisher is written
  vector<string> param_tags;
  vector<double> param_values;
  vector<double> redshift_list;

  int hierarchy_flag;  //neutrino hierarchy
  int neutrino_mass_flag;
  int neutrino_mass_flag_alt;
  int transfer_spline_flag;

  //distance fisher matrix
  int distance_flag;   //treat distances as independent parameters
  int slowroll_flag;   //use slow roll inflationary parameters and CAMB
  int lensing_flag;    //apply weak gravitational lensing to CMB
  int omnu_flag;       //use omnu instead of omnuhh - stupid thing to do
  int omk_flag;        //handle derivatives w.r.t omk carefully
  int reionization_flag;  //use tanh reionization model
};

//misc utilities
char *concat(char *str1, char *str2);
void printMatrix(double **matrixin, int nmat);

//distance integration functions
double setCoordDistanceKernelF(double z, CosmParam *c1, int iflag);
double getCoordDistanceKernelF(double z);

//growth function integration code
double setGrowth(double a,CosmParam *c1,int iflag);
double intGrowth(double a);

//fix angular diameter distance
double getFixDA(double oml);
double setFixDA(double oml, CosmParam *c1, Fisher *ff1, int iflag);

double getSoundHorizonKernel(double x);
double setSoundHorizonKernel(double x, CosmParam *c1, Fisher *ff1, int iflag);
double setFixThetas(double oml, CosmParam *c1, Fisher *ff1, int iflag);
double getFixThetas(double oml);

//tau CMB
double setTauIntF(double z, CosmParam *c1, Fisher *ff1,int iflag);
double getTauIntF(double z);

//matrix manipulation
void invertMatrix(double **fisher,int nparam, double **ifisher);
void gaussj(double **a, int n, double **b, int m);
void writeMatrix(double **M,int nparam, char *file);
void writeMatrixS(double **M,int nparam, string file);
double detMatrix(double **M, int nparam);
void multMatrix(double **A, double **B, double **C, int N);
void invertMatrix2(double **fisher,int nparam, double **ifisher);
void invertMatrixSVD(double **fisher,int nparam, double **ifisher);

#endif
