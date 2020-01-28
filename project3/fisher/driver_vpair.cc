/* driver.cc
 *
 * Program for computing the fisher matrix for a galaxy survey.
 * Call by
 *   ./driver.x -z6 > (output filename)
 */

#include <iostream>
#include <fstream>
#include "dnumrecipes.h"
#include "math.h"
#include "dcosmology2.h"
#include "fisherCMB.h"
#include "fisherGAL.h"
#include "fisherVPAIR.h"
#include "fisher.h"


int bubble_flag(0);
int bubble_use(0);

using namespace std;

const double INT_ACCURACY(1.0e-6);

int main(int argc, char *argv[])
{
    int i, survey;
    double frvalue;
    string frname;
  int reset_flag(0);
  //char *file;
  ofstream fout;
  /////////////////////////////////////////////////////////////
  // Handle command line input
  ////////////////////////////////////////////////////////////

  if ((argc ==4) && (argv[1][0]=='-')) {
    switch (argv[1][1]) {
    case 'r':
      reset_flag = atoi(&argv[1][2]);
      break;
    default:
      cerr << "Bad option " << argv[1] << "\n";
    }
    //--argc;
    //++argv;
    survey = (int)atof(argv[2]);
    frvalue = (double)atof(argv[3]);
    frname = argv[3] ;
  } else if (argc == 3 && argv[1][0] != '-') {
    survey = (int)atof(argv[1]);
    frvalue = (double)atof(argv[2]);
    frname = argv[2] ;
  }
  else {
    cerr << endl;
    cerr << "Usage: ./driver.x <survey> <frvalue> " << endl;
    cerr << "<survery> 0 for Planck; " << endl;
    cerr << "          1 for SPT; " << endl;
    cerr << "          2 for ACT; " << endl;
    cerr << "<frvalue> fR0 value; " << endl;
    cerr << endl;
    return 0;
  }

  /////////////////////////////////////////////////////////////////
  // Specify fiducial cosmology
  /////////////////////////////////////////////////////////////////
  FisherCMB WMAPF;
  FisherCMB PLANCK, PLANCKT;
  string name;

  CosmParam fiducial;
  double omm,omb,oml,omnu,tau,ascal2,nscal,alpha,ts,w0,w1, sigma8,h;
  double Mnu;
  string namec="_fiducial", files;
  string dirbase;

  //dirbase="./checkdir";
  dirbase="./data";
  files="mkdir "+dirbase;
  system(files.c_str());

  //Fabian Fiducial Values
  omm=0.27;
  omb=0.046;
  oml=1.0-omm;
  Mnu=0.3;
  h=0.72;
  //omnu=0.2/h/h;
  cout << WMAPF.omnuhhFromMnu(Mnu) << endl;
  omnu=WMAPF.omnuhhFromMnu(Mnu)/h/h;
  tau=0.1;
  //nscal=1;
  nscal=0.963;
  alpha=0.0;
  ascal2=29.5;
  ascal2*=0.9;
  ts=0.0;
  w0=-1.0;
  w1=0.0;
  //sigma8=0.9;
  sigma8=0.809;
  if(reset_flag==1){
    fiducial=WMAPF.setCosmParam(h,omm,omb,oml,omnu,tau,ascal2,nscal,alpha,ts,w0,w1,-1.0,namec);
  }else{
    fiducial=WMAPF.setCosmParam(h,omm,omb,oml,omnu,tau,ascal2,nscal,alpha,ts,w0,w1,sigma8,namec);
  }

  //Cosmology cosm(omm,oml,omb,h,fiducial.sigma8,nscal,omnu,w0,w1);

// specify parameters for reset
  int *choiceVec;
  choiceVec=ivector(1,23);
  for(i=1;i<=23;i++) choiceVec[i]=1;
  choiceVec[1]=1;   //ommhh
  choiceVec[2]=1;   //ombhh
  choiceVec[3]=1;   //oml
  choiceVec[4]=1;   //nscal
  choiceVec[5]=1;   //ascal2
  choiceVec[6]=1;   //cosmological constant
  choiceVec[7]=0;   // tau
  choiceVec[8]=0;   //tensors
  choiceVec[9]=0;   //flat Universe
  choiceVec[10]=0;   // helium
  choiceVec[11]=0;  //running
  choiceVec[12]=1;  //Mnu
  choiceVec[13]=0;  // bias
  choiceVec[14]=0;  //numass1
  choiceVec[15]=0;  //numass2
  choiceVec[16]=0;  //numass3
  choiceVec[17]=0;  //no lda
  choiceVec[18]=0;  //no lh
  choiceVec[19]=0;  //no lg
  choiceVec[20]=0;  //epsilon
  choiceVec[21]=0;  //eta
  choiceVec[22]=0;  //xi
  choiceVec[23]=0;  //w1

  WMAPF.specifyParams(choiceVec,fiducial);

  cout<<"initialising experiments"<<endl;
  //initialise different cosmologies for derivatives
  if(reset_flag==1) WMAPF.setCosmFiles();

  /////////////////////////////////////////////////////////////////
  // CMB Experiment Parameters
  ////////////////////////////////////////////////////////////////
  choiceVec[1]=1;   //ommhh
  choiceVec[2]=1;   //ombhh
  choiceVec[3]=1;   //oml
  choiceVec[4]=1;   //nscal
  choiceVec[5]=1;   //ascal2
  choiceVec[6]=1;   //cosmological constant
  choiceVec[7]=0;   // tau
  choiceVec[8]=0;   //tensors
  choiceVec[9]=0;   //flat Universe
  choiceVec[10]=0;   // helium
  choiceVec[11]=0;  //running
  choiceVec[12]=0;  //Mnu
  choiceVec[14]=0;  //numass1
  choiceVec[15]=0;  //numass2
  choiceVec[16]=0;  //numass3
  choiceVec[16]=0;  //numass3
  choiceVec[17]=0;  //no lda
  choiceVec[18]=0;  //no lh
  choiceVec[19]=0;  //no lg
  choiceVec[20]=0;  //epsilon
  choiceVec[21]=0;  //eta
  choiceVec[22]=0;  //xi
  choiceVec[23]=1;  //w1
  ////////////////////////////////////////////////////////////
  int lmax(1000);
  //int lmax(2500);

  cout<<"calculating CMB Fisher matrix"<<endl;
  //FisherCMB WMAP;
  //WMAP.setDirbase(dirbase);
  //WMAP.specifyParams(choiceVec,fiducial);
  //name=dirbase+"/wmap";
  //WMAP.initExperiment(name,0.8,3,lmax,26,0.68*60,4.3,6.1);  //Ka
  //WMAP.addChannelExperiment(33,0.53*60.0,7.1,10.0);   //Q
  //WMAP.addChannelExperiment(49,0.35*60.0,16.3,23.0);   //V

  //WMAP.fisherMatrixCMB();

  //WMAP.printErrors();

  //WMAP.latexTable();

  //WMAP8 Ka+Q+V - channels used for analysis - only Q+V used in WMAP3 (Dunkley 2008)
  /*FisherCMB WMAP8;
  WMAP8.setDirbase(dirbase);
  WMAP8.specifyParams(choiceVec,fiducial);
  name=dirbase+"/wmap8";

  WMAP8.initExperiment(name,0.8,3,lmax,26,0.68*60,3.4,4.8);  //Ka
  WMAP8.addChannelExperiment(33,0.53*60.0,5.6,7.9);   //Q
  WMAP8.addChannelExperiment(49,0.35*60.0,12.9,18.2);   //V

  WMAP8.fisherMatrixCMB();

  WMAP8.printErrors();
  WMAP8.latexTable();*/


  PLANCK.setDirbase(dirbase);
  PLANCK.specifyParams(choiceVec,fiducial);
  name=dirbase+"/planck_wmnu0.1_wlensing_l1000";
  PLANCK.initExperiment(name,0.65,3,lmax,143,7.1,6.0,11.4);
  PLANCK.addChannelExperiment(100,9.5,6.8,10.9);
  PLANCK.addChannelExperiment(217,5.0,13.1,26.7);
  PLANCK.addChannelExperiment(353,5.0,40.1,10.9);
  PLANCK.fisherMatrixCMB_DE(-1);
  PLANCKT.copyFisherMatrix(PLANCK.returnFisherData(),&fiducial);

  PLANCK.projectAngularScaleOntoDE(3);

  PLANCK.fisherMatrixCMB();

  PLANCKT.printErrors();

  PLANCK.latexTable();


  /////////////////////////////////////////////////////////////////
  // BAO Experiment Parameters
  ////////////////////////////////////////////////////////////////
  /*cout<<"calculating BAO Fisher matrix"<<endl;
  FisherBAO BAO;
  double skycov, zmin, zmax, sigmaz;
  string label;
  // specify the Cluster survey: Planck, SPT, ACT
  name=dirbase+"/bao";
  if (survey ==0) {// Planck
      skycov=fsdegree*0.75;
      label="Planck";
  } else if (survey ==1) { // SPT
      skycov=fsdegree*0.097;
      label="SPT";
  } else if (survey ==2) { // ACT
      skycov=fsdegree*0.087;
      label="ACT";
  }
  sigmaz=0.01; // 0.01 optimistic 0.05 pessimistic 0 spectroscopic
  zmin=0;
  zmax=0.15;
  BAO.setCosmology(fiducial)
  BAO.initSurvey(name,skycov,zmin,zmax,sigmaz);
  BAO.fisherMatrixBAO();
  BAO.addRedshiftBin(0.15,0.3);
  BAO.addRedshiftBin(0.30,0.45);
  BAO.addRedshiftBin(0.45,0.60);
  BAO.addRedshiftBin(0.60,0.75);
  BAO.addRedshiftBin(0.75,0.90);
  BAO.addRedshiftBin(0.90,1.0);
  BAO.printErrors();*/

  /////////////////////////////////////////////////////////////////
  // Galaxy Experiment Parameters
  ////////////////////////////////////////////////////////////////
  choiceVec[7]=0;   // tau
  choiceVec[8]=0;   //tensors
  choiceVec[10]=0;   //helium
  choiceVec[11]=0;  //no running
  choiceVec[13]=1;  //bias
  ////////////////////////
  /*FisherGAL SDSS;
  double volume, density, kmin, kmax, sigma8gal,zuse, sigmaz;

  sigmaz=0.0;

  name=dirbase+"/sdss";
  volume=1e9*pow(h,-3.0);
  density=1e-4*pow(h,3.0);
  kmin=0.0003;
  kmax=0.1*h;
  sigma8gal=1.8;
  zuse=0.3;
  SDSS.initSurvey(name, zuse, volume, density, kmin, kmax, sigma8gal, sigmaz);

	SDSS.specifyParams(choiceVec,fiducial);

	 SDSS.addParam("_bias");
	 SDSS.removeParam("_tau");
	 SDSS.removeParam("_ts");
	SDSS.fisherMatrixGAL();

  SDSS.printErrors();
  SDSS.latexTable();
  SDSS.outPowerErrors(); */

    /*  Fisher SDSS_PLANCK;
  name=dirbase+"/sdss_planck";
  SDSS_PLANCK.setName(name);

  SDSS_PLANCK.combineFisherMatrices(SDSS.returnFisherData(), PLANCK.returnFisherData(), &fiducial);
  SDSS_PLANCK.latexTable(); */


    /*FisherGAL LSST;

  name=dirbase+"/lsst";
  volume=1.8e9*pow(fiducial.h,-3.0);
  density=5.0e-4*pow(fiducial.h,3.0);
  kmin=0.0003;
  kmax=0.2*fiducial.h;
  sigma8gal=1.0;
  zuse=1.0;
  LSST.initSurvey(name, zuse, volume, density, kmin, kmax, sigma8gal, sigmaz);

	LSST.specifyParams(choiceVec,fiducial);

	 LSST.addParam("_bias");
	 LSST.removeParam("_tau");
	 LSST.removeParam("_ts");
	LSST.fisherMatrixGAL();

  LSST.printErrors();
  LSST.latexTable();
  LSST.outPowerErrors(); */


  //Fisher LSST_PLANCK;
  //name=dirbase+"/lsst_planck";
  //LSST_PLANCK.setName(name);

  //LSST_PLANCK.combineFisherMatrices(LSST.returnFisherData(), PLANCK.returnFisherData(), &fiducial);
  //LSST_PLANCK.latexTable();
  //LSST_PLANCK.printErrors();



  /////////////////////////////////////////////////////////////////
  // Galaxy Cluster Experiment Parameters
  ////////////////////////////////////////////////////////////////
  FisherVPAIR VPAIR;
  double skycov=0;
  string label;
  if (survey ==0) {// Planck
      skycov=fsdegree*0.75*pow(PI/180.,2);
      label="Planck";
  } else if (survey ==1) { // SPT
      skycov=2500*pow(PI/180.,2);
      label="SPT";
  } else if (survey ==2) { // ACTpol
      skycov=4000*pow(PI/180.,2);
      label="ACT";
  } else if (survey ==3) { // SPTpol
      skycov=625*pow(PI/180.,2);
      label="SPTpol";
  } else if (survey ==4) { // CCAT
      skycov=2000*pow(PI/180.,2);
      label="CCAT";
  } else if (survey ==6) { // CCATwide
      skycov=33000*pow(PI/180.,2);
      label="CCATwide";
  } else if (survey ==7) {
      label="CCATwide_10k";
      skycov=10000.*pow(PI/180.,2);
  } else if (survey ==8) {
      label="CCATwide_5k";
      skycov=5000.*pow(PI/180.,2);
  } else if (survey ==9) {
      label="CCATwide_1k";
      skycov=1000.*pow(PI/180.,2);
  }

  // cout<<endl << "calculating the VPAIR Fisher Matrix" << endl;
  // VPAIR.specifyParams(choiceVec,fiducial);
  //
  // VPAIR.initExperiment(survey,label,skycov , 1,0.1 ,20,100,2, 1);
  //
  // cout<<endl << "calculating the VPAIR covariance" << endl;
  // VPAIR.calccov();
  // cout<<endl << "Finished calculating the VPAIR Fisher Matrix" << endl;
  // cout<<endl << "calculating the VPAIR Derivatives" << endl;
  // VPAIR.derivatives();
  // VPAIR.fisherMatrix();

////////////////////////////////////////////////////////////////
  //clean up memory and end
  return 0;
}
