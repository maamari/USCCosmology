/* fisherBAO.cc
 * Contains functions for calculating fisher matrix for a Galaxy experiment
  */

#include <math.h>
#include <iostream>
#include <fstream>
#include "fisherBAO.h"
#include "fisher.h"
#include "dcosmology2.h"
#include "dnumrecipes.h"
#include "spline.h"

using namespace std;

//References:
//Spitzer:  Spitzer, "Physical Processes in the ISM".


/*********************************************************************
 **************** Constructor/Destructor *****************************
 ********************************************************************/

 FisherBAO::FisherBAO()
{
	//cout<<"FisherBAO constructor called"<<endl;
	fisher_flag=0;
       
}

 FisherBAO::~FisherBAO()
{
	//cout<<"FisherBAO destructor called"<<endl;

}

//////////////////////////////////////////////////////////////////////
// Initialisation routine
//////////////////////////////////////////////////////////////////////

//initialise Galaxy survey parameters
void FisherBAO::initSurvey(string name1, double obsarea1, double zmin1, double zmax1, double sigmaz1)
{
	exp_name=name1;
        obsarea=obsarea1;
	zmin=zmin1;	
	zmax=zmax1;
	exp_z=(zmin+zmax)/2.0;
	sigmaz=sigmaz1;
        ncal=1000.0;

	//volume=volumeFromArea(zmin,zmax,obsarea);

        double cH;
        cH=SPEEDOFLIGHT_CGS/H0_CMSMPC;
        volume=pow(coAngDiamDistance(exp_z,&fiducialFC),2.0)/hubbleZ(exp_z,&fiducialFC);
        volume*=(zmax-zmin)*(obsarea/4.0e4);
        volume*=cH*cH*cH;


	sigmar=0.0;
	dirbase=dirbase+"_z"+numToString(exp_z);


	cout<<"bao survey: "<<exp_name<<" initialised"<<endl;
}

//////////////////////////////////////////////////////////////////////
// Member functions
//////////////////////////////////////////////////////////////////////

string FisherBAO::getSurveyName()
{
	return exp_name;
}

double FisherBAO::getSurveyZ()
{
	return survey_z;
}

double FisherBAO::getZMin()
{
	return zmin;
}

double FisherBAO::getZMax()
{
	return zmax;
}

double FisherBAO::getVolume()
{
	return volume;
}

void FisherBAO::setCosmology(CosmParam fiducial)
{
	fiducialFC=fiducial;
}

////////////////////////////////////////////////////////////////////////
// Calculate Fisher matrix
////////////////////////////////////////////////////////////////////////

void FisherBAO::fisherMatrixBAO()
{
   double sigmaD, sigmaH, sigmaDZ;
   double dDdDZ, dHdDZ, dz, cH;
  double vol0(2.16e9);  //volume in (h^{-1}Mpc)^3 
  double x0d(0.0085), x0h(0.0148);
  double **fisher, **ifisher;
  int nmat;
  string tag;
  int deltaz_flag(0);

  vol0/=pow(fiducialFC.ommhh/(1.0-fiducialFC.oml-fiducialFC.omk),1.5);

  //establish parameters
  param_tags.clear();
  param_values.clear();
  param_tags.push_back("_fiducial");
  param_values.push_back(0.0);

  fiducialFC.lda=log(coAngDiamDistance(exp_z,&fiducialFC));
  fiducialFC.lh=log(hubbleZ(exp_z,&fiducialFC));

  nparamFC=0;
  addParam("_lh");
  addParam("_lda");

  if(sigmaz>1.0e-5 && deltaz_flag) addParam("_deltaz");

  nmat=nparamFC;
  fisher=dmatrix(1,nmat,1,nmat);
  ifisher=dmatrix(1,nmat,1,nmat);

  cout<<nmat<<endl;

  saveCosmParam(&fiducialFC,fiducialFC.name);

  //photo-z performance
  if(sigmaz>1.0e-5 && deltaz_flag){
     cH=SPEEDOFLIGHT_CGS/H0_CMSMPC;
     x0d=0.0123*sqrt(sigmaz*(1.0+exp_z));
     x0d/=sqrt(34.1/fiducialFC.h*exp(fiducialFC.lh)/cH);
     x0h=1.0e40;  //numerically infinity
     sigmaDZ=sigmaz*(1.0+exp_z)/sqrt(ncal);
  }

  sigmaD=x0d*4.0/3.0*sqrt(vol0/volume)*getFnl(exp_z);
  sigmaH=x0h*4.0/3.0*sqrt(vol0/volume)*getFnl(exp_z);

  sigmaH=sqrt(sigmaH*sigmaH+pow(0.01*sqrt(0.5/(zmax-zmin)),2.0));
  sigmaD=sqrt(sigmaD*sigmaD+pow(0.01*sqrt(0.5/(zmax-zmin)),2.0));

  fisher[1][1]=1.0/sigmaH/sigmaH;
  fisher[2][2]=1.0/sigmaD/sigmaD;
  fisher[1][2]=0.0;
  fisher[2][1]=0.0;

  if(sigmaz>1.0e-5 && deltaz_flag){
     dz=0.01;
     dDdDZ=(log(coAngDiamDistance(exp_z+dz,&fiducialFC))-log(coAngDiamDistance(exp_z-dz,&fiducialFC)))/(2.0*dz);
     dHdDZ=(log(hubbleZ(exp_z+dz,&fiducialFC))-log(hubbleZ(exp_z-dz,&fiducialFC)))/(2.0*dz);
     cout<<exp_z<<"\t"<<dDdDZ<<"\t"<<dHdDZ<<endl;
     fisher[1][3]=1.0/sigmaH/sigmaH*dHdDZ;
     fisher[3][1]=fisher[1][3];
     fisher[2][3]=1.0/sigmaD/sigmaD*dDdDZ;
     fisher[3][2]=fisher[2][3];
     fisher[3][3]=pow(1.0/sigmaH*dHdDZ,2.0)+pow(1.0/sigmaD*dDdDZ,2.0);
     //add delta z prior
     fisher[3][3]+=1.0/sigmaDZ/sigmaDZ;
  }

  //save information about fisher matrix and inverse into class
  setFisherMatrix(fisher,nmat,&fiducialFC);

  //write result to file
  tag="_z"+numToString(exp_z)+".dat";
  saveFisher(tag);
  free_dmatrix(fisher,1,nmat,1,nmat);
  free_dmatrix(ifisher,1,nmat,1,nmat);

  fisher_flag=1;

  return;
}

////////////////////////////////////////////////////////////////////////
// Non-linear smoothing 
////////////////////////////////////////////////////////////////////////
double FisherBAO::getFnl(double z)
{
	double zm(1.4);
	double gamma(0.5);

	if(z>zm) return 1.0;

	return pow(zm/z,gamma);
}

////////////////////////////////////////////////////////////////////////
// Add systematic errors to distance measurements
////////////////////////////////////////////////////////////////////////
void FisherBAO::addSystematicErrors(double deltaz)
{
	int i;	
	string tag;
	double **ifishernew, **fishernew;
	double sigmaH, sigmaD;

	sigmaD=0.01*sqrt(0.5/deltaz);
	sigmaH=sigmaD;

	ifishernew=dmatrix(1,nparamFC,1,nparamFC);
	fishernew=dmatrix(1,nparamFC,1,nparamFC);

	//copy fisherFC into ifishernew
	copyMatrix(ifisherFC,nparamFC,ifishernew);

	//add systematic error in quadrature to inverse fisher matrix
	for(i=1;i<=param_tags.size()-1;i++){
		 tag=param_tags[i];
		   if(tag.find("_lh",0)!=string::npos){
			 ifishernew[i][i]=sigmaH*sigmaH+ifisherFC[i][i];
		   }
		   if(tag.find("_lda",0)!=string::npos){
			 ifishernew[i][i]=sigmaD*sigmaD+ifisherFC[i][i];
		   }
	}

	//now invert inverse fisher matrix to get new fisher matrix w/sys err
	invertMatrix(ifishernew,nparamFC,fishernew);	

	setFisherMatrix(fishernew,nparamFC,&fiducialFC);

	free_dmatrix(fishernew,1,nparamFC,1,nparamFC);
	free_dmatrix(ifishernew,1,nparamFC,1,nparamFC);
}

////////////////////////////////////////////////////////////////////////
// Add redshift bin
////////////////////////////////////////////////////////////////////////
void FisherBAO::addRedshiftBin(double zmin1, double zmax1)
{
	fisherData oldData, newData;

	//save old fisherData somewhere safe
	oldData=returnFisherData();
	saveFisher("_outs");
	oldData.fisher_file=fisher_fileFC;

	//update redshift information
	zmin=zmin1;	
	zmax=zmax1;
	exp_z=(zmin+zmax)/2.0;
	volume=volumeFromArea(zmin,zmax,obsarea);

	//calculate new redshift bin data
	fisherMatrixBAO();
	newData=returnFisherData();
	
	combineFisherMatrices(oldData,newData,&fiducialFC);

}

void FisherBAO::addDistanceMeasurement(double zbin, double sigmaD, double sigmaH)
{
	fisherData oldData, newData;
        double **fisher;
        int nmat(2);
        string tag;
        int old_data_flag(0);

        cout<<fisher_flag<<endl;

	//save old fisherData somewhere safe
        if(fisher_flag==1){
           oldData=returnFisherData();
           saveFisher("_outs");
           oldData.fisher_file=fisher_fileFC;
           old_data_flag=1;
        }

	//update redshift information
	exp_z=zbin;

	//calculate new redshift bin data
	//fisherMatrixBAO();

        //establish parameters
        param_tags.clear();
        param_values.clear();
        param_tags.push_back("_fiducial");
        param_values.push_back(0.0);

        fiducialFC.lda=log(coAngDiamDistance(exp_z,&fiducialFC));
        fiducialFC.lh=log(hubbleZ(exp_z,&fiducialFC));

        nparamFC=0;
        addParam("_lh");
        addParam("_lda");

        nmat=nparamFC;
        fisher=dmatrix(1,nmat,1,nmat);
 
        cout<<"no. param: "<<nmat<<endl;

        saveCosmParam(&fiducialFC,fiducialFC.name);

        fisher[1][1]=1.0/sigmaH/sigmaH;
        fisher[2][2]=1.0/sigmaD/sigmaD;
        fisher[1][2]=0.0;
        fisher[2][1]=0.0;

        //save information about fisher matrix and inverse into class
        setFisherMatrix(fisher,nmat,&fiducialFC);

        //write result to file
        tag="_z"+numToString(exp_z)+".dat";
        cout<<tag<<endl;
        saveFisher(tag);
        free_dmatrix(fisher,1,nmat,1,nmat);

        //combine into new fisher matrix
	newData=returnFisherData();
        if(old_data_flag==1){
           combineFisherMatrices(oldData,newData,&fiducialFC);
        }

        fisher_flag=1;
}
