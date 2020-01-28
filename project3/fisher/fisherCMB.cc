/* fisherCMB.cc
 * Contains functions for calculating fisher matrix for a CMB experiment
  */

#include <math.h>
#include <iostream>
#include <fstream>
#include "fisherCMB.h"
#include "fisher.h"
#include "dcosmology2.h"
#include <gsl/gsl_sf_coupling.h>
#include <gsl/gsl_sf_result.h>
#include "dnumrecipes.h"
#include "spline.h"

using namespace std;

//References:
//Spitzer:  Spitzer, "Physical Processes in the ISM".


/*********************************************************************
 **************** Constructor/Destructor *****************************
 ********************************************************************/

 FisherCMB::FisherCMB()
{
	//cout<<"FisherCMB constructor called"<<endl;
	nchannel=0;
	fisher_flag=0;
        nparamFC=nparamCORE+nparamCMBX;
        lminPol=2;

}

 FisherCMB::~FisherCMB()
{
	//cout<<"FisherCMB destructor called"<<endl;

	if(nchannel>0){
		free_dvector(beam,1,nchannel);
		free_dvector(sigmat,1,nchannel);
		free_dvector(sigmap,1,nchannel);
	}
}

//////////////////////////////////////////////////////////////////////
// Initialisation routine
//////////////////////////////////////////////////////////////////////

//initialise CMB experiment parameters
void FisherCMB::initExperiment(string name1, double fsky1, int xmax1, int lmax1, double nu1, double beam1, double sigmat1, double sigmap1)
{

	exp_name=name1;
	fsky=fsky1;
	xmax=xmax1;
	lmax=lmax1;
	addChannelExperiment(nu1,beam1,sigmat1,sigmap1);

}

void FisherCMB::addChannelExperiment(double nu1, double beamin, double sigmat1, double sigmap1)
{
	double **beamsave;
        double beam1,arcmin(2.90888e-4); //conversion between arcmin and radians
	int i;

	beam1=beamin*arcmin;

	//increment channel counter
	nchannel++;

	if(nchannel>1){
		//store old beam data
		beamsave=dmatrix(1,nchannel-1,1,3);
		for(i=1;i<=nchannel-1;i++){
			beamsave[i][1]=beam[i];
			beamsave[i][2]=sigmat[i];
			beamsave[i][3]=sigmap[i];
		}
		//free up old beam vectors
		free_dvector(beam,1,nchannel-1);
		free_dvector(sigmat,1,nchannel-1);
		free_dvector(sigmap,1,nchannel-1);

		//resize and enter data
		beam=dvector(1,nchannel);
		sigmat=dvector(1,nchannel);
		sigmap=dvector(1,nchannel);

		//enter old data
		for(i=1;i<=nchannel-1;i++){
			beam[i]=beamsave[i][1];
			sigmat[i]=beamsave[i][2];
			sigmap[i]=beamsave[i][3];
		}

		//enter new data
		beam[nchannel]=beam1;
		sigmat[nchannel]=sigmat1;
		sigmap[nchannel]=sigmap1;

		//clean up
		free_dmatrix(beamsave,1,nchannel-1,1,3);
	}else{
		beam=dvector(1,nchannel);
		sigmat=dvector(1,nchannel);
		sigmap=dvector(1,nchannel);
		beam[nchannel]=beam1;
		sigmat[nchannel]=sigmat1;
		sigmap[nchannel]=sigmap1;
	}

}

// change info for a channel
int FisherCMB::replaceChannel(int nreplace, double nu1, double beam1, double sigmat1, double sigmap1)
{
	if(nreplace>nchannel){
		cout<<"optimist!  Not enough channels"<<endl;
		return 0;
	}

	beam[nreplace]=beam1;
	sigmat[nreplace]=sigmat1;
	sigmap[nreplace]=sigmap1;

	return 1;
}

// change polarization details
void FisherCMB::changePolarization(int xmax1)
{
	xmax=xmax1;

	 if(xmax==1)  cout <<"CMB Errors: Temperature only"<<endl;
 	 if(xmax==2)  cout <<"CMB Errors: T+TE"<<endl;
 	 if(xmax==3)  cout <<"CMB Errors: T+TE+EE"<<endl;
  	 if(xmax==4)  cout <<"CMB Errors: T+TE+EE+BB"<<endl;

}

///////////////////////////////////////////////////////////////////////
// Member functions
//////////////////////////////////////////////////////////////////////
void FisherCMB::setLMax(int lmax1)
{
	lmax=lmax1;
}

void FisherCMB::setLMinPol(int lmin1)
{
	lminPol=lmin1;
}



////////////////////////////////////////////////////////////////////////
// CMB Fisher matrix code
///////////////////////////////////////////////////////////////////////

//Evaluate the CMB fisher matrix
void FisherCMB::fisherMatrixCMB()
{
  double fM;
  string file_dCli;
  string file_dClj;

  string namet, base, tag;
  int i,j,l,m;
  int lmin(2);
  double luse;
  ifstream fini,finj;
  ifstream fin_ref;
  ofstream fout;
  double *dCli, *dClj,*Cl;
  string name=dirbase;
  string file_fid;
  CosmParam *c;

  int nmat;
  double **fisher, **ifisher;

  c=&fiducialFC;

  file_fid=c->cl_file;

  nparamFC=param_tags.size()-1;
  string file_deriv[nparamFC];
//for(i=0;i<=nparamFC;i++) cout<<param_tags[i]<<endl;
  nmat=nparamFC;
  fisher=dmatrix(1,nmat,1,nmat);
  ifisher=dmatrix(1,nmat,1,nmat);

  dCli=dvector(1,4);
  dClj=dvector(1,4);
  Cl=dvector(1,4);

  //specify derivative files
  tag="_cl_deriv.dat";
  for(i=0;i<=nparamFC-1;i++) { file_deriv[i]=name+param_tags[i+1]+tag;
  cout<<file_deriv[i]<<endl;}

  //everything is done using files to create the fisher matrix
  //I'll exploit this to automate the procedure

  //parameter loops
  for(i=0;i<nmat;i++){
    for(j=0;j<=i;j++){
      fM=0.0;
      file_dCli=file_deriv[i];
      file_dClj=file_deriv[j];
      fini.open(file_dCli.c_str());
      finj.open(file_dClj.c_str());
      fin_ref.open(file_fid.c_str());
      //sum over l
      for(l=2;l<=lmax;l++){
	//use Cl in (\mu K)^2
	//reorder to better suit experimental possibilities
	//x=(T,E,B,C) -> x=(T,C,E,B)
	fini>>luse>>dCli[1]>>dCli[3]>>dCli[4]>>dCli[2];
	finj>>luse>>dClj[1]>>dClj[3]>>dClj[4]>>dClj[2];
	fin_ref>>luse>>Cl[1]>>Cl[3]>>Cl[4]>>Cl[2];
	//recover the Cl from l(l+1)Cl/2*PI
	for(m=1;m<=4;m++){
	  Cl[m]/=double(luse*(luse+1))/(2.0*PI);
	  dCli[m]/=double(luse*(luse+1))/(2.0*PI);
	  dClj[m]/=double(luse*(luse+1))/(2.0*PI);
	}

	if(l>=lmin) fM+=formFisherCMB(l,Cl,dCli,dClj);
	//cout <<i<<"\t"<<j<<"\t"<<luse<<"\t"<<fM<<endl;
      }
      fisher[i+1][j+1]=fM;
      if(i!=j) fisher[j+1][i+1]=fM;  //fill in matrix by symmetry
      fini.close();
      finj.close();
      fin_ref.close();
           // cout <<i<<"\t"<<j<<"\t"<<fM<<endl;
    }
  }

//cout << " Okay " << endl;
  free_dvector(dCli,1,4);
  free_dvector(dClj,1,4);
  free_dvector(Cl,1,4);

  //save information about fisher matrix and inverse into class
  setFisherMatrix(fisher,nmat,c);

  //now write results to file
  if(xmax==1) tag="_T.dat";
  if(xmax==2) tag="_TC.dat";
  if(xmax==3) tag="_TCE.dat";
  if(xmax==4) tag="_TCEB.dat";

  saveFisher(tag);

  cout<<exp_name<<endl;
  cout<<nparamFC<<endl;
  if(xmax==1)  cout <<"CMB Errors: Temperature only"<<endl;
  if(xmax==2)  cout <<"CMB Errors: T+TE"<<endl;
  if(xmax==3)  cout <<"CMB Errors: T+TE+EE"<<endl;
  if(xmax==4)  cout <<"CMB Errors: T+TE+EE+BB"<<endl;

  for(i=1;i<=nparamFC;i++){
	cout<<"delta"<<param_tags[i]<<"\t"<<param_values[i]<<"\t"<<sqrt(ifisherFC[i][i])<<"\t"<<sqrt(ifisherFC[i][i])/param_values[i]<<endl;
  }

  free_dmatrix(fisher,1,nmat,1,nmat);
  free_dmatrix(ifisher,1,nmat,1,nmat);
}


double FisherCMB::formFisherCMB(int luse,double *Cl,double *dCli,double *dClj)
{
  double fFC(0.0);
  int x,y;
  double **cov,**icov;

  cov=dmatrix(1,xmax,1,xmax);
  icov=dmatrix(1,xmax,1,xmax);

  //form covariance matrix
  for(x=1;x<=xmax;x++){
    for(y=1;y<=xmax;y++){
      cov[x][y]=covCMB(luse,x,y,Cl);
    }
  }
  //obtain inverse covariance matrix
  invertMatrix(cov,xmax,icov);

  //loop over polarisation types
  for(x=1;x<=xmax;x++){
    for(y=1;y<=xmax;y++){
      fFC+=dCli[x]*dClj[y]*icov[x][y];
    }
  }

  free_dmatrix(cov,1,xmax,1,xmax);
  free_dmatrix(icov,1,xmax,1,xmax);

  return fFC;
}

//covariance matrix for a CMB experiment
//x=(T,C,E,B)
double FisherCMB::covCMB(int l,int x,int y,double *Cl)
{
  double cov;
  double bl2;
  double bl2temp;
  double ll;
  double wtbl2, wpbl2;
  static double iwtbl2,iwpbl2;
  static double prefactor;
  static int lsave;
	int i;

  ll=double(l);

  //calculate weighting functions anew for each value of l
  if(lsave !=l){
    lsave=l;
	//sum contribution from all channels
	wtbl2=0.0;
	wpbl2=0.0;
	for(i=1;i<=nchannel;i++){
    		bl2temp=pow(beam[i],2.0)*ll*(ll+1.0);
    		bl2temp/=-8.0*log(2.0);
    		bl2=exp(bl2temp);
    		wtbl2+=pow(sigmat[i]*beam[i],-2.0)*bl2;
    		wpbl2+=pow(sigmap[i]*beam[i],-2.0)*bl2;
	}
    iwtbl2=1.0/wtbl2;
    iwpbl2=1.0/wpbl2;
    prefactor=2.0/(2.0*ll+1.0)/fsky;
  }

  //impose cutoff for low l polarization
  if(l<lminPol && (x!=1 || y!=1)) prefactor*=1.0e30;

  //evaluate covariance matrices
  //x=(T,C,E,B)  note this is different from CAMB output ordering

  cov=prefactor;
  if(x==1 && y==1) cov*=pow(Cl[1]+iwtbl2,2.0);
  if(x==1 && y==3) cov*=pow(Cl[2],2.0);
  if(x==1 && y==2) cov*=Cl[2]*(Cl[1]+iwtbl2);
  if(x==1 && y==4) return 0.0;
  if(x==3 && y==1) cov*=pow(Cl[2],2.0);
  if(x==3 && y==3) cov*=pow(Cl[3]+iwpbl2,2.0);
  if(x==3 && y==2) cov*=Cl[2]*(Cl[3]+iwpbl2);
  if(x==3 && y==4) return 0.0;
  if(x==2 && y==1) cov*=Cl[2]*(Cl[1]+iwtbl2);
  if(x==2 && y==3) cov*=Cl[2]*(Cl[3]+iwpbl2);
  if(x==2 && y==2) cov*=(pow(Cl[2],2.0)+(Cl[1]+iwtbl2)*(Cl[3]+iwpbl2))/2.0;
  if(x==2 && y==4) return 0.0;
  if(x==4 && y==1) return 0.0;
  if(x==4 && y==3) return 0.0;
  if(x==4 && y==2) return 0.0;
  if(x==4 && y==4) cov*=pow(Cl[4]+iwpbl2,2.0);

  return cov;
}



////////////////////////////////////////////////////////////////////////
// Dark energy parameters from CMB
///////////////////////////////////////////////////////////////////////
//
// There is a strict degeneracy in the CMB between \Omega_k, \Omega_\Lambda
// and dark energy parameters, since these enter only through the angular
// diameter distance to the surface of last scattering.  This is to some extent
// broken by the late time ISW effect at small l.  It is useful to be able
// to enforce the degeneracy exactly.
//
// 1) Calculate fisher for flat \Lambda CDM Universe
// 2) Convert constraint on \Omega_\Lambda to constraint on dark energy
//    parameters via dependence on angular size of the sound horizon
//    \theta_s=r_s/D_A
//

// Project from logH, logD_A onto dark energy parameters
//
// de_flag=1 oml,  2: oml, w0  3: oml, w0, w1
//
void FisherCMB::projectOmegaL2AngularScale()
{
   int nparamNew;
   double **fisherNew;
   int i,j, n, m;
   vector<string> tags_new;
   vector<double> values_new;
   double fnm;

//create new fisher matrix
  tags_new=param_tags;
  values_new=param_values;

//erase references to oml from new parameter list
   for(j=0;j<=param_tags.size()-1;j++){
	for(i=0;i<tags_new.size();i++){
  	   if(tags_new[i].find("_oml",0)!=string::npos){
		 tags_new.erase(tags_new.begin()+i,tags_new.begin()+i+1);
		 values_new.erase(values_new.begin()+i,values_new.begin()+i+1);
	   }
	}
   }

  //add thetas parameter
  tags_new.push_back("_thetas");
  values_new.push_back(angularScale(&fiducialFC));

  nparamNew=tags_new.size()-1;
  fisherNew=dmatrix(1,nparamNew,1,nparamNew);

	for(i=0;i<=tags_new.size()-1;i++) cout<<tags_new[i]<<endl;

	//form new fisher matrix from the old
	for(i=1;i<=nparamNew;i++){
		for(j=1;j<=nparamNew;j++){
			fisherNew[i][j]=0.0;
			for(n=1;n<=nparamFC;n++){
				for(m=1;m<=nparamFC;m++){
					fnm=thetaDerivatives(param_tags[n],tags_new[i]);
					fnm*=thetaDerivatives(param_tags[m],tags_new[j]);
					fnm*=fisherFC[n][m];
					fisherNew[i][j]+=fnm;
				}
			}

		}
	}

	setFisherMatrix(fisherNew, nparamNew, &fiducialFC);
	param_tags.clear();
	param_values.clear();

	for(i=0;i<=nparamNew;i++){
		param_tags.push_back(tags_new[i]);
		param_values.push_back(values_new[i]);
	}

	free_dmatrix(fisherNew,1,nparamNew,1,nparamNew);

}

double FisherCMB::thetaDerivatives(string tagn1, string tagi)
{
   double domldtheta, dthetadoml;
   CosmParam cm, cp;
   double pstep;
   string tagn;

   tagn=tagn1;

   if(tagn=="_oml" && tagi=="_thetas"){
	cm=fiducialFC;
	cp=fiducialFC;
	pstep=0.01*(1.0-cp.oml);
	cm.oml-=pstep;
	cp.oml+=pstep;

	dthetadoml=(angularScale(&cp)-angularScale(&cm))/(2.0*pstep);
	domldtheta=1.0/dthetadoml;
	return domldtheta;
   }

  if(tagn==tagi) return 1.0;
  if(tagn!=tagi) return 0.0;

	cout<<"Error: in thetaDerivative for tags: "<<tagn<<"\t"<<tagi<<"\t"<<endl;
	return -1.0;
}

///////////////////////////////////////////////////////////////////////////
// Now project from theta onto dark energy parameters
//////////////////////////////////////////////////////////////////////////

void FisherCMB::projectAngularScaleOntoDE(int de_flag)
{
   int nparamNew;
   double **fisherNew;
   int i,j, n, m;
   vector<string> tags_new;
   vector<double> values_new;
   double fnm;

//create new fisher matrix
  tags_new=param_tags;
  values_new=param_values;

//erase references to oml from new parameter list
   for(j=0;j<=param_tags.size()-1;j++){
	for(i=0;i<tags_new.size();i++){
  	   if(tags_new[i].find("_thetas",0)!=string::npos){
		 tags_new.erase(tags_new.begin()+i,tags_new.begin()+i+1);
		 values_new.erase(values_new.begin()+i,values_new.begin()+i+1);
	   }
	}
   }

  //add dark energy parameters
  tags_new.push_back("_oml");
  values_new.push_back(fiducialFC.oml);
	if(de_flag>1){
		 tags_new.push_back("_w0");
		values_new.push_back(fiducialFC.w0);
	}
	if(de_flag>2 && de_flag<5){
		 tags_new.push_back("_w1");
		values_new.push_back(fiducialFC.w1);
	}
	if(de_flag>3){
		 tags_new.push_back("_omk");
		values_new.push_back(fiducialFC.omk);
	}


  nparamNew=tags_new.size()-1;
  fisherNew=dmatrix(1,nparamNew,1,nparamNew);

	for(i=0;i<=tags_new.size()-1;i++) cout<<tags_new[i]<<endl;

	//form new fisher matrix from the old
	for(i=1;i<=nparamNew;i++){
		for(j=1;j<=nparamNew;j++){
			fisherNew[i][j]=0.0;
			for(n=1;n<=nparamFC;n++){
				for(m=1;m<=nparamFC;m++){
					fnm=cmbDerivatives(param_tags[n],tags_new[i]);
					fnm*=cmbDerivatives(param_tags[m],tags_new[j]);
					fnm*=fisherFC[n][m];
					fisherNew[i][j]+=fnm;
				}
			}

		}
	}

	setFisherMatrix(fisherNew, nparamNew, &fiducialFC);
	param_tags.clear();
	param_values.clear();

	for(i=0;i<=nparamNew;i++){
		param_tags.push_back(tags_new[i]);
		param_values.push_back(values_new[i]);
	}

	free_dmatrix(fisherNew,1,nparamNew,1,nparamNew);

}

double FisherCMB::cmbDerivatives(string tagn1, string tagi)
{
   double dthetadoml;
   CosmParam cm, cp;
   double pstep;
   string tagn;

   tagn=tagn1;

   if(tagn=="_thetas" && tagi=="_oml"){
	cm=fiducialFC;
	cp=fiducialFC;
	pstep=0.01*(1.0-cp.oml);
	cm.oml-=pstep;
	cp.oml+=pstep;

	dthetadoml=(angularScale(&cp)-angularScale(&cm))/(2.0*pstep);
	return dthetadoml;
   }else if(tagn=="_thetas" && tagi=="_w0"){
	cm=fiducialFC;
	cp=fiducialFC;
	pstep=0.01*cp.w0;
	cm.w0-=pstep;
	cp.w0+=pstep;

	dthetadoml=(angularScale(&cp)-angularScale(&cm))/(2.0*pstep);
	return dthetadoml;
   }else if(tagn=="_thetas" && tagi=="_w1"){
	cm=fiducialFC;
	cp=fiducialFC;
	pstep=0.001;
	cm.w1-=pstep;
	cp.w1+=pstep;

	dthetadoml=(angularScale(&cp)-angularScale(&cm))/(2.0*pstep);
	return dthetadoml;
   }else if(tagn=="_thetas" && tagi=="_omk"){
	cm=fiducialFC;
	cp=fiducialFC;
	pstep=0.001;
	cm.omk-=pstep;
	cp.omk+=pstep;

	dthetadoml=(angularScale(&cp)-angularScale(&cm))/(2.0*pstep);
	return dthetadoml;
   }else if(tagn=="_thetas" && tagi=="_ommhh"){
	cm=fiducialFC;
	cp=fiducialFC;
	pstep=0.01*cm.ommhh;
	cm.ommhh-=pstep;
	cp.ommhh+=pstep;

	dthetadoml=(angularScale(&cp)-angularScale(&cm))/(2.0*pstep);
	return dthetadoml;
   }else if(tagn=="_thetas" && tagi=="_ombhh"){
	cm=fiducialFC;
	cp=fiducialFC;
	pstep=0.01*cm.ombhh;
	cm.ombhh-=pstep;
	cp.ombhh+=pstep;

	dthetadoml=(angularScale(&cp)-angularScale(&cm))/(2.0*pstep);
	return dthetadoml;
   }else if(tagn=="_thetas" && tagi=="_Mnu"){
      //assume neutrinos go non-relativistic after CMB decoupling
      //then varying Mnu is like changing ommhh early on
      //since this is dominated by the cold dark matter fraction then
	cm=fiducialFC;
	cp=fiducialFC;
	pstep=0.01*cm.omnuhh;
	cm.omnuhh-=pstep;
	cp.omnuhh+=pstep;

	dthetadoml=(angularScale(&cp)-angularScale(&cm))/(2.0*pstep);
	return dthetadoml*omnuhhFromMnu(1.0);
   }



  if(tagn==tagi) return 1.0;
  if(tagn!=tagi) return 0.0;

	cout<<"Error: in deDerivative for tags: "<<tagn<<"\t"<<tagi<<"\t"<<endl;
	return -1.0;
}

///////////////////////////////////////////////////////////////////////////
// Wrapper to take automatically calculate fisher matrix using above technique
//////////////////////////////////////////////////////////////////////////

//  de_flag=-1  thetas
//  de_flag=0
//
void FisherCMB::fisherMatrixCMB_DE(int de_flag)
{
	int i, oml_flag(0);
 	double thetas;

	oml_flag=1;

	//need to look at fiducial parameter set and see if dark energy parameters
	//are being used

	for(i=1;i<=param_tags.size()-1;i++){
		if(param_tags[i]=="_oml"){
			if(1>de_flag) {
				de_flag=1;
				oml_flag=1;
			}
		}
		if(param_tags[i]=="_w0"){if(2>de_flag) de_flag=2;}
		if(param_tags[i]=="_w1"){if(3>de_flag) de_flag=3;}
		if(param_tags[i]=="_omk"){if(4>de_flag) de_flag=4;}
	}

	//Conversion to thetas based on having oml about
	if(oml_flag==0){
		cout<<"Error must have oml as parameter"<<endl;
		return;
	}

	//empty parameter set of non-Lambda CDM parameters
	removeParam("_oml");
	if(de_flag>1) removeParam("_w0");
	if(de_flag>2) removeParam("_w1");
	if(de_flag>3) removeParam("_omk");

	//add in angular scale as parameter
	thetas=angularScale(&fiducialFC);

	addParamWithValue("_thetas",thetas);

	//calculate fisher matrix
	fisherMatrixCMB();

	//project onto dark energy parameters
        if(de_flag>0){
           projectAngularScaleOntoDE(de_flag);
        }

	printErrors();
}



////////////////////////////////////////////////////////////////////
//  Set Fisher matrix data files
///////////////////////////////////////////////////////////////////
//
// Initialise files for use with theta_s as a replacement parameter for any involving
// dark energy and curvature.

double FisherCMB::setCosmFilesCMB()
{
  int i;
  double z, thetas;

   //remove dark energy and curvature parameters
	for(i=1;i<=param_tags.size()-1;i++){
		if(param_tags[i]=="_oml") removeParam("_oml");
		if(param_tags[i]=="_w0") removeParam("_w0");
		if(param_tags[i]=="_w1") removeParam("_w1");
		if(param_tags[i]=="_omk") removeParam("_omk");
	}


  //add in angular scale as parameter
  thetas=angularScale(&fiducialFC);
  addParamWithValue("_thetas",thetas);

  //set redshifts to calculate matter power spectrum

	redshift_list.clear();
	redshift_list.push_back(0.0);

  //parameters for calculating derivatives
  for(i=0;i<=nparamFC;i++) setParamPairCMB(&fiducialFC,i,1);
  return 1.0;
}

double FisherCMB::setParamPairCMB(CosmParam *fiducial, int iflag, int clflag)
{
  double p_step;
  CosmParam c_p;
  CosmParam c_m;
  Neutrinos NuClass;
  double Mnu;
  double vary(0.05);  //control derivative step sizes- optimal value??
  double vary2(0.1);
  string name, file;
  int lmax(2500);
  vector<double> sigma8list_p,  sigma8list_m;
  double thetas;

  //initialise with fiducial values
  c_p=loadCosmParam(fiducial->name);
  c_m=loadCosmParam(fiducial->name);

  name=dirbase;
  file="default";
  p_step=0.0;


	cout<<param_tags[iflag]<<endl;

  //use if statements to specify which parameter to modify
  if(param_tags[iflag]=="_fiducial"){
    file=name+param_tags[iflag];
  }else if(param_tags[iflag]=="_ommhh"){
    cout<<"creating ommhh"<<endl;
    file=name+param_tags[iflag];
  p_step=vary*fiducial->ommhh/4.0;

    thetas=angularScale(&c_p);

  c_p.ommhh+=p_step;
  c_m.ommhh-=p_step;

  c_p.oml=fixThetas(&c_p,thetas);
  c_m.oml=fixThetas(&c_m,thetas);

  }else if(param_tags[iflag]=="_ombhh"){
    file=name+param_tags[iflag];
  p_step=vary*fiducial->ombhh;

  thetas=angularScale(&c_p);
  c_p.ombhh+=p_step;
  c_m.ombhh-=p_step;

  c_p.oml=fixThetas(&c_p,thetas);
  c_m.oml=fixThetas(&c_m,thetas);

  }else if(param_tags[iflag]=="_thetas"){
    file=name+param_tags[iflag];

  thetas=angularScale(&c_p);
  p_step=vary*thetas;

  c_p.oml=fixThetas(&c_p,thetas+p_step);
  c_m.oml=fixThetas(&c_m,thetas-p_step);

  }else if(param_tags[iflag]=="_ascal2"){
    file=name+param_tags[iflag];
  p_step=vary*fiducial->ascal2;
  c_p.ascal2+=p_step;
  c_m.ascal2-=p_step;

  }else if(param_tags[iflag]=="_nscal"){
    file=name+param_tags[iflag];
  p_step=vary*fiducial->nscal;
  p_step=0.005/2.0;               // note fixed step size
  c_p.nscal+=p_step;
  c_m.nscal-=p_step;

  }else if(param_tags[iflag]=="_tau"){
    file=name+param_tags[iflag];
  p_step=vary*fiducial->tau;
  c_p.tau+=p_step;
  c_m.tau-=p_step;

  }else if(param_tags[iflag]=="_ts"){
    file=name+param_tags[iflag];
  p_step=vary*fiducial->ts/4.0;
  if(p_step<1.0e-6) p_step=0.001;
  c_p.ts+=p_step;
  c_m.ts-=p_step;

  if(c_m.ts<0.0){
	  p_step=0.001;               // note fixed step size
  	c_p.ts=fiducial->ts+2.0*p_step;        //one-sided derivative
  	c_m.ts=0.0;
  }

  }else if(param_tags[iflag]=="_w0"){
    file=name+param_tags[iflag];
  p_step=vary2*fiducial->w0;
  c_p.w0+=p_step;
  c_m.w0-=p_step;

  }else if(param_tags[iflag]=="_alpha"){
    file=name+param_tags[iflag];
  //p_step=vary*fiducial->alpha;
  p_step=0.0005;                //note fixed step size
  c_p.alpha+=p_step;
  c_m.alpha-=p_step;

  }else if(param_tags[iflag]=="_h"){
    file=name+param_tags[iflag];
  p_step=vary*fiducial->h;
  c_p.h+=p_step;
  c_m.h-=p_step;
  //MIGHT NEED TO MODIFY OMMHH AND OMBHH WHEN USING THIS

  }else if(param_tags[iflag]=="_omk"){
  //being sloppy with this as should treat curvature more carefully
    file=name+param_tags[iflag];
  p_step=vary2*fiducial->omm;
  c_p.omk+=p_step;
  c_m.omk-=p_step;
  c_p.omm=1.0-c_p.oml-c_p.omk;   //ensure consistency
  c_m.omm=1.0-c_m.oml-c_m.omk;   //ensure consistency
  c_p.h=sqrt(c_p.ommhh/c_p.omm);
  c_m.h=sqrt(c_m.ommhh/c_m.omm);

  }else if(param_tags[iflag]=="_omnuhh"){
    file=name+param_tags[iflag];
  p_step=vary2*fiducial->omnuhh;
  thetas=angularScale(&c_p);

  c_p.omnuhh+=p_step;
  c_m.omnuhh-=p_step;

  c_p.oml=fixThetas(&c_p,thetas);
  c_m.oml=fixThetas(&c_m,thetas);

  }else if(param_tags[iflag]=="_yHe"){
    file=name+param_tags[iflag];
  p_step=vary2*fiducial->yHe;
  c_p.yHe+=p_step;
  c_m.yHe-=p_step;
  }else if(param_tags[iflag]=="_omnu"){
    file=name+param_tags[iflag];
  p_step=vary2*fiducial->omnu;
    if(p_step<1.0e-6) p_step=0.005;
  c_p.omnu+=p_step;
  c_m.omnu-=p_step;
  if(c_m.omnu<0.0){
	c_p.omnu=fiducial->omnu+2.0*p_step;    //one sided derivative
	c_m.omnu=fiducial->omnu;
  }
  c_p.omnuhh=c_p.omnu*c_p.h*c_p.h;
  c_m.omnuhh=c_m.omnu*c_m.h*c_m.h;

  //PROBABLY NEED TO MODIFY CALL_CAMB FOR OMNU TO WORK WITH OMMHH AND OML
  }else if(param_tags[iflag]=="_Mnu"){
    file=name+param_tags[iflag];
   p_step=vary2*fiducial->omnuhh;

   if(p_step<1.0e-6) p_step=0.00005;

   Mnu=NuClass.mnuFromOmnuhh(fiducial->omnuhh);
   NuClass.initNeutrinos(Mnu,3,3,hierarchy_flag);  // normal hierarchy
   Mnu=NuClass.minimumNeutrinoMass();

  c_p.omnuhh+=p_step;
  c_m.omnuhh-=p_step;
  if(mnuFromOmnuhh(c_m.omnuhh)<Mnu || c_m.omnuhh<0.0){  //one sided derivative
	c_m.omnuhh=fiducial->omnuhh;
	c_p.omnuhh=fiducial->omnuhh+2.0*p_step;
  }
   p_step=(mnuFromOmnuhh(c_p.omnuhh)-mnuFromOmnuhh(c_m.omnuhh))/2.0;

	//Mnu, m2, m3 parametreization
	if(neutrino_mass_flag_alt==1){
 		Mnu=NuClass.mnuFromOmnuhh(c_p.omnuhh);
		c_p.numass1=Mnu-c_p.numass2-c_p.numass3;
 		Mnu=NuClass.mnuFromOmnuhh(c_m.omnuhh);
		c_m.numass1=Mnu-c_m.numass2-c_m.numass3;
	}

  }else if(param_tags[iflag]=="_numass1"){
        file=name+param_tags[iflag];
	p_step=vary*c_p.numass1*2.0;
	p_step=0.008;
	c_p.numass1+=p_step;
	c_p.omnuhh=omnuhhFromMnu(c_p.numass1+c_p.numass2+c_p.numass3);
	c_m.numass1-=p_step;
	c_m.omnuhh=omnuhhFromMnu(c_m.numass1+c_m.numass2+c_m.numass3);
  }else if(param_tags[iflag]=="_numass2"){
        file=name+param_tags[iflag];
	p_step=vary*c_p.numass2;
	p_step=0.008;
	c_p.numass2+=p_step;
	c_p.omnuhh=omnuhhFromMnu(c_p.numass1+c_p.numass2+c_p.numass3);
	c_m.numass2-=p_step;
	c_m.omnuhh=omnuhhFromMnu(c_m.numass1+c_m.numass2+c_m.numass3);

	//Mnu, m2, m3 parametreization
	if(neutrino_mass_flag_alt==1){
		c_p=fiducialFC;
		c_m=fiducialFC;
		c_p.numass2+=p_step;
		c_p.numass1-=p_step;
		c_m.numass2-=p_step;
		c_m.numass1+=p_step;
	}

  }else if(param_tags[iflag]=="_numass3"){
        file=name+param_tags[iflag];
	p_step=vary*c_p.numass3;
	p_step=0.008;
	c_p.numass3+=p_step;
	c_p.omnuhh=omnuhhFromMnu(c_p.numass1+c_p.numass2+c_p.numass3);
	c_m.numass3-=p_step;
	c_m.omnuhh=omnuhhFromMnu(c_m.numass1+c_m.numass2+c_m.numass3);

	//Mnu, m2, m3 parametreization
	if(neutrino_mass_flag_alt==1){
		c_p=fiducialFC;
		c_m=fiducialFC;
		c_p.numass3+=p_step;
		c_p.numass1-=p_step;
		c_m.numass3-=p_step;
		c_m.numass1+=p_step;
	}

  }else if(param_tags[iflag]=="_epsilon"){
        file=name+param_tags[iflag];
	p_step=vary*c_p.ts/8.0;
	p_step=0.0005;
	if(fabs(p_step)<1.0e-6) p_step=0.0005;
	c_p.ts+=p_step;
	c_m.ts-=p_step;
  }else if(param_tags[iflag]=="_eta"){
        file=name+param_tags[iflag];
	p_step=vary*c_p.nscal/2.0;
	if(fabs(p_step)<1.0e-6) p_step=0.0005;
	c_p.nscal+=p_step;
	c_m.nscal-=p_step;
  }else if(param_tags[iflag]=="_xi"){
        file=name+param_tags[iflag];
	p_step=vary*c_p.alpha;
	if(fabs(p_step)<1.0e-6) p_step=0.0005;
	c_p.alpha+=p_step;
	c_m.alpha-=p_step;
  }else if(param_tags[iflag]=="_w1"){
        file=name+param_tags[iflag];
	p_step=vary*c_p.w1;
	if(p_step<1.0e-6) p_step=0.001;   //NEED TO CHECK THIS VALUE WORKS
	c_p.w1+=p_step;
	c_m.w1-=p_step;
  }else if(param_tags[iflag]=="_ntens"){
        file=name+param_tags[iflag];
	p_step=vary*c_p.ntens;
	if(p_step<1.0e-6) p_step=0.001;   //NEED TO CHECK THIS VALUE WORKS
	c_p.ntens+=p_step;
	c_m.ntens-=p_step;
  }

  //now initialise the two cosmology files
	c_p.file=file+"_p.dat";
	c_m.file=file+"_m.dat";
	c_p.cl_file=file+"_cl_p.dat";
	c_m.cl_file=file+"_cl_m.dat";
	c_p.filebase=param_tags[iflag];
	c_m.filebase=param_tags[iflag];

  c_p.p_step=p_step;
  c_m.p_step=p_step;
  sigma8list_p=callCAMB(&c_p);
  sigma8list_m=callCAMB(&c_m);

  c_p.sigma8=sigma8list_p[sigma8list_p.size()-1];
  c_m.sigma8=sigma8list_m[sigma8list_m.size()-1];

  c_p.name=file+"_param_p.dat";
  c_m.name=file+"_param_m.dat";
  saveCosmParam(&c_m,c_m.name);
  saveCosmParam(&c_p,c_p.name);

  // calculate the derivatives of the Cl
  if(clflag){
	  file=file+"_cl_deriv.dat";
	  derivCMB(lmax,p_step,&c_m,&c_p,file);
  }

  //since the parameters are redshift independent copy them
  //into files marked with the redshift - need for biases
  int i;
  for(i=0;i<=redshift_list.size()-1;i++){
	file=dirbase+"_z"+numToString(redshift_list[i])+param_tags[iflag];
	c_p.exp_z=redshift_list[i];
	c_m.exp_z=redshift_list[i];
	c_p.sigma8=sigma8list_p[i];
	c_m.sigma8=sigma8list_m[i];
	c_p.name=file+"_param_p.dat";
  	c_m.name=file+"_param_m.dat";
  	saveCosmParam(&c_m,c_m.name);
  	saveCosmParam(&c_p,c_p.name);
  }

  return 1.0;
}
