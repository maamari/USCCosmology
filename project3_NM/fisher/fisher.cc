/* fisher.cc
 * Contains functions for calculating fisher matrix
  */

#include <math.h>
#include <iostream>
#include <fstream>
#include "fisher.h"
#include "dcosmology2.h"
#include "dnumrecipes.h"
#include "spline.h"
#include <sstream>
#include "neutrinos.h"
#include "inflation.h"

using namespace std;

//References:
//Spitzer:  Spitzer, "Physical Processes in the ISM".


/*********************************************************************
 **************** Constructor/Destructor *****************************
 ********************************************************************/

 Fisher::Fisher()
{
	//cout<<"Fisher constructor called"<<endl;
   fisher_flag=0;    // flag to indicate whether fisher matrix exists
	hierarchy_flag=1; //default to normal hierarchy
	neutrino_mass_flag=0; //use individual neutrino masses
	neutrino_mass_flag_alt=0; //use m1,m2,m3 parametrization
	distance_flag=0;      //treat distances as independent parameters
	slowroll_flag=0;      //use slow roll inflationary parametrization
	lensing_flag=1;         //use weak lensing of c_l
	omnu_flag=0;
        omk_flag=1;           //use d_A for calculating omk derivatives
                              //NOT YET FIXED FOR P(K) ONLY FOR C_L
        reionization_flag=0;  //use tanh function for ionization history
 	transfer_spline_flag=1;  //calculate transfer functions ???
	//dirbase="./checkdir/unlensed";
	dirbase="./data/cosmology";
}

 Fisher::~Fisher()
{
	//cout<<"Fisher destructor called"<<endl;

}


///////////////////////////////////////////////////////////////////////
// Member functions
//////////////////////////////////////////////////////////////////////
//set neutrino hierarchy flag;
int Fisher::setHierarchyFlag(int hflag)
{
	hierarchy_flag=hflag;
	return hierarchy_flag;
}

//set experiment redshift
void Fisher::setExpZ(double z1)
{
	exp_z=z1;
}

//set slow roll flag
int Fisher::setSlowrollFlag(int srflag)
{
	slowroll_flag=srflag;
	//if(slowroll_flag==1) dirbase="./data_sr/cosmology";
	if(slowroll_flag==1) transfer_spline_flag=0;
	return slowroll_flag;
}

//set lensing flag
int Fisher::setLensingFlag(int lflag)
{
	lensing_flag=lflag;
	return lensing_flag;
}

//set lensing flag
int Fisher::setAltNeutrinoParam()
{
	neutrino_mass_flag_alt=1;
	removeParam("_Mnu");
	removeParam("_numass1");
	removeParam("_numass2");
	removeParam("_numass3");
	addParam("_Mnu");
	addParam("_numass2");
	addParam("_numass3");

	return neutrino_mass_flag_alt;
}

void Fisher::assignSpline(Spline *spline, string name, string tag)
{
	string file;

	file=name;
	if(transfer_spline_flag==1){
	//transfer function spline
          //cout << name+"_transfer"+tag << endl;
		spline->loadFileSpline(name+"_transfer"+tag,7,1,2);
	}else{
	//matter power spline
          //cout << name+tag << endl;
		spline->loadFileSpline(name+tag,2,1,2);
	}
	return;
}

void Fisher::assignSpline2(Spline *spline, string name, string tag, int flag)
{
	string file;

	file=name;
	if(flag==1){
	//transfer function spline
          //cout << name+"_transfer"+tag << endl;
		spline->loadFileSpline(name+"_transfer"+tag,7,1,2);
	}else{
	//matter power spline
          //cout << name+tag << endl;
		spline->loadFileSpline(name+tag,2,1,2);
	}
	return;
}

// set base string for directory to save files to
//
void Fisher::setDirbase(string name)
{
	string temp;

	temp=name+"/cosmology";
	dirbase=temp;
	return;
}


//set internal fisher matrix
// - automatically calculates inverse fisher matrix
//
void Fisher::setFisherMatrix(double **fisher, int nparam, CosmParam *fiducial)
{
	int i,j;

	if(fisher_flag){
		free_dmatrix(fisherFC,1,nparamFC,1,nparamFC);
		free_dmatrix(ifisherFC,1,nparamFC,1,nparamFC);
		fisher_flag=0;
	}

	fisher_flag=1;
	nparamFC=nparam;
	fisherFC=dmatrix(1,nparamFC,1,nparamFC);
	ifisherFC=dmatrix(1,nparamFC,1,nparamFC);
	for(i=1;i<=nparamFC;i++){
		for(j=1;j<=nparamFC;j++){
			fisherFC[i][j]=fisher[i][j];
		}
	}

	// invertMatrix(fisherFC,nparamFC,ifisherFC);

	fiducialFC=loadCosmParam(fiducial->name);

}

// Function to set up CosmParam in sensible way
CosmParam Fisher::setCosmParam(double h, double omm, double omb, double oml, double omnu, double tau, double ascal2, double nscal, double alpha, double ts, double w0, double w1, double sigma8, string name)
{
  CosmParam c;
  Neutrinos NuClass;
  double Mnu;
  vector<double> mass_splittings;
  string base,file,tag;
  vector<double> sigma8list;

  //set parameters
  c.physical_flag=1;
  c.h=h;
  c.omm=omm;
  c.omb=omb;
  c.oml=oml;
  c.omnu=omnu;
  c.tau=tau;
  c.ascal2=ascal2;
  c.nscal=nscal;
  c.alpha=alpha;
  c.ts=ts;
  c.w0=w0;
  c.w1=w1;

 //default redshift
  c.exp_z=0.3;

 //derived quantities
  c.omk=1.0-omm-oml;
  c.yHe=0.24;          //default value of yHe=0.24
  c.ommhh=c.omm*c.h*c.h;
  c.ombhh=c.omb*c.h*c.h;
  c.omnuhh=c.omnu*c.h*c.h;
  c.p_step=0.0;
  c.bias=1.0;

   //set up neutrino information
   Mnu=NuClass.mnuFromOmnuhh(c.omnuhh);
   NuClass.initNeutrinos(Mnu,3,3,hierarchy_flag);  // normal hierarchy
   mass_splittings=NuClass.massSplittings();
   c.numass1=Mnu*mass_splittings[0];
   c.numass2=Mnu*mass_splittings[1];
   c.numass3=Mnu*mass_splittings[2];

   // Cluster parameters
   c.sigma0= 0.25;
   c.Bm0 = -0.15;
   c.alpha= 0;
   c.beta = 0;
   c.fR0= 0;

  base=dirbase;
  tag=".dat";
  file=base+name+tag;

  c.file=file;
  c.filebase=name;

  tag="_cl.dat";
  file=base+name+tag;
  c.cl_file=file;

  if(sigma8<0.0){
	  sigma8list=callCAMB(&c);
	c.sigma8=sigma8list[sigma8list.size()-1];
  }else{
	c.sigma8=sigma8;
  }

  base=dirbase;
  tag="_param.dat";
  file=base+name+tag;
  c.name=file;
  saveCosmParam(&c,c.name);

  return c;
}

//return fisher matrix
int Fisher::getFisherMatrix(double **fisher)
{
	int i,j;

	for(i=1;i<=nparamFC;i++){
		for(j=1;j<=nparamFC;j++){
			fisher[i][j]=fisherFC[i][j];
		}
	}

	return nparamFC;
}

//return fisher matrix
int Fisher::getInverseFisherMatrix(double **ifisher)
{
	int i,j;

	for(i=1;i<=nparamFC;i++){
		for(j=1;j<=nparamFC;j++){
			ifisher[i][j]=ifisherFC[i][j];
		}
	}

	return nparamFC;
}

int Fisher::getFisherSize()
{
	return nparamFC;
}

void Fisher::setName(string name)
{
	exp_name=name;
}


vector<string> Fisher::returnParam_Tags()
{
	return param_tags;
}

vector<double> Fisher::returnParam_Values()
{
	return param_values;
}

fisherData Fisher::returnFisherData()
{
	fisherData fisherDataFC;

	saveFisher("_out");

	fisherDataFC.nparam=nparamFC_file;
	fisherDataFC.param_tags=param_tags;
	fisherDataFC.param_values=param_values;

	fisherDataFC.fisher_file=fisher_fileFC;
	fisherDataFC.exp_z=exp_z;

	return fisherDataFC;
}

string Fisher::returnDirbase()
{
	return dirbase;
}

// Convert double to string with xxxx.xx fixed precision
//
string Fisher::numToString(double zin)
{
	stringstream ss;
	string str;

	ss.precision(3);
	ss<<zin;
	ss>>str;

	return str;
}

// Convert double to string with xxxx.xx fixed precision
//
string Fisher::numToStringLong(double zin, int np)
{
	stringstream ss;
	string str;

	ss.precision(np);
	ss<<zin;
	ss>>str;

	return str;
}

//copy matrix A into matrix B
void Fisher::copyMatrix(double **A, int nmat, double **B)
{
	int i,j;

	for(i=1;i<=nmat;i++){
		for(j=1;j<=nmat;j++){
			B[i][j]=A[i][j];
		}
	}

}



//calculate redshift of CMB from Hu & White fitting function
double Fisher::getZCMB(CosmParam *fiducial)
{
   double b1,b2;
   double zcmb;
   double temp;

   b1=0.0783*pow(fiducial->ombhh,-0.238);
   //cout<<b1<<"\t"<<fiducial->ombhh<<endl;
   b1/=1.0+39.5*pow(fiducial->ombhh,0.763);

   b2=0.56/(1.0+21.1*pow(fiducial->ombhh,1.81));

   zcmb=1048.0*(1.0+0.00124*pow(fiducial->ombhh,-0.738));
   zcmb*=1.0+b1*pow(fiducial->ommhh,b2);

   //cout<<b1<<"\t"<<b2<<"\t"<<zcmb<<endl;
   return zcmb;
}


///////////////////////////////////////////////////////////////////////
// Distance calculations
//////////////////////////////////////////////////////////////////////
//numerical calculation of  coordDistance

/* Returns angular diameter distance for given cosmology, to redshift
 * z, in units of c/H0.  */
double Fisher::angDiamDistance(double zCall, CosmParam *c)
{
   double omegak,cDist,h;
  omegak=c->omk;
  h=sqrt(c->ommhh/(1.0-c->omk-c->oml));

  cDist=coordDistanceNum(zCall,c);

  if(omegak>1.0e-5){
	 //open
	 return sinh(sqrt(omegak)*h*cDist)/(1.0+zCall)/sqrt(omegak)/h;
  }else if(omegak<-1.0e-5){
	 //closed
	 return sin(sqrt(-omegak)*h*cDist)/(1.0+zCall)/sqrt(-omegak)/h;
  }else{
	//flat
	 return cDist/(1.0+zCall);
  }

  cout<<"Error in angDiamDistance"<<endl;
  return coordDistanceNum(zCall,c)/(1.0+zCall);
}

/* Returns comoving angular diameter distance for given cosmology, to redshift
 * z, in units of c/H0.  */
double Fisher::coAngDiamDistance(double zCall, CosmParam *c)
{
   return angDiamDistance(zCall,c)*(1.0+zCall);
}

/* Returns angular diameter distance for given cosmology, to redshift
 * z, in units of c/H0.  */
double Fisher::lumDistance(double zCall, CosmParam *c)
{
   return angDiamDistance(zCall,c)*pow(1.0+zCall,2.0);
}


double Fisher::coordDistanceNum(double zCall, CosmParam *c)
{
	double cDist;
	double tol(1.0e-5);

	setCoordDistanceKernelF(0.0,c,1);
  	cDist=qromb(getCoordDistanceKernelF,0.0,zCall,tol);

	return cDist;
}

double setCoordDistanceKernelF(double z, CosmParam *c1, int iflag)
{
	static CosmParam *c;
	double kernel;
	double hubble, temp;
	double omlZ;
	double omegarhh(4.15e-5), omegar;
	double h;
        double zNU(1100.0);

	if(iflag==1){
		c=c1;
		return 0.0;
	}

	h=sqrt(c->ommhh/(1.0-c->omk-c->oml));
	omegar=omegarhh/h/h;

	temp=(1.0-c->omk-c->oml-omegar)*pow(1.0+z,3.0);
        //if massive neutrinos need to only have nr matter
        if(z>zNU) temp=(c->ommhh-c->omnuhh-omegarhh)/h/h*pow(1.0+z,3.0);
	temp+=omegar*pow(1.0+z,4.0);
	temp+=(c->omk)*pow(1.0+z,2.0);

	omlZ=c->oml*pow(1.0+z,3.0*(1.0+c->w0+c->w1));
	omlZ*=exp(-3.0*c->w1*z/(1.0+z));
	temp+=omlZ;

	hubble=sqrt(temp)*h;

	kernel=1.0/hubble;

	return kernel;
}

double getCoordDistanceKernelF(double z)
{
	return setCoordDistanceKernelF(z,NULL,0);
}

double Fisher::hubbleZ(double z, CosmParam *c)
{
	double temp, hubble;
	double omlZ;
	double omegarhh(4.15e-5), omegar;
	double h;
        double zNU(1100.0);

	h=sqrt(c->ommhh/(1.0-c->omk-c->oml));
	omegar=omegarhh/h/h;

	temp=(1.0-c->omk-c->oml-omegar)*pow(1.0+z,3.0);
        //if massive neutrinos need to only have nr matter
        if(z>zNU) temp=(c->ommhh-c->omnuhh-omegarhh)/h/h*pow(1.0+z,3.0);
	temp+=c->omk*pow(1.0+z,2.0);
	temp+=omegar*pow(1.0+z,4.0);

	omlZ=c->oml*pow(1.0+z,3.0*(1.0+c->w0+c->w1));
	omlZ*=exp(-3.0*c->w1*z/(1.0+z));
	temp+=omlZ;

	hubble=sqrt(temp);
	return hubble*h;
}

double Fisher::omegaMZ(double z, CosmParam *c)
{
	double temp,ommz;
	temp=c->omm*pow(1.0+z,3.0);
	temp+=c->omk*pow(1.0+z,2.0);
	temp+=c->oml*pow(1.0+z,3.0*(1.0+c->w0));

	ommz=c->omm*pow(1.0+z,3.0)/temp;

	return ommz;
}

double Fisher::growthZ(double z, CosmParam *c)
{
	double Dz;
        //Cosmology cosm(c->omm,c->oml,c->omb,c->h,c->sigma8,c->nscal,c->omnu);
        //Spline sp;
        Cosmology cosm(c->omm,c->oml,c->omb,c->h,c->sigma8,c->nscal,c->omnu, c->w0, c->w1);
	Dz=cosm.growthFac(z);

	return Dz;
}

//
// Units : dVdO
//
double Fisher::getDVolumeDOmega(double z, CosmParam *c)
{
	double dVdO;
        //Spline sp;
        //Cosmology cosm(c->omm,c->oml,c->omb,c->h,c->sigma8,c->nscal,c->omnu);
        Cosmology cosm(c->omm,c->oml,c->omb,c->h,c->sigma8,c->nscal,c->omnu,c->w0, c->w1);
	dVdO=cosm.getVolume(z);  //Mpc/sterradian

	dVdO/=pow(180.0/PI,2.0);  //convert to sq. deg.

	return dVdO;
}

//calculate survey volume from observed area on sky
//
// Units:: area  sq. deg.
//	   volume  Mpc^-3
double Fisher::volumeFromArea(double zmin, double zmax, double obsarea)
{
	double volume;
	double V1, V2;

	V2=getDVolumeDOmega(zmax,&fiducialFC);
	V1=getDVolumeDOmega(zmin,&fiducialFC);
  	volume=(V2-V1)*obsarea;

	return volume;
}


//Growth factor defined in EH97 (8,9)
//
// INVALID CALCULATION NEED TO UPDATE TO INTEGRATE GROWTH NUMERICALLY
//
double Fisher::getGrowth(double ause,CosmParam *c)
{
  double D1;
  double tol(1.0e-5);

  setGrowth(1.0,c,1);
  D1=qromb(intGrowth,1.0e-3,ause,tol);
  D1*=5.0*c->omm*setGrowth(ause,c,0)/2.0;

  // cout <<"D1: "<<D1<<endl;
  return D1;
}

double setGrowth(double a,CosmParam *c1,int iflag)
{
  static CosmParam *c;
  double growth;
  double omlZ,z;
  double zmax(21);
  double amin(1.0/zmax);

  if(iflag==1){
    c=c1;
    return 0.0;
  }

  z=1.0/a;
  //limit growth of Oml to redshift z=20
  if(a<amin){
    omlZ=c->oml*pow(zmax,3.0*(1.0+c->w0-c->w1))*exp(3.0*c->w1*(zmax-1.0));
  }else{
    omlZ=c->oml*pow(z,3.0*(1.0+c->w0-c->w1))*exp(3.0*c->w1*(z-1.0));
  }
  growth=sqrt(c->omm*z*z*z+(1.0-c->omm-c->oml)*z*z+omlZ);
  return growth;
}

double intGrowth(double a)
{
  return pow(a*setGrowth(a,NULL,0),-3.0);
}

///////////////////////////////////////////////////////////////////////
// Reionization code for using CAMB
///////////////////////////////////////////////////////////////////////

//tanh model for ionization history
double Fisher::getXE(double z, CosmParam *c)
{
   double fHe(0.08);
   double xe;
   double y, yri, deltay;

   //cout<<c->zri<<"\t"<<c->delta_zri<<endl;

   fHe=c->yHe/4.0/(1.0-c->yHe);
   y=pow(1.0+z,1.5);
   yri=pow(1.0+c->zri,1.5);
   deltay=1.5*pow(1.0+c->zri,0.5)*c->delta_zri;

   xe=(1.0+fHe)/2.0*(1.0-tanh((y-yri)/deltay));

   return xe;
}

//get neutral fraction based on electron fraction model
double Fisher::getXH(double z, CosmParam *c)
{
   double xe, xH;
   double fHe(0.08);

  fHe=c->yHe/4.0/(1.0-c->yHe);
   xe=getXE(z,c);
   xH=1.0-xe/(1.0+fHe);  //divide out effect of HeI ionization

   return xH;
}

//calculate optical depth to CMB using above ionization history
double Fisher::getTauCMB(CosmParam *c)
{
  double tau;
  double zmin(0.01),zmax(40);
  double tol(1.0e-4);

  setTauIntF(zmin,c,this,0);

  tau=qromb(getTauIntF,zmin,zmax,tol);

  return tau;
}

//integrand for tau integral
double setTauIntF(double z, CosmParam *c1, Fisher *ff1,int iflag)
{
  double dtaudz;
  static Fisher *ff;
  static CosmParam *c;
  double zHe(3.5);
  double xi, xe, xmean;
  double h;
  double fHe(0.08);

  if(iflag==0){
    c=c1;
    ff=ff1;
    return 0.0;
  }

  h=sqrt(c->ommhh/(1.0-c->oml-c->omk));

  fHe=c->yHe/4.0/(1.0-c->yHe);
  dtaudz=CRITDENSITY/PROTONMASS*(1.0-c->yHe)*c->ombhh*pow(1.0+z,3.0); //nH
  dtaudz*=SPEEDOFLIGHT_CGS*SIGMATHOMSON;
  dtaudz/=(1.0+z)*ff->hubbleZ(z,c)*UNH;
  dtaudz*=ff->getXE(z,c);
  if(z<zHe) dtaudz*=(1.0+2.0*fHe)/(1.0+fHe);  //account for Helium reionization

  return dtaudz;
}

//dummy function for tau integral
double getTauIntF(double z)
{
  return setTauIntF(z,NULL,NULL,1);
}

//set reionization parameters independently
void Fisher::setZRI(double zri, double delta_zri)
{
   fiducialFC.zri=zri;
   fiducialFC.delta_zri=delta_zri;

   return;
}


///////////////////////////////////////////////////////////////////////
// File I/O for fisher matrix
///////////////////////////////////////////////////////////////////////

// save fisher and ifisher to files
//
// Filename = exp_name + "_fisher" + tag
void Fisher::saveFisher(string tag)
{
	ofstream fout;
	string file, base;
	int i,j;

	//first output parameter names and values
	base="_param";
	file=exp_name+base+tag;
	fout.open(file.c_str());
		fout<<nparamFC<<endl;
	for(i=1;i<=nparamFC;i++){
		fout<<param_tags[i]<<"\t"<<param_values[i]<<"\t"<<sqrt(ifisherFC[i][i])<<endl;
	}
	fout.close();

	//save fiducial values
	base="_fiducial";
	file=exp_name+base+tag;
	saveCosmParam(&fiducialFC,file);

	//save fisher matrix to file
	base="_fisher";
	file=exp_name+base+tag;
	fisher_fileFC=file;
	nparamFC_file=nparamFC;
	writeMatrixS(fisherFC,nparamFC,file);

	//save ifisher matrix to file
	base="_ifisher";
	file=exp_name+base+tag;
	writeMatrixS(ifisherFC,nparamFC,file);

}

//read fisher matrix from file
void Fisher::loadFisherFile(string name, double **fisher, int nmat)
{
	ifstream fin;
	int i,j;

	fin.open(name.c_str());
	for(i=1;i<=nmat;i++){
		for(j=1;j<=nmat;j++){
			fin>>fisher[i][j];
		}
	}
	fin.close();

	return;
}

// load fisher class details from file
void Fisher::loadFisher(string name)
{
	string file, tag;
	double **fisher;
	CosmParam cosm;

	tag="_out";

	//load params
	file=name+"_param"+tag;
	loadParamValues(file);

	//read in fisher matrix
	fisher=dmatrix(1,nparamFC,1,nparamFC);
	file=name+"_fisher"+tag;
	loadFisherFile(file,fisher,nparamFC);

	//read in fiducial values
	file=name+"_fiducial"+tag;
	cosm=loadCosmParam(file);

	setFisherMatrix(fisher,nparamFC,&cosm);

	free_dmatrix(fisher,1,nparamFC,1,nparamFC);
}


// load fisher variables from file
//
// Filename = exp_name + "_fisher" + tag
void Fisher::loadParamValues(string file)
{
	ifstream fin;
	int i;
	string tag;
	double param;
	string err;

	param_tags.clear();
	param_values.clear();

	//first output parameter names and values
	fin.open(file.c_str());
		fin>>nparamFC;
		param_tags.push_back("_fiducial");
		param_values.push_back(0.0);
	for(i=1;i<=nparamFC;i++){
		fin>>tag>>param>>err;
		param_tags.push_back(tag);
		param_values.push_back(param);
	}
	fin.close();

}

//force new cosmology - use with care
void Fisher::setCosmology(CosmParam fiducial)
{
   fiducialFC=fiducial;
}

// Specify parameters to calculate Fisher matrix for
//
// There are 12 parameter choices so that choiceVec[12]
//
void Fisher::specifyParams(int *choiceVec, CosmParam fiducial)
{
	int i,nmax;
	vector<string> tags;

	fiducialFC=fiducial;

	param_tags.clear();
	param_values.clear();

	//parameter list to choose from
	tags.push_back("_fiducial");   //0
	tags.push_back("_ommhh");      //1
	tags.push_back("_ombhh");      //2
	tags.push_back("_oml");        //3
	tags.push_back("_nscal");      //4
	tags.push_back("_ascal2");     //5
	tags.push_back("_w0");         //6
	tags.push_back("_tau");        //7
	tags.push_back("_ts");         //8
	tags.push_back("_omk");        //9
	tags.push_back("_yHe");	       //10
	tags.push_back("_alpha");      //11
	tags.push_back("_Mnu");	       //12
	tags.push_back("_bias");       //13
	tags.push_back("_numass1");    //14
	tags.push_back("_numass2");    //15
	tags.push_back("_numass3");    //16
	tags.push_back("_lda");        //17
	tags.push_back("_lh");         //18
	tags.push_back("_lg");         //19
	tags.push_back("_epsilon");    //20
	tags.push_back("_eta");	       //21
	tags.push_back("_xi");         //22
	tags.push_back("_w1");         //23

        //tags.push_back("_zri");      //23
        //tags.push_back("_dzri");     //24

	nmax=tags.size()-1;

	//apply the parameter choices that apply
	param_tags.push_back("_fiducial");
	param_values.push_back(0.0);

	//consistency with distances
	//if use lda, lh don't use w, oml
	if(choiceVec[17]==1 || choiceVec[18]==1){
		distance_flag=1;
		choiceVec[17]=1;
		choiceVec[18]=1;
		choiceVec[3]=0;
		choiceVec[6]=0;
	}

	//consistency with inflationary parameters
  	//either use n_s, alpha, ts   - spectral
        //or use epsilon, eta, xi    - slow roll
	//can't use mixture as repurpose parameters in CosmParam
	if(choiceVec[20]==1 || choiceVec[21]==1 || choiceVec[22]==1){
		setSlowrollFlag(1);
		choiceVec[4]=0;
		choiceVec[8]=0;
		choiceVec[11]=0;
	}else if(choiceVec[4]==1 || choiceVec[8]==1 || choiceVec[11]==1){
		setSlowrollFlag(0);
		choiceVec[20]=0;
		choiceVec[21]=0;
		choiceVec[22]=0;
	}

	if(slowroll_flag==1) cout<<"Using slowroll inflation parameters"<<endl;
	if(slowroll_flag==0) cout<<"Using spectral inflation parameters"<<endl;

	//assign parameters
	for(i=1;i<=nmax;i++){
		if(choiceVec[i]){
			param_tags.push_back(tags[i]);
			param_values.push_back(getParamFromTag(tags[i],&fiducialFC));
		}
	}
	nparamFC=param_values.size()-1;
}

// Use string tag to assign parameter value from CosmParam structure
//
double Fisher::getParamFromTag(string tag, CosmParam *c)
{
  int i;

  if(tag=="_fiducial") return 0.0;
  if(tag=="_ommhh") return c->ommhh;
  if(tag=="_ombhh") return c->ombhh;
  if(tag=="_oml") return c->oml;
  if(tag=="_nscal") return c->nscal;
  if(tag=="_ascal2") return c->ascal2;
  if(tag=="_w0") return c->w0;
  if(tag=="_w1") return c->w1;
  if(tag=="_tau") return c->tau;
  if(tag=="_ts") return c->ts;
  if(tag=="_omk") return c->omk;
  if(tag=="_h") return c->h;
  if(tag=="_bias") return c->bias;
  if(tag=="_omb") return c->omb;
  if(tag=="_omm") return c->omm;
  if(tag=="_omnu") return c->omnu;
  if(tag=="_alpha") return c->alpha;
  if(tag=="_omnuhh") return c->omnuhh;
  if(tag=="_yHe") return c->yHe;
  if(tag=="_Mnu") return mnuFromOmnuhh(c->omnuhh);
  if(tag=="_numass1") return c->numass1;
  if(tag=="_numass2") return c->numass2;
  if(tag=="_numass3") return c->numass3;
  if(tag=="_lda") return c->lda;
  if(tag=="_lh") return c->lh;
  if(tag=="_lg") return c->lg;
  if(tag=="_epsilon") return c->ts;    //double use of these parameters
  if(tag=="_eta") return c->nscal;
  if(tag=="_xi") return c->alpha;
  if(tag=="_ntens") return c->ntens;
  if(tag=="_omchh") return c->ommhh-c->ombhh-c->omnuhh;
  if(tag=="_zri") return c->zri;
  if(tag=="_dzri") return c->delta_zri;

  //nonstandard label search in param_values
  for(i=0;i<=nparamFC;i++){
	if(param_tags[i].find(tag,0)) return param_values[i];
  }

	cout<<"Error tag not known: "<<tag<<endl;
	return -1.0;
}

// Add parameter and value to vectors
//
void Fisher::addParam(string tag)
{
	int i;

	for(i=0;i<=nparamFC;i++){
		if(param_tags[i]==tag) {
			cout<<"Already using parameter"<<endl;
			return;
		}
	}

	nparamFC++;
	param_tags.push_back(tag);
	param_values.push_back(getParamFromTag(tag,&fiducialFC));
	if(tag=="_omnu") omnu_flag=1;
}


// Add parameter and value to vectors
//
void Fisher::addParamWithValue(string tag, double value)
{
	int i;

	for(i=0;i<=nparamFC;i++){
		if(param_tags[i]==tag) {
			cout<<"Already using parameter"<<endl;
			return;
		}
	}

	nparamFC++;
	param_tags.push_back(tag);
	param_values.push_back(value);
	if(tag=="_omnu") omnu_flag=1;
}

// Remove parameter and value from vectors
//
void Fisher::removeParam(string tag)
{
	int i,indx;

	if(tag=="_omnu") omnu_flag=0;

	for(i=0;i<=nparamFC;i++){
		if(param_tags[i]==tag) {
			indx=i;
			param_tags.erase(param_tags.begin()+indx, param_tags.begin()+indx+1);
			param_values.erase(param_values.begin()+indx, param_values.begin()+indx+1);
			nparamFC--;
			return;
		}
	}

	cout<<"parameter not found"<<endl;
}

// Calculate the energy density in massive neutrinos omnuhh
// from total neutrino mass Mnu (in eV)
double Fisher::omnuhhFromMnu(double Mnu)
{
	return Mnu/93.14;
}

double Fisher::mnuFromOmnuhh(double omnuhh)
{
	return omnuhh*93.14;
}


double Fisher::getErrorFromTag(string tag)
{
	int i,found_flag(0);

	for(i=1;i<=param_tags.size()-1;i++){
           //if(param_tags[i]==tag){   //more restrictive
           if(param_tags[i].find(tag,0)!=string::npos){
			 found_flag=1;
			 return sqrt(ifisherFC[i][i]);
		}
	}

	if(found_flag==0) cout<<"Tag not found in getErrorFromTag \t"<<tag<<endl;
	return 0.0;
}

double Fisher::getFisherFromTag(string tag)
{
	int i,found_flag(0);

	for(i=1;i<=param_tags.size()-1;i++){
		if(param_tags[i]==tag){
			 found_flag=1;
			 return fisherFC[i][i];
		}
	}

	if(found_flag==0) cout<<"Tag not found in getErrorFromTag"<<endl;
	return 0.0;
}

///////////////////////////////////////////////////////////////////////
// Translate between formalisms
//////////////////////////////////////////////////////////////////////

//NEED TO THINK THROUGH THIS CONVERSION PROCESS IN MORE DETAIL
void Fisher::CosmParamFromCosmology(CosmParam cosm, Cosmology c, string name)
{
	string cltag="_cl.dat";
	string paramtag="_param.dat";
	string file;

	cosm.physical_flag=1;
	cosm.h=c.getH();
	cosm.omb=c.getOmegab();
	cosm.omm=c.getOmega0();
	cosm.oml=c.getLambda0();
	cosm.omk=1.0-cosm.omm-cosm.oml;
	cosm.omnu=0.0;  //no massive neutrinos
	cosm.nscal=c.getnSpec();
	cosm.alpha=0.0;  //no spectral running
	cosm.ts=0.0;  //no tensor modes
	cosm.tau=0.0; //need astrophysics to get this value
	//normalisations used are somewhat different, since Cosmology uses sigma8
	//while CosmParam (and CAMB) use A_S^2
	cosm.ascal2=1.0;
	cosm.sigma8=1.0;

	//cosmological constant
	cosm.w0= -1.0;
	cosm.w1=0.0;

	//derived quantities
	cosm.ombhh=cosm.omb*cosm.h*cosm.h;
	cosm.ommhh=cosm.omm*cosm.h*cosm.h;
	cosm.omnuhh=cosm.omnu*cosm.h*cosm.h;

	//file names
	file=name+cltag;
	file=name+paramtag;

}

double Fisher::daratio(CosmParam *c)
{
  double da, daFid, z;

	z=c->exp_z;
	da=angDiamDistance(z,c)/c->h;
	daFid=angDiamDistance(z,&fiducialFC)/fiducialFC.h;
	return da/daFid;
}

double Fisher::hratio(CosmParam *c)
{
  double Hc, HFid,z;

	z=c->exp_z;

	Hc=hubbleZ(z,c)*c->h;
	HFid=hubbleZ(z,&fiducialFC)*fiducialFC.h;
	return Hc/HFid;
}

////////////////////////////////////////////////////////////////////////
// Save/Load Cosmology files
/////////////////////////////////////////////////////////////////////////

double Fisher::saveCosmParam(CosmParam *c, string file)
{
  ofstream fout;

  fout.open(file.c_str());
  fout.precision(NPRECISION);

  fout <<c->physical_flag<<endl;
  fout <<c->h<<endl;
  fout <<c->omb<<endl;
  fout <<c->omm<<endl;
  fout <<c->omk<<endl;
  fout <<c->oml<<endl;
  fout <<c->omnu<<endl;
  fout <<c->ombhh<<endl;
  fout <<c->ommhh<<endl;
  fout <<c->omnuhh<<endl;
  fout <<c->tau<<endl;
  fout <<c->ts<<endl;
  fout <<c->nscal<<endl;
  fout <<c->alpha<<endl;
  fout <<c->ascal2<<endl;
  fout <<c->yHe<<endl;
  fout <<c->sigma8<<endl;
  fout <<c->bias<<endl;
  fout <<c->w0<<endl;
  fout <<c->w1<<endl;
  fout <<c->name<<endl;
  fout <<c->file<<endl;
  fout <<c->cl_file<<endl;
  fout <<c->p_step<<endl;
  fout <<c->filebase<<endl;
  fout <<c->exp_z<<endl;
  fout <<c->numass1<<endl;
  fout <<c->numass2<<endl;
  fout <<c->numass3<<endl;
  fout <<c->lda<<endl;
  fout <<c->lh<<endl;
  fout <<c->lg<<endl;
  fout <<c->Bm0<<endl;
  fout <<c->sigma0<<endl;
  fout <<c->alpha<<endl;
  fout <<c->beta<<endl;
  fout <<c->fR0<<endl;
  fout.close();

  return 0.0;
}

CosmParam Fisher::loadCosmParam(string file)
{
  ifstream fin;
  CosmParam c;

  fin.open(file.c_str());

  fin >>c.physical_flag;
  fin >>c.h;
  fin >>c.omb;
  fin >>c.omm;
  fin >>c.omk;
  fin >>c.oml;
  fin >>c.omnu;
  fin >>c.ombhh;
  fin >>c.ommhh;
  fin >>c.omnuhh;
  fin >>c.tau;
  fin >>c.ts;
  fin >>c.nscal;
  fin >>c.alpha;
  fin >>c.ascal2;
  fin >>c.yHe;
  fin >>c.sigma8;
  fin >>c.bias;
  fin >>c.w0;
  fin >>c.w1;
  fin >>c.name;
  fin >>c.file;
  fin >>c.cl_file;
  fin >>c.p_step;
  fin >>c.filebase;
  fin >>c.exp_z;
  fin >>c.numass1;
  fin >>c.numass2;
  fin >>c.numass3;
  fin >>c.lda;
  fin >>c.lh;
  fin >>c.lg;
  fin >>c.Bm0;
  fin >>c.sigma0;
  fin >>c.alpha;
  fin >>c.beta;
  fin >>c.fR0;

  fin.close();

  return c;
}

////////////////////////////////////////////////////////////////////////
// concatenation routine
////////////////////////////////////////////////////////////////////////
//concatenate two strings str1 and str2
char *concat(char *str1, char *str2){
  char *str3;

  str3 = (char *)calloc(strlen(str1) + strlen(str2) + 1,
                        sizeof(char));

  strcpy(str3, str1);
  strcat(str3, str2);

  return str3;

  //I worry about not cleaning up properly here
  free(str3);

}

/*
///////////////////////////////////////////////////////////////////
// Matrix manipulation routines
///////////////////////////////////////////////////////////////////
//return inverse of matrix leaving original untouched
//
//Warning: this routine has problems when the matrix is nearly singular
//this is a common situation when two parameters are nearly degenerate
void invertMatrix(double **fisher,int nparam, double **ifisher)
{
  int i,j,m(1);
  double **b;

  b=dmatrix(1,nparam,1,m);

  for(i=1;i<=nparam;i++){
    for(j=1;j<=nparam;j++){
      ifisher[i][j]=fisher[i][j];
    }
  }

  for(i=1;i<=nparam;i++){
    b[i][1]=1.0;
  }
  //need this routine from NR
  gaussj(ifisher,nparam,b,m);

  free_dmatrix(b,1,nparam,1,m);
}

//save a matrix to file
void writeMatrix(double **M,int nparam, char *file)
{
  ofstream fout;
  int i,j;

  fout.open(file);
  fout.precision(NPRECISION);
  for(i=1;i<=nparam;i++){
    for(j=1;j<=nparam;j++){
      if(j==nparam){
	fout<<M[i][j]<<endl;
      }else{
	fout<<M[i][j]<<"\t";
      }
    }
  }

  fout.close();
}

//save a matrix to file
void writeMatrixS(double **M,int nparam, string file)
{
  ofstream fout;
  int i,j;

  fout.open(file.c_str());
  fout.precision(NPRECISION);
  for(i=1;i<=nparam;i++){
    for(j=1;j<=nparam;j++){
      if(j==nparam){
	fout<<M[i][j]<<endl;
      }else{
	fout<<M[i][j]<<"\t";
      }
    }
  }

  fout.close();
}

////////////////////////////////////////////////////////////////////////
// Numerical recipes gaussj routine ///////////////////////////////////
///////////////////////////////////////////////////////////////////////

#define NRANSI
#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}

void gaussj(double **a, int n, double **b, int m)
{
	int *indxc,*indxr,*ipiv;
	int i,icol,irow,j,k,l,ll;
	double big,dum,pivinv,temp;

	indxc=ivector(1,n);
	indxr=ivector(1,n);
	ipiv=ivector(1,n);
	for (j=1;j<=n;j++) ipiv[j]=0;
	for (i=1;i<=n;i++) {
		big=0.0;
		for (j=1;j<=n;j++)
			if (ipiv[j] != 1)
				for (k=1;k<=n;k++) {
					if (ipiv[k] == 0) {
						if (fabs(a[j][k]) >= big) {
							big=fabs(a[j][k]);
							irow=j;
							icol=k;
						}
					}
				}
		++(ipiv[icol]);
		if (irow != icol) {
			for (l=1;l<=n;l++) SWAP(a[irow][l],a[icol][l])
			for (l=1;l<=m;l++) SWAP(b[irow][l],b[icol][l])
		}
		indxr[i]=irow;
		indxc[i]=icol;
		if (a[icol][icol] == 0.0) nrerror("gaussj: Singular Matrix");
		pivinv=1.0/a[icol][icol];
		a[icol][icol]=1.0;
		for (l=1;l<=n;l++) a[icol][l] *= pivinv;
		for (l=1;l<=m;l++) b[icol][l] *= pivinv;
		for (ll=1;ll<=n;ll++)
			if (ll != icol) {
				dum=a[ll][icol];
				a[ll][icol]=0.0;
				for (l=1;l<=n;l++) a[ll][l] -= a[icol][l]*dum;
				for (l=1;l<=m;l++) b[ll][l] -= b[icol][l]*dum;
			}
	}
	for (l=n;l>=1;l--) {
		if (indxr[l] != indxc[l])
			for (k=1;k<=n;k++)
				SWAP(a[k][indxr[l]],a[k][indxc[l]]);
	}
	free_ivector(ipiv,1,n);
	free_ivector(indxr,1,n);
	free_ivector(indxc,1,n);
}
#undef SWAP
#undef NRANSI

/////////////////////////////////////////////////////////////////
*/

/**********************************************************************
 ********************* Utilities              *************************
 *********************************************************************/
//function to call CAMB with particular parameters
//creates a .ini file and calls CAMB from the command line.
vector<double> Fisher::callCAMB(CosmParam *c){

  static int ncount;
  int precision_flag(1);
  char *param_file="./CAMB/useparam.ini";
  double sigma8;
  ifstream fin;
  ofstream fout;
  int fyes;
  Neutrinos NuClass;
  double Mnu;
  vector<double> nuFrac;
  int nu_degenerate(0);
  string slowroll_call;

  double Omb,Omc,Oml,Omnu,Omm;
  double Ombhh, Omchh,Om0hh;
  double Omk,Omnuhh,w0,w1;
  double h,hubble,tau,TS,nscal,ntens;
  double alpha;
  string matterpower_file="matterpower.dat";
  double ascal2;
  int phys_flag;
  double matter_z;
  double yHe(0.24);
  double zri(11.0), delta_zri(1.5);
  int number_massive_nu(0);
  vector<double> sigma8list;
  int use_tanh_reionization_flag(0);

  string transfer_file, cl_file;

  transfer_file="transfer_out.dat";
  matterpower_file=c->file;
  cl_file=c->cl_file;
  phys_flag=c->physical_flag;
  ncount++;
  cout <<"ncount= "<<ncount<<endl;

  //specify list of redshifts to calculate matter power spectrum
  string matter_file,tag;
  vector<string> redshift_matter_files;
  vector<string> redshift_transfer_files;
  double zuse;
  int i;

  //redshift_list runs from high to low
	cout<<matterpower_file<<endl;
	if(matterpower_file.find("_p.dat")!=std::string::npos) tag="_p.dat";
	if(matterpower_file.find("_m.dat")!=std::string::npos) tag="_m.dat";
	cout<<tag<<endl;

  if(redshift_list.size()>0){
    for(i=0;i<=redshift_list.size()-1;i++){
	zuse=redshift_list[i];
	matter_file=dirbase+"_z"+numToString(zuse)+c->filebase+tag;
	transfer_file=dirbase+"_z"+numToString(zuse)+"_transfer"+c->filebase+tag;
	redshift_matter_files.push_back(matter_file);
	redshift_transfer_files.push_back(transfer_file);

    }
  }

  //assign variables from CosmParam structure
  if(phys_flag==1){
    //input parameters
    Ombhh=c->ombhh;
    Om0hh=c->ommhh;
    if(omnu_flag==1){
	    Omnuhh=c->omnu*sqrt(c->ommhh/(1.0-c->oml-c->omk));
    }else{
	    Omnuhh=c->omnuhh;
    }
    Oml=c->oml;
    Omk=c->omk;
    tau=c->tau;
    TS=c->ts;
    nscal=c->nscal;
    ascal2=c->ascal2*1.0e-10;
    alpha=c->alpha;
    w0=c->w0;
    w1=c->w1;
    //derived parameters
    Omm=1.0-Oml-Omk;
    h=sqrt(Om0hh/Omm);
    hubble=100.0*h;
    Omchh=Om0hh-Ombhh-Omnuhh;
    Omb=Ombhh/h/h;
    Omc=Omchh/h/h;
    Omnu=Omnuhh/h/h;
	matter_z=c->exp_z;
    yHe=c->yHe;
  }else{
    Omb=c->omb;
    Omm=c->omm;
    Omc=c->omm-Omb-c->omnu;
    Oml=c->oml;
    Omk=c->omk;
    Omnu=c->omnu;
    h=c->h;
    Ombhh=Omb*h*h;
    Om0hh=c->omm*h*h;
    Omchh=Omc*h*h;
    Omnuhh=c->omnu*h*h;
    hubble=100.0*h;
    tau=c->tau;
    TS=c->ts;
    nscal=c->nscal;
    ascal2=c->ascal2*1.0e-10;
    alpha=c->alpha;
    w0=c->w0;
    matter_z=c->exp_z;
    yHe=c->yHe;
  }

  //for neutrinos need to specify mass splittings
  if(Omnuhh>0.0){
	 Mnu=NuClass.mnuFromOmnuhh(Omnuhh);
	 NuClass.initNeutrinos(Mnu,3,3,hierarchy_flag);  // normal hierarchy
	 if(Mnu<NuClass.minimumNeutrinoMass()){
		nu_degenerate=1;
		number_massive_nu=3;
		nuFrac.clear();
		nuFrac.push_back(1);
	}else{
		 number_massive_nu=3;
		// nuFrac=NuClass.massSplittings();
		nuFrac.clear();
		nuFrac.push_back(c->numass1/Mnu);
		nuFrac.push_back(c->numass2/Mnu);
		nuFrac.push_back(c->numass3/Mnu);

		//if we're varying Mnu reset splittings to match
		if(c->filebase=="_Mnu" || c->filebase=="_omnu"){
		 nuFrac.clear();
		 nuFrac=NuClass.massSplittings();
		}

		if(neutrino_mass_flag_alt==1){
		//funky parametreziation
		nuFrac.clear();
		nuFrac.push_back(c->numass1/Mnu);
		nuFrac.push_back(c->numass2/Mnu);
		nuFrac.push_back(c->numass3/Mnu);
		}

	}

  }

   w1=0.0;
   ntens=0.0;

   //set reionization parameters - only actually used if reionization_flag==1
   zri=c->zri;
   delta_zri=c->delta_zri;

  //some parameters require calculation to high precision as small derivatives
  if(c->filebase=="_numass1") precision_flag=2;
  if(c->filebase=="_numass2") precision_flag=2;
  if(c->filebase=="_numass3") precision_flag=2;

  //output parameters to see they make sense
  if(phys_flag==1){
     cout<<Ombhh<<"\t"<<Om0hh<<"\t"<<Oml<<"\t"<<Omk<<"\t"<<hubble<<"\t"<<ascal2<<"\t"<<nscal<<"\t"<<alpha<<"\t"<<TS<<endl;
  }else{
    cout<<Omb<<"\t"<<Omc<<"\t"<<Oml<<"\t"<<hubble<<endl;
  }

  //if using slowroll parameters need to precompute primordial ps
  if(slowroll_flag==1){
  	param_file="./CAMB_inf/useparam.ini";
	slowroll_call="./CAMB_inf/hsr.x ps.dat -10 7 0.01";
	slowroll_call=slowroll_call+" "+numToStringLong(ascal2,8);
	slowroll_call=slowroll_call+" "+numToStringLong(TS,8);      //epsilon
	slowroll_call=slowroll_call+" "+numToStringLong(nscal,8);	//eta
	slowroll_call=slowroll_call+" "+numToStringLong(alpha,8);	//xi
	slowroll_call=slowroll_call+" "+numToStringLong(0,8);	//lambda4
	cout<<slowroll_call<<endl;
	system(slowroll_call.c_str());
  }


  //write parameters to a CAMB .ini file
  fout.open(param_file);

fout<<"#Parameters for CAMB"<<endl;
fout<<endl;
fout<<"#output_root is prefixed to output file names"<<endl;
fout<<"output_root = "<<endl;
fout<<endl;
fout<<"#What to do"<<endl;
fout<<"get_scalar_cls = T"<<endl;
fout<<"get_vector_cls = F"<<endl;
fout<<"get_tensor_cls = T"<<endl;
fout<<"get_transfer   = T"<<endl;
fout<<endl;
fout<<"#if do_lensing then scalar_output_file contains additional columns of l^4 C_l^{pp} and l^3 C_l^{pT}"<<endl;
fout<<"#where p is the projected potential. Output lensed CMB Cls (without tensors) are in lensed_output_file below."<<endl;

if(lensing_flag==1){
fout<<"do_lensing     = T"<<endl;
}else{
fout<<"do_lensing     = F"<<endl;
}

fout<<endl;
fout<<"# 0: linear, 1: non-linear matter power (HALOFIT), 2: non-linear CMB lensing (HALOFIT)"<<endl;
fout<<"do_nonlinear = 0"<<endl;
fout<<endl;
fout<<"#Maximum multipole and k*eta. "<<endl;
fout<<"#  Note that C_ls near l_max are inaccurate (about 5%), go to 50 more than you need"<<endl;
fout<<"#  Lensed power spectra are computed to l_max_scalar-250 where accurate at %-level"<<endl;
fout<<"#  For high accuracy lensed spectra set l_max_scalar = (l you need) + 500"<<endl;
fout<<"#  To get accurate lensed BB need to have l_max_scalar>2000, k_eta_max_scalar > 10000"<<endl;
fout<<"#  Otherwise k_eta_max_scalar=2*l_max_scalar usually suffices"<<endl;
fout<<"l_max_scalar      = 3000"<<endl;
fout<<"k_eta_max_scalar  = 12000"<<endl;
fout<<endl;
fout<<"#  Tensor settings should be less than or equal to the above"<<endl;
fout<<"l_max_tensor      = 3000"<<endl;
fout<<"k_eta_max_tensor  = 8000"<<endl;
fout<<endl;
fout<<"#Main cosmological parameters, neutrino masses are assumed degenerate"<<endl;
fout<<"# If use_phyical set phyiscal densities in baryone, CDM and neutrinos + Omega_k"<<endl;
  if(phys_flag==1){
fout<<"use_physical   = T"<<endl;
fout<<"ombh2          = "<<Ombhh<<endl;
fout<<"omch2          = "<<Omchh<<endl;
fout<<"omnuh2         = "<<Omnuhh<<endl;
fout<<"omk            = "<<Omk<<endl;
fout<<"hubble         = "<<hubble<<endl;
fout<<"#effective equation of state parameter for dark energy, assumed constant"<<endl;
fout<<"w              = "<<w0<<endl;
fout<<"#constant comoving sound speed of the dark energy (1=quintessence)"<<endl;
fout<<"cs2_lam        = 1"<<endl;
fout<<endl;
fout<<"#if use_physical = F set parameters as here"<<endl;
fout<<"#omega_baryon   = "<<Omb<<endl;
fout<<"#omega_cdm      = "<<Omc<<endl;
fout<<"#omega_lambda   = "<<Oml<<endl;
fout<<"#omega_neutrino = "<<Omnu<<endl;
  }else{
fout<<"use_physical   = F"<<endl;
fout<<"#ombh2          = "<<Ombhh<<endl;
fout<<"#omch2          = "<<Omchh<<endl;
fout<<"#omnuh2         = "<<Omnuhh<<endl;
fout<<"omk            = "<<Omk<<endl;
fout<<"hubble         = "<<hubble<<endl;
fout<<"#effective equation of state parameter for dark energy, assumed constant"<<endl;
fout<<"w              = "<<w0<<endl;
fout<<"#constant comoving sound speed of the dark energy (1=quintessence)"<<endl;
fout<<"cs2_lam        = 1"<<endl;
fout<<endl;
fout<<"#if use_physical = F set parameters as here"<<endl;
fout<<"omega_baryon   = "<<Omb<<endl;
fout<<"omega_cdm      = "<<Omc<<endl;
fout<<"omega_lambda   = "<<Oml<<endl;
fout<<"omega_neutrino = "<<Omnu<<endl;
}

fout<<endl;
fout<<"#massless_neutrinos is the effective number (for QED + non-instantaneous decoupling)"<<endl;
fout<<"temp_cmb           = 2.726"<<endl;
fout<<"helium_fraction    = "<<yHe<<endl;
fout<<"massless_neutrinos = 3.04"<<endl;
fout<<"massive_neutrinos  = "<<number_massive_nu<<endl;
fout<<endl;

//handle neutrino mass splittings
if(number_massive_nu>1){
fout<<"#Neutrino mass splittings"<<endl;
fout<<"nu_mass_eigenstates = "<<NuClass.getNoMassStates()<<endl;
fout<<"#nu_mass_degeneracies = 0 sets nu_mass_degeneracies = massive_neutrinos"<<endl;
fout<<"#otherwise should be an array"<<endl;
fout<<"#e.g. for 3 neutrinos with 2 non-degenerate eigenstates, nu_mass_degeneracies = 2 1"<<endl;
fout<<"nu_mass_degeneracies =";
	if(nu_degenerate==1){
		fout<<0;
	}else{
		for(i=0;i<=nuFrac.size()-1;i++) fout<<1<<" ";
	}
fout<<endl;
fout<<"#Fraction of total omega_nu h^2 accounted for by each eigenstate, eg. 0.5 0.5"<<endl;
fout<<"nu_mass_fractions = ";
	if(nu_degenerate==1){
		fout<<1;
	}else{
		for(i=0;i<=nuFrac.size()-1;i++) fout<<nuFrac[i]<<" ";
	}
fout<<endl;
}else{
fout<<"#Neutrino mass splittings"<<endl;
fout<<"nu_mass_eigenstates = 1"<<endl;
fout<<"#nu_mass_degeneracies = 0 sets nu_mass_degeneracies = massive_neutrinos"<<endl;
fout<<"#otherwise should be an array"<<endl;
fout<<"#e.g. for 3 neutrinos with 2 non-degenerate eigenstates, nu_mass_degeneracies = 2 1"<<endl;
fout<<"nu_mass_degeneracies = 0  "<<endl;
fout<<"#Fraction of total omega_nu h^2 accounted for by each eigenstate, eg. 0.5 0.5"<<endl;
fout<<"nu_mass_fractions = 1"<<endl;
}

fout<<endl;
fout<<"#Initial power spectrum, amplitude, spectral index and running. Pivot k in Mpc^{-1}."<<endl;
fout<<"initial_power_num         = 1"<<endl;
fout<<"pivot_scalar              = 0.05"<<endl;
fout<<"pivot_tensor              = 0.05"<<endl;
fout<<"scalar_amp(1)             = "<<ascal2<<endl;
fout<<"scalar_spectral_index(1)  = "<<nscal<<endl;
fout<<"scalar_nrun(1)            = "<<alpha<<endl;
fout<<"tensor_spectral_index(1)  = "<<ntens<<endl;
fout<<"#ratio is that of the initial tens/scal power spectrum amplitudes"<<endl;
fout<<"initial_ratio(1)          = "<<TS<<endl;
fout<<"#note vector modes use the scalar settings above"<<endl;
fout<<endl;

fout<<endl;
fout<<"#Reionization, ignored unless reionization = T, re_redshift measures where x_e=0.5"<<endl;
fout<<"reionization         = T"<<endl;
fout<<endl;

if(reionization_flag==1){
fout<<"re_use_optical_depth = F"<<endl;
fout<<"re_optical_depth     ="<<tau<<endl;
fout<<"#If re_use_optical_depth = F then use following, otherwise ignored"<<endl;
fout<<"re_redshift          = "<<zri<<endl;
fout<<"#width of reionization transition. CMBFAST model was similar to re_delta_redshift~0.5."<<endl;
fout<<"re_delta_redshift    = "<<delta_zri<<endl;
fout<<"#re_ionization_frac=-1 sets to become fully ionized using YHe to get helium contribution"<<endl;
fout<<"#Otherwise x_e varies from 0 to re_ionization_frac"<<endl;
fout<<"re_ionization_frac   = -1"<<endl;
}else{
fout<<"re_use_optical_depth = T"<<endl;
fout<<"re_optical_depth     ="<<tau<<endl;
fout<<"#If re_use_optical_depth = F then use following, otherwise ignored"<<endl;
fout<<"re_redshift          = 11"<<endl;
fout<<"#width of reionization transition. CMBFAST model was similar to re_delta_redshift~0.5."<<endl;
fout<<"re_delta_redshift    = 1.5"<<endl;
fout<<"#re_ionization_frac=-1 sets to become fully ionized using YHe to get helium contribution"<<endl;
fout<<"#Otherwise x_e varies from 0 to re_ionization_frac"<<endl;
fout<<"re_ionization_frac   = -1"<<endl;
}

fout<<endl;
fout<<endl;
fout<<"#RECFAST 1.4 recombination parameters"<<endl;
fout<<"RECFAST_fudge = 1.14"<<endl;
fout<<"RECFAST_fudge_He = 0.86"<<endl;
fout<<"RECFAST_Heswitch = 6"<<endl;
fout<<endl;
fout<<endl;
fout<<"#Initial scalar perturbation mode (adiabatic=1, CDM iso=2, Baryon iso=3, "<<endl;
fout<<"# neutrino density iso =4, neutrino velocity iso = 5) "<<endl;
fout<<"initial_condition   = 1"<<endl;
fout<<"#If above is zero, use modes in the following (totally correlated) proportions"<<endl;
fout<<"#Note: we assume all modes have the same initial power spectrum"<<endl;
fout<<"initial_vector = -1 0 0 0 0"<<endl;
fout<<endl;
fout<<"#For vector modes: 0 for regular (neutrino vorticity mode), 1 for magnetic"<<endl;
fout<<"vector_mode = 0"<<endl;
fout<<endl;
fout<<"#Normalization"<<endl;
fout<<"COBE_normalize = F"<<endl;
fout<<"##CMB_outputscale scales the output Cls"<<endl;
fout<<"#To get MuK^2 set realistic initial amplitude (e.g. scalar_amp(1) = 2.3e-9 above) and"<<endl;
fout<<"#otherwise for dimensionless transfer functions set scalar_amp(1)=1 and use"<<endl;
fout<<"#CMB_outputscale = 1"<<endl;
fout<<"CMB_outputscale = 7.4311e12"<<endl;
fout<<endl;


 //handle multiple redshifts for output
  if(redshift_list.size()>0){
fout<<"#Transfer function settings, transfer_kmax=0.5 is enough for sigma_8"<<endl;
fout<<"#transfer_k_per_logint=0 sets sensible non-even sampling; "<<endl;
fout<<"#transfer_k_per_logint=5 samples fixed spacing in log-k"<<endl;
fout<<"#transfer_interp_matterpower =T produces matter power in regular interpolated grid in log k; "<<endl;
fout<<"# use transfer_interp_matterpower =F to output calculated values (e.g. for later interpolation)"<<endl;
fout<<"transfer_high_precision = T"<<endl;
fout<<"transfer_kmax           = 5"<<endl;
if(precision_flag>1){
fout<<"transfer_k_per_logint   = 25"<<endl;
}else{
fout<<"transfer_k_per_logint   = 0"<<endl;
}

fout<<"transfer_num_redshifts  = "<<redshift_list.size()<<endl;
fout<<"transfer_interp_matterpower = F"<<endl;
  for(i=0;i<=redshift_list.size()-1;i++){
     fout<<"transfer_redshift("<<i+1<<")    = "<<redshift_list[i]<<endl;
     fout<<"transfer_filename("<<i+1<<")    = "<<redshift_transfer_files[i]<<endl;
     fout<<"#Matter power spectrum output against k/h in units of h^3 Mpc^{-3}"<<endl;
     fout<<"transfer_matterpower("<<i+1<<") = "<<redshift_matter_files[i]<<endl;
  }

fout<<endl;
}else{
fout<<"#Transfer function settings, transfer_kmax=0.5 is enough for sigma_8"<<endl;
fout<<"#transfer_k_per_logint=0 sets sensible non-even sampling; "<<endl;
fout<<"#transfer_k_per_logint=5 samples fixed spacing in log-k"<<endl;
fout<<"#transfer_interp_matterpower =T produces matter power in regular interpolated grid in log k; "<<endl;
fout<<"# use transfer_interp_matterpower =F to output calculated values (e.g. for later interpolation)"<<endl;
fout<<"transfer_high_precision = T"<<endl;
fout<<"transfer_kmax           = 15"<<endl;
fout<<"transfer_k_per_logint   = 0"<<endl;
fout<<"transfer_num_redshifts  = 1"<<endl;
fout<<"transfer_interp_matterpower = F"<<endl;
fout<<"transfer_redshift(1)    = 0"<<endl;
fout<<"transfer_filename(1)    = "<<transfer_file<<endl;
fout<<"#Matter power spectrum output against k/h in units of h^{-3} Mpc^3"<<endl;
fout<<"transfer_matterpower(1) = "<<matterpower_file<<endl;
fout<<endl;
}

fout<<endl;
fout<<"#Output files not produced if blank. make camb_fits to use use the FITS setting."<<endl;
fout<<"scalar_output_file = scalCls.dat"<<endl;
fout<<"vector_output_file = vecCls.dat"<<endl;
fout<<"tensor_output_file = tensCls.dat"<<endl;

if(lensing_flag==0){
fout<<"total_output_file  = "<<cl_file<<endl;
fout<<"lensed_output_file = lensedCls.dat"<<endl;
fout<<"lensed_total_output_file  =lensedtotCls.dat"<<endl;
}else{
fout<<"total_output_file  = totalCls.dat"<<endl;
fout<<"lensed_output_file = lensedCls.dat"<<endl;
fout<<"lensed_total_output_file  ="<<cl_file<<endl;
}

fout<<"FITS_filename      = scalCls.fits"<<endl;
fout<<endl;
fout<<"##Optional parameters to control the computation speed,accuracy and feedback"<<endl;
fout<<endl;
fout<<"#If feedback_level > 0 print out useful information computed about the model"<<endl;
fout<<"feedback_level = 0"<<endl;
fout<<endl;
fout<<"# 1: curved correlation function, 2: flat correlation function, 3: inaccurate harmonic method"<<endl;
fout<<"lensing_method = 1"<<endl;
fout<<"accurate_BB = T"<<endl;
fout<<endl;
fout<<endl;
fout<<"#massive_nu_approx: 0 - integrate distribution function"<<endl;
fout<<"#                   1 - switch to series in velocity weight once non-relativistic"<<endl;
fout<<"#                   2 - use fast approximate scheme (CMB only- accurate for light neutrinos)"<<endl;
fout<<"#                   3 - intelligently use the best accurate method"<<endl;
if(precision_flag>1){
fout<<"massive_nu_approx = 0"<<endl;
}else{
fout<<"massive_nu_approx = 3"<<endl;
}

fout<<endl;
fout<<"#Whether you are bothered about polarization. "<<endl;
fout<<"accurate_polarization   = T"<<endl;
fout<<endl;
fout<<"#Whether you are bothered about percent accuracy on EE from reionization"<<endl;
fout<<"accurate_reionization   = T"<<endl;
fout<<endl;
fout<<"#whether or not to include neutrinos in the tensor evolution equations"<<endl;
fout<<"do_tensor_neutrinos     = T"<<endl;
fout<<endl;
fout<<"#Whether to turn off small-scale late time radiation hierarchies (save time,v. accurate)"<<endl;
fout<<"do_late_rad_truncation   = T"<<endl;
fout<<endl;
fout<<"#Computation parameters"<<endl;
fout<<"#if number_of_threads=0 assigned automatically"<<endl;
fout<<"number_of_threads       = 0"<<endl;
fout<<endl;
fout<<"#Default scalar accuracy is about 0.3% (except lensed BB). "<<endl;
fout<<"#For 0.1%-level try accuracy_boost=2, l_accuracy_boost=2."<<endl;
fout<<endl;
fout<<"#Increase accuracy_boost to decrease time steps, use more k values,  etc."<<endl;
fout<<"#Decrease to speed up at cost of worse accuracy. Suggest 0.8 to 3."<<endl;
fout<<"accuracy_boost          ="<<precision_flag<<endl;
fout<<endl;
fout<<"#Larger to keep more terms in the hierarchy evolution. "<<endl;
fout<<"l_accuracy_boost        ="<<precision_flag<<endl;
fout<<endl;
fout<<"#Increase to use more C_l values for interpolation."<<endl;
fout<<"#Increasing a bit will improve the polarization accuracy at l up to 200 -"<<endl;
fout<<"#interpolation errors may be up to 3%"<<endl;
fout<<"#Decrease to speed up non-flat models a bit"<<endl;
fout<<"l_sample_boost          ="<<precision_flag<<endl;

  fout.close();

  if(slowroll_flag==0){
	  cout <<"Calling CAMB"<<endl;
	  //fyes=system("/usr/local/CAMB/camb /usr/local/CAMB/useparam.ini > temp.dat");
	  fyes=system("./CAMB/camb ./CAMB/useparam.ini > temp.dat");
  }else{
	  cout <<"Calling CAMB_inf"<<endl;
	  fyes=system("./CAMB_inf/camb ./CAMB_inf/useparam.ini > temp.dat");
  }
  if (fyes==0) cout <<"Call to CAMB successful"<<endl;

  //get value of sigma 8
  //Need modified version of CAMB to output this information
  fin.open("temp.dat");

  if(redshift_list.size()>0){
	  for(i=1;i<=redshift_list.size();i++){
		  fin>>sigma8;
		 cout<<redshift_list[i-1]<<"\t"<<sigma8<<endl;
		sigma8list.push_back(sigma8);
 	 }
  }else{
		  fin>>sigma8;
		  sigma8list.push_back(sigma8);
  }

  fin.close();

  cout<<sigma8<<endl;

  return sigma8list;
}

////////////////////////////////////////////////////////////////////
//  Set Fisher matrix data files
///////////////////////////////////////////////////////////////////
//
// Note: step size for oml and w0 need to be larger than other parameters
// if I decrease the step size the errors first get larger then
// become smaller and smaller.  These are probably from numerical noise
// values chosen seem to produce stable derivatives.

double Fisher::setCosmFiles()
{
  int i;
  double z;

  //set redshifts to calculate matter power spectrum

	redshift_list.clear();
	/*redshift_list.push_back(4.0);
	redshift_list.push_back(3.5);
	redshift_list.push_back(3.0);
	redshift_list.push_back(2.5);
	redshift_list.push_back(2.0);
	redshift_list.push_back(1.5);
	redshift_list.push_back(1.0);
	redshift_list.push_back(0.5);
	redshift_list.push_back(0.3);
	redshift_list.push_back(0.0);	*/

   // generate all redshift files

	redshift_list.clear();
	z=1.0;
	while(z>=0.0){
		redshift_list.push_back(z);
		z-=0.05;
	}
	redshift_list.push_back(0.0);


  //parameters for calculating derivatives
  for(i=0;i<=nparamFC;i++) setParamPair(&fiducialFC,i,1);
  return 1.0;
}

double Fisher::setParamPair(CosmParam *fiducial, int iflag, int clflag)
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
  double zcmb(1100.0);

  zcmb=getZCMB(&fiducialFC);

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
  c_p.ommhh+=p_step;
  c_m.ommhh-=p_step;
  c_p.h=sqrt(c_p.ommhh/c_p.omm);
  c_m.h=sqrt(c_m.ommhh/c_m.omm);

  }else if(param_tags[iflag]=="_ombhh"){
    file=name+param_tags[iflag];
  p_step=vary*fiducial->ombhh;
  c_p.ombhh+=p_step;
  c_m.ombhh-=p_step;

  }else if(param_tags[iflag]=="_oml"){
    file=name+param_tags[iflag];
  p_step=vary*fiducial->omm;
  c_p.oml+=p_step;
  c_m.oml-=p_step;
  c_p.omm=1.0-c_p.oml-c_p.omk;
  c_m.omm=1.0-c_m.oml-c_m.omk;
  c_p.h=sqrt(c_p.ommhh/c_p.omm);
  c_m.h=sqrt(c_m.ommhh/c_m.omm);

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

  }else if(param_tags[iflag]=="_w0" && omk_flag==0){
     file=name+param_tags[iflag];
     p_step=vary2*fiducial->w0;
     c_p.w0+=p_step;
     c_m.w0-=p_step;

  }else if(param_tags[iflag]=="_w0" && omk_flag==1){
     file=name+param_tags[iflag];
     p_step=vary2*fiducial->w0;
     c_p.w0+=p_step;
     c_m.w0-=p_step;

     //want derivatives at fixed D_A - will convert later
     c_p.oml=fixDA(&c_p,coAngDiamDistance(zcmb,fiducial));
     c_m.oml=fixDA(&c_m,coAngDiamDistance(zcmb,fiducial));

     c_p.omm=1.0-c_p.oml-c_p.omk;   //ensure consistency
     c_m.omm=1.0-c_m.oml-c_m.omk;   //ensure consistency
     c_p.h=sqrt(c_p.ommhh/c_p.omm);
     c_m.h=sqrt(c_m.ommhh/c_m.omm);

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

  }else if(param_tags[iflag]=="_omk" && omk_flag==0){
  //being sloppy with this as should treat curvature more carefully
    file=name+param_tags[iflag];
  p_step=vary2*fiducial->omm;
  c_p.omk+=p_step;
  c_m.omk-=p_step;
  c_p.omm=1.0-c_p.oml-c_p.omk;   //ensure consistency
  c_m.omm=1.0-c_m.oml-c_m.omk;   //ensure consistency
  c_p.h=sqrt(c_p.ommhh/c_p.omm);
  c_m.h=sqrt(c_m.ommhh/c_m.omm);

  }else if(param_tags[iflag]=="_omk" && omk_flag==1){
  //handle curvature more carefully
    file=name+param_tags[iflag];
  p_step=0.02;

  c_p.omk+=0.0;
  c_m.omk-=2.0*p_step;

  //want derivatives at fixed D_A - will convert later
  c_p.oml=fixDA(&c_p,coAngDiamDistance(zcmb,fiducial));
  c_m.oml=fixDA(&c_m,coAngDiamDistance(zcmb,fiducial));

  //c_p.oml=fixThetas(&c_p,angularScale(fiducial));
  //c_m.oml=fixThetas(&c_m,angularScale(fiducial));

  c_p.omm=1.0-c_p.oml-c_p.omk;   //ensure consistency
  c_m.omm=1.0-c_m.oml-c_m.omk;   //ensure consistency
  c_p.h=sqrt(c_p.ommhh/c_p.omm);
  c_m.h=sqrt(c_m.ommhh/c_m.omm);

  }else if(param_tags[iflag]=="_omnuhh"){
    file=name+param_tags[iflag];
  p_step=vary2*fiducial->omnuhh;
  c_p.omnuhh+=p_step;
  c_m.omnuhh-=p_step;
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
          if(param_tags[iflag]=="_omk" && omk_flag==1){
             derivCMB_OMK(lmax,fiducial,name);
          }
          if(param_tags[iflag]=="_w0" && omk_flag==1){
             derivCMB_W0(lmax,fiducial,name);
          }
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

void Fisher::derivCMB(int lmax, double p_step,CosmParam *cm, CosmParam *cp,string file)
{
  string files;
  ifstream fin_m;
  ifstream fin_p;
  ofstream fout;
  int lm,lp;
  double clTTm,clTTp;
  double clTEm,clTEp;
  double clEEm,clEEp;
  double clBBm,clBBp;
  double dclTT,dclTE,dclEE,dclBB;
  int i;

  files=cm->cl_file;
  fin_m.open(files.c_str());
  files=cp->cl_file;
  fin_p.open(files.c_str());

  fout.open(file.c_str());

  for(i=2;i<=lmax;i++){
    fin_m>>lm>>clTTm>>clEEm>>clBBm>>clTEm;
    fin_p>>lp>>clTTp>>clEEp>>clBBp>>clTEp;
    dclTT=(clTTp-clTTm)/(2.0*p_step);
    dclTE=(clTEp-clTEm)/(2.0*p_step);
    dclEE=(clEEp-clEEm)/(2.0*p_step);
    dclBB=(clBBp-clBBm)/(2.0*p_step);
    fout <<lm<<"\t"<<dclTT<<"\t"<<dclEE<<"\t"<<dclBB<<"\t"<<dclTE<<endl;
  }

  fin_m.close();
  fin_p.close();
  fout.close();
}

// Calculate derivative of matter power spectrum
//
void Fisher::derivGAL(string tag,string file)
{
	CosmParam cp, cm, c;
	ofstream fout;
	double dkmin, dkmax;
 	string base, tagm,tagp;
	double p_step;
	int i,nstep;
	double k,kstep;
	double Pp, Pm, dP, dlP, Ps, dlP2;

	Spline splineP, splineM, splineFID;

	tagm="_m.dat";
  	tagp="_p.dat";
	base=dirbase+"_z"+numToString(0.0);
	base=dirbase+"_z"+numToString(0.3);
	base=dirbase+"_z"+numToString(1.0);
//	base=dirbase+"_z"+numToString(4.0);
	base=dirbase+"_z"+numToString(8.0);


	cp=loadCosmParam(base+tag+"_param"+tagp);
	cm=loadCosmParam(base+tag+"_param"+tagm);
	c=loadCosmParam(base+"_fiducial"+"_param"+tagm);

        p_step=getParamFromTag(tag,&cp)-getParamFromTag(tag,&cm);

	assignSpline(&splineFID,base,"_fiducial"+tagm);
	assignSpline(&splineP,base,tag+tagp);
	assignSpline(&splineM,base,tag+tagm);


	dkmin=fmax(splineP.getXMin(),splineM.getXMin())*0.8;
	dkmax=fmin(splineP.getXMax(),splineM.getXMax())*0.7;
	nstep=splineP.getN();

	kstep=exp(log(dkmax/dkmin)/(double)(nstep));

	fout.open(file.c_str());
	k=dkmin;
	for(i=2;i<nstep-1;i++){
		k=splineP.getXelement(i)*0.74;
		Pp=powerSpectrum(k,&cp,&splineP);
		Pm=powerSpectrum(k,&cm,&splineM);
		Ps=powerSpectrum(k,&c,&splineFID);
		dP=(Pp-Pm)/p_step;
		dlP=(log(Pp)-log(Pm))/p_step;
		dlP2=dP/Ps;
		fout<<k<<"\t"<<Pp<<"\t"<<Pm<<"\t"<<dP<<"\t"<<dlP<<"\t"<<Ps<<"\t"<<dlP2<<endl;
		cout<<k<<"\t"<<Pp<<"\t"<<Pm<<"\t"<<dP<<"\t"<<dlP<<endl;
	}
	fout.close();

}

/////////////////////////////////////////////////////////////////
// Code for handling derivatives w.r.t curvature more carefully
/////////////////////////////////////////////////////////////////
// Convert derivative w.r.t omk at fixed da to
// derivative w.r.t. omk at fixed oml
void Fisher::derivCMB_OMK(int lmax, CosmParam *c, string name)
{
   string fileoml, fileomk;
   string file_out;
   double pstep_oml, pstep_omk;
   CosmParam ck_m, ck_p;
   CosmParam cl_m, cl_p;
   double dDAdomk, dDAdoml;
   double zcmb(1100.0);
   int clip_flag(1);

  zcmb=getZCMB(&fiducialFC);

   cout<<name<<endl;
   fileomk=name+"_omk_cl_deriv.dat";
   fileoml=name+"_oml_cl_deriv.dat";
   file_out=name+"_omkDA_cl_deriv.dat";

   ck_p=loadCosmParam(name+"_fiducial_param_p.dat");
   ck_m=ck_p;
   cl_p=ck_p;
   cl_m=ck_p;

   pstep_omk=0.001;
   ck_p.omk=ck_p.omk+pstep_omk;
   ck_m.omk=ck_m.omk-pstep_omk;
   dDAdomk=(coAngDiamDistance(zcmb,&ck_p)-coAngDiamDistance(zcmb,&ck_m))/(2.0*pstep_omk);

   pstep_oml=0.001;
   cl_p.oml=cl_p.oml+pstep_oml;
   cl_m.oml=cl_m.oml-pstep_oml;
   dDAdoml=(coAngDiamDistance(zcmb,&cl_p)-coAngDiamDistance(zcmb,&cl_m))/(2.0*pstep_oml);

   cout<<"derivatives: "<<dDAdomk<<"\t"<<dDAdoml<<endl;

  ifstream fin_l;
  ifstream fin_k;
  ofstream fout;
  int ll,lk;
  double clTTl,clTTk;
  double clTEl,clTEk;
  double clEEl,clEEk;
  double clBBl,clBBk;
  double dclTT,dclTE,dclEE,dclBB;
  int i;
  int lclip(30);
  double clipTTk(0.0), clipTEk(0.0), clipEEk(0.0), clipBBk(0.0);

  fin_k.open(fileomk.c_str());
  fin_l.open(fileoml.c_str());

  fout.open(file_out.c_str());

  for(i=2;i<=lmax;i++){
    fin_l>>ll>>clTTl>>clEEl>>clBBl>>clTEl;
    fin_k>>lk>>clTTk>>clEEk>>clBBk>>clTEk;

    //noise clipping
    if(ll==lclip && clip_flag==1){
       clipTTk=clTTk/(double)(ll*(ll+1));
       clipTEk=clTEk/(double)(ll*(ll+1));;
       clipEEk=clEEk/(double)(ll*(ll+1));;
       clipBBk=clBBk/(double)(ll*(ll+1));;
    }else if(ll>lclip && clip_flag==1){
       //have to be careful about int values here
       clTTk=clipTTk*pow((double)(lclip)/ll,4.0-c->nscal)*(double)(ll*(ll+1));
       clTEk=clipTEk*pow((double)(lclip)/ll,4.0-c->nscal)*(double)(ll*(ll+1));
       clEEk=clipEEk*pow((double)(lclip)/ll,4.0-c->nscal)*(double)(ll*(ll+1));
       clBBk=clipBBk*pow((double)(lclip)/ll,4.0-c->nscal)*(double)(ll*(ll+1));
    }

    cout<<ll<<"\t"<<clTTk<<"\t"<<clipTTk<<endl;

    dclTT=clTTk+dDAdomk/dDAdoml*clTTl;
    dclTE=clTEk+dDAdomk/dDAdoml*clTEl;
    dclEE=clEEk+dDAdomk/dDAdoml*clEEl;
    dclBB=clBBk+dDAdomk/dDAdoml*clBBl;

    fout <<lk<<"\t"<<dclTT<<"\t"<<dclEE<<"\t"<<dclBB<<"\t"<<dclTE<<endl;
    //cout <<lk<<"\t"<<dclTT<<"\t"<<dclEE<<"\t"<<dclBB<<"\t"<<dclTE<<endl;
  }

  fin_k.close();
  fin_l.close();
  fout.close();

  //now overwrite the file containing derivatives w.r.t omk
  string savestr;
  savestr="cp "+fileomk+" "+name+"_omkDAs_cl_deriv.dat";
  system(savestr.c_str());

  fin_k.open(file_out.c_str());
  fout.open(fileomk.c_str());
  for(i=2;i<=lmax;i++){
    fin_k>>lk>>clTTk>>clEEk>>clBBk>>clTEk;
    fout<<lk<<"\t"<<clTTk<<"\t"<<clEEk<<"\t"<<clBBk<<"\t"<<clTEk<<endl;
  }
  fin_k.close();
  fout.close();

}

// Convert derivative w.r.t omk at fixed da to
// derivative w.r.t. omk at fixed oml
void Fisher::derivCMB_W0(int lmax, CosmParam *c, string name)
{
   string fileoml, filew0;
   string file_out;
   double pstep_oml, pstep_w0;
   CosmParam ck_m, ck_p;
   CosmParam cl_m, cl_p;
   double dDAdw0, dDAdoml;
   double zcmb(1100.0);
   int clip_flag(1);

  zcmb=getZCMB(&fiducialFC);

   cout<<name<<endl;
   filew0=name+"_w0_cl_deriv.dat";
   fileoml=name+"_oml_cl_deriv.dat";
   file_out=name+"_w0DA_cl_deriv.dat";

   ck_p=loadCosmParam(name+"_fiducial_param_p.dat");
   ck_m=ck_p;
   cl_p=ck_p;
   cl_m=ck_p;

   pstep_w0=0.01;
   ck_p.w0=ck_p.w0+pstep_w0;
   ck_m.w0=ck_m.w0-pstep_w0;
   dDAdw0=(coAngDiamDistance(zcmb,&ck_p)-coAngDiamDistance(zcmb,&ck_m))/(2.0*pstep_w0);

   pstep_oml=0.001;
   cl_p.oml=cl_p.oml+pstep_oml;
   cl_m.oml=cl_m.oml-pstep_oml;
   dDAdoml=(coAngDiamDistance(zcmb,&cl_p)-coAngDiamDistance(zcmb,&cl_m))/(2.0*pstep_oml);

   cout<<"derivatives: "<<dDAdw0<<"\t"<<dDAdoml<<endl;

  ifstream fin_l;
  ifstream fin_k;
  ofstream fout;
  int ll,lk;
  double clTTl,clTTk;
  double clTEl,clTEk;
  double clEEl,clEEk;
  double clBBl,clBBk;
  double dclTT,dclTE,dclEE,dclBB;
  int i;
  int lclip(100);
  double clipTTk, clipTEk, clipEEk, clipBBk;

  fin_k.open(filew0.c_str());
  fin_l.open(fileoml.c_str());

  fout.open(file_out.c_str());

  for(i=2;i<=lmax;i++){
    fin_l>>ll>>clTTl>>clEEl>>clBBl>>clTEl;
    fin_k>>lk>>clTTk>>clEEk>>clBBk>>clTEk;

    //noise clipping
    if(ll==lclip && clip_flag==1){
       clipTTk=clTTk/ll/(ll+1);
       clipTEk=clTEk/ll/(ll+1);
       clipEEk=clEEk/ll/(ll+1);
       clipBBk=clBBk/ll/(ll+1);
    }else if(ll>lclip && clip_flag==1){
       clTTk=clipTTk*pow((double)(lclip)/ll,4.0-c->nscal)*ll*(ll+1);
       clTEk=clipTEk*pow((double)(lclip)/ll,4.0-c->nscal)*ll*(ll+1);
       clEEk=clipEEk*pow((double)(lclip)/ll,4.0-c->nscal)*ll*(ll+1);
       clBBk=clipBBk*pow((double)(lclip)/ll,4.0-c->nscal)*ll*(ll+1);
    }

    dclTT=clTTk+dDAdw0/dDAdoml*clTTl;
    dclTE=clTEk+dDAdw0/dDAdoml*clTEl;
    dclEE=clEEk+dDAdw0/dDAdoml*clEEl;
    dclBB=clBBk+dDAdw0/dDAdoml*clBBl;

    fout <<lk<<"\t"<<dclTT<<"\t"<<dclEE<<"\t"<<dclBB<<"\t"<<dclTE<<endl;
    //cout <<lk<<"\t"<<dclTT<<"\t"<<dclEE<<"\t"<<dclBB<<"\t"<<dclTE<<endl;
  }

  fin_k.close();
  fin_l.close();
  fout.close();

  //now overwrite the file containing derivatives w.r.t w0
  string savestr;
  savestr="cp "+filew0+" "+name+"_w0DAs_cl_deriv.dat";
  system(savestr.c_str());

  fin_k.open(file_out.c_str());
  fout.open(filew0.c_str());
  for(i=2;i<=lmax;i++){
    fin_k>>lk>>clTTk>>clEEk>>clBBk>>clTEk;
    fout<<lk<<"\t"<<clTTk<<"\t"<<clEEk<<"\t"<<clBBk<<"\t"<<clTEk<<endl;
  }
  fin_k.close();
  fout.close();

}



// Routine to return the value of oml required to obtain specifed value of D_A
// given other cosmological parameters specified in c.  Obtained by root finding
//
// Need this for calculating derivatives of C_l at fixed D_A
//
// NOT VERY ROBUST TO CHANGES IN OMK STEP SIZE < 0.02 needed to work
//
double Fisher::fixDA(CosmParam *c, double thetas)
{
	double oml;
	double thetasmin, thetasmax;
	double omlmin(0.001);
	double omlmax(0.99);
	double tol(1.0e-5);

        omlmax=omlmax-c->omk;

	setFixDA(0.0,c,this,1);

	thetasmax=getFixDA(omlmin);
	thetasmin=getFixDA(omlmax);
        cout<<"limits: "<<thetasmin<<"\t"<<thetasmax<<endl;

	//check thetas is within bounds for rootfinding
	if(thetas>thetasmax || thetas<thetasmin){
           cout<<"Error d_A outside of bounds: "<<"\t"<<thetas<<"\t"<<thetasmin<<"\t"<<thetasmax<<endl;
		return -1.0;
	}

	oml=zriddrConst(getFixDA,thetas,omlmin,omlmax,tol);
	cout<<oml<<endl;
	return oml;
}

double getFixDA(double oml)
{
	return setFixDA(oml,NULL,NULL,0);
}

double setFixDA(double oml, CosmParam *c1, Fisher *ff1, int iflag)
{
	static CosmParam c;
	static Fisher *ff;
	double da;
        double zcmb(1100.0);
	string file;

	if(iflag==1){
		ff=ff1;
		file="c_temp.dat";
		ff->saveCosmParam(c1,file);
		c=ff->loadCosmParam(file);
		return 0.0;
	}

        zcmb=ff->getZCMB(&c);
	c.oml=oml;
	da=ff->coAngDiamDistance(zcmb,&c);
    	return da;
}

////////////////////////////////////////////////////////////////////////
// Calculate angular scale of CMB - position of first peak
///////////////////////////////////////////////////////////////////////
//

//calculate CMB angular scale - effectively multipole of first peak
double Fisher::angularScale(CosmParam *c)
{
	double rs, dC, lpeak;
	double zcmb(1100.0);
        zcmb=getZCMB(c);

	rs=soundHorizon(c);
	dC=coordDistanceNum(zcmb,c);  //comoving angular diameter distance

	lpeak=PI*dC/rs;

	return lpeak;
}

//calculate sound horizon at CMB redshift for specified cosmology
//use variable x=1/z to handle integral to infinity
//
// Units :  c/H0
double Fisher::soundHorizon(CosmParam *c)
{
	double rs;
	double eta;
	double tol(1.0e-5);
	double zmin, zmax, xmin, xmax;

	setSoundHorizonKernel(0.0,c,this,1);

	zmin=1100.0;  //really this should depend on cosmology
        zmin=getZCMB(c);
	zmax=1.0e8;

	xmin=1.0/zmax;
	xmax=1.0/zmin;
	rs=qromb(getSoundHorizonKernel,xmin,xmax,tol);

	return rs;
}

double getSoundHorizonKernel(double x)
{
	return setSoundHorizonKernel(x,NULL,NULL,0);
}

double setSoundHorizonKernel(double x, CosmParam *c1, Fisher *ff1, int iflag)
{
	static Fisher *ff;
	static CosmParam *c;
	double drs, cs;
	double z, R, H;
	double omegaghh(2.47e-5);

	if(iflag==1){
		ff=ff1;
		c=c1;
		return 0.0;
	}

	z=1.0/x;

	R=3.0*c->ombhh/4.0/(omegaghh)/(1.0+z);
	cs=1.0/sqrt(3.0*(1.0+R));
	H=ff->hubbleZ(z,c);
	drs=cs/H;

	drs/=x*x;

	return drs;
}


//////////////////////
//
/////////////////////

// Routine to return the value of oml required to obtain specifed value of thetas
// given other cosmological parameters specified in c.  Obtained by root finding
//
// Need this for calculating derivatives of C_l at fixed thetas
//
double Fisher::fixThetas(CosmParam *c, double thetas)
{
	double oml;
	double thetasmin, thetasmax;
	double omlmin(0.001);
	double omlmax(0.999);
	double tol(1.0e-5);

	setFixThetas(0.0,c,this,1);

	thetasmax=getFixThetas(omlmin);
	thetasmin=getFixThetas(omlmax);

	//check thetas is within bounds for rootfinding
	if(thetas>thetasmax || thetas<thetasmin){
		cout<<"Error thetas outside of bounds: "<<"\t"<<thetasmin<<"\t"<<thetasmax<<endl;
		return -1.0;
	}

	oml=zriddrConst(getFixThetas,thetas,omlmin,omlmax,tol);

	return oml;
}

double getFixThetas(double oml)
{
	return setFixThetas(oml,NULL,NULL,0);
}

double setFixThetas(double oml, CosmParam *c1, Fisher *ff1, int iflag)
{
	static CosmParam c;
	static Fisher *ff;
	double thetas;
	string file;

	if(iflag==1){
		ff=ff1;
		file="c_temp.dat";
		ff->saveCosmParam(c1,file);
		c=ff->loadCosmParam(file);
		return 0.0;
	}

	c.oml=oml;
	thetas=ff->angularScale(&c);

	return thetas;
}

/////////////////////////////////////////////////////////////////
//  Generate error ellipse data
/////////////////////////////////////////////////////////////////

//calculate the ellipse data for my fisher matrix
//iflag=1 CMB
//iflag=2 GALAXY
//iflag=3 CMB +GAL
//iflag=4 CMB +GAL1+GAL2
//iflag=5 GAL1+GAL2
int Fisher::fisherContour()
{
  ofstream fout;
  int i,j;
  double param[12];
  double **A,**iA;
  double A11,A22,A12;
  string name;
  A=dmatrix(1,2,1,2);
  iA=dmatrix(1,2,1,2);

  double **ifisher;
  int nparam;
  CosmParam *c;

  ifisher=ifisherFC;
  nparam=nparamFC;
  c=&fiducialFC;

  name=exp_name+"_contour.dat";

  fout.open(name.c_str());
  for(i=1;i<=nparam;i++){
    for(j=1;j<=nparam;j++){
       if(j>i){
	 A[1][1]=ifisher[i][i];
	 A[1][2]=ifisher[i][j];
	 A[2][1]=ifisher[j][i];
	 A[2][2]=ifisher[j][j];
	 invertMatrix(A,2,iA);

	 A11=iA[1][1];
	 A12=iA[1][2];
	 A22=iA[2][2];
	 fout<<getLatexFromTag(param_tags[i])<<"\t"<<getLatexFromTag(param_tags[j])<<"\t"<<param_values[i]<<"\t"<<param_values[j]<<"\t"<<A11<<"\t"<<A22<<"\t"<<A12<<endl;
      }
    }
  }
  fout.close();

  free_dmatrix(A,1,2,1,2);
  free_dmatrix(iA,1,2,1,2);

  return 0;
}

//calculate the ellipse data for my fisher matrix - dont marginalise over other variables
//iflag=1 CMB
//iflag=2 GALAXY
//iflag=3 CMB +GAL
//iflag=4 CMB +GAL1+GAL2
//iflag=5 GAL1+GAL2
int Fisher::fisherContourNotMarg()
{
  ofstream fout;
  int i,j;
  double param[12];
  double **A,**iA;
  double A11,A22,A12;
  string name;
  A=dmatrix(1,2,1,2);
  iA=dmatrix(1,2,1,2);

  double **fisher;
  int nparam;
  CosmParam *c;

  fisher=fisherFC;
  nparam=nparamFC;
  c=&fiducialFC;

  name=exp_name+"_contour.dat";

  fout.open(name.c_str());
  for(i=1;i<=nparam;i++){
    for(j=1;j<=nparam;j++){
       if(j>i){
	 A[1][1]=fisher[i][i];
	 A[1][2]=fisher[i][j];
	 A[2][1]=fisher[j][i];
	 A[2][2]=fisher[j][j];
	 //invertMatrix(A,2,iA);

	 A11=A[1][1];
	 A12=A[1][2];
	 A22=A[2][2];
	 fout<<getLatexFromTag(param_tags[i])<<"\t"<<getLatexFromTag(param_tags[j])<<"\t"<<param_values[i]<<"\t"<<param_values[j]<<"\t"<<A11<<"\t"<<A22<<"\t"<<A12<<endl;
      }
    }
  }
  fout.close();

  free_dmatrix(A,1,2,1,2);
  free_dmatrix(iA,1,2,1,2);

  return 0;
}

///////////////////////////////////////////////////////////////////////
// print matrix
///////////////////////////////////////////////////////////////////////

void printMatrix(double **matrixin, int nmat)
{
	int i,j;

	cout<<endl;
	for(i=1;i<=nmat;i++){
		for(j=1;j<=nmat;j++){
			cout<<matrixin[i][j]<<"\t";
		}
	cout<<endl;
	}
	cout<<endl;

}

///////////////////////////////////////////////////////////////////////
// regulate size of fisher matrix and which parameters to include
///////////////////////////////////////////////////////////////////////

//specify ivector containing which elements of full fisher matrix
//to use in calculating errors
void Fisher::reduceFisherMatrix(int ichoice[])
{
	double **fisherNew;
	vector<string> param_tags_new;
	vector<double> param_values_new;

	int nchoice,i,j,l,k;

	//find how many non-zero entries we're after
	nchoice=0;
	for(i=1;i<=nparamFC;i++) if(ichoice[i]==1) nchoice++;
	cout<<nchoice<<endl;

	//store vectors
	param_tags_new.push_back("_fiducial");
	param_values_new.push_back(0.0);

	for(i=1;i<=nparamFC;i++){
		if(ichoice[i]){
			param_tags_new.push_back(param_tags[i]);
			param_values_new.push_back(param_values[i]);
		}
	}
	param_tags=param_tags_new;
	param_values=param_values_new;

	fisherNew=dmatrix(1,nchoice,1,nchoice);

	//copy relevant entries of fisher into fisherNew
	l=0;
	for(i=1;i<=nparamFC;i++){
		if(ichoice[i]) l++;
		k=1;
		for(j=1;j<=nparamFC;j++){
			if(ichoice[i] && ichoice[j]){
				fisherNew[l][k]=fisherFC[i][j];
				k++;
			}
		}
	}

	free_dmatrix(fisherFC,1,nparamFC,1,nparamFC);
	nparamFC=nchoice;
	fisherFC=dmatrix(1,nparamFC,1,nparamFC);

	//save reduced matrix
	for(i=1;i<=nparamFC;i++){
		for(j=1;j<=nparamFC;j++){
			fisherFC[i][j]=fisherNew[i][j];
		}
	}

	// clean up memory
	free_dmatrix(fisherNew,1,nchoice,1,nchoice);

}



///////////////////////////////////////////////////////////////
// Return errors
//////////////////////////////////////////////////////////////
void Fisher::printErrors()
{
	int i;
	invertMatrix(fisherFC,nparamFC,ifisherFC);

	cout<<"Errors"<<endl;
	for(i=1;i<=nparamFC;i++){
		cout<<"\t"<<param_tags[i]<<"\t"<<param_values[i]<<"\t"<<sqrt(ifisherFC[i][i])<<endl;
	}

}

/////////////////////////////////////////////////////////////////////
// copy and combine fisher classes
/////////////////////////////////////////////////////////////////////

void Fisher::copyFisherMatrix(fisherData fisherData1, CosmParam *fiducial)
{
   double **fisher;
	exp_z=fisherData1.exp_z;

	param_tags=fisherData1.param_tags;
	param_values=fisherData1.param_values;

	nparamFC=param_tags.size()-1;

	fisher=dmatrix(1,nparamFC,1,nparamFC);
	loadFisherFile(fisherData1.fisher_file,fisher,nparamFC);
	setFisherMatrix(fisher, nparamFC, fiducial);
        fisher_flag=1;
        free_dmatrix(fisher,1,nparamFC,1,nparamFC);
}


//combine two fisher matrices
void Fisher::combineFisherMatrices(fisherData fisherData1, fisherData fisherData2, CosmParam *fiducial)
{
//	double **fisher1in, **fisher2in;
	vector<double> param_values1, param_values2;
	vector<string> param_tags1, param_tags2;
 	double exp_z1, exp_z2;

	string temp;

	int nparam1, nparam2;
	double **fisher1, **fisher2;
	double **fisherTot;
	int i,j, nparamTot(0), nexcess1, nexcess2;
	int *indx1, *indx2, indxMark;
	vector<string> param_tags_new;
	vector<double> param_values_new;
	vector<string> redshift_tags;
	string rtag;

	exp_z1=fisherData1.exp_z;
	exp_z2=fisherData2.exp_z;

	param_tags1=fisherData1.param_tags;
	param_values1=fisherData1.param_values;

	param_tags2=fisherData2.param_tags;
	param_values2=fisherData2.param_values;

	nparam1=param_tags1.size()-1;
	nparam2=param_tags2.size()-1;

	fisher1=dmatrix(1,nparam1,1,nparam1);
	fisher2=dmatrix(1,nparam2,1,nparam2);
	loadFisherFile(fisherData1.fisher_file,fisher1,nparam1);
	loadFisherFile(fisherData2.fisher_file,fisher2,nparam2);

		 param_tags_new.push_back("_fiducial");
		 param_values_new.push_back(0.0);


	//relabel parameters that can't compare between redshifts
	redshift_tags.push_back("_bias");
	redshift_tags.push_back("_lda");
	redshift_tags.push_back("_lh");
	redshift_tags.push_back("_lg");
	redshift_tags.push_back("_mag");
	redshift_tags.push_back("_deltaz");

	for(j=0;j<=redshift_tags.size()-1;j++){
		rtag=redshift_tags[j];

		for(i=1;i<=nparam1;i++){
			 if(param_tags1[i]==rtag){
				temp=rtag+"_z";
				temp=temp+numToString(exp_z1);
				 param_tags1[i]=temp;
			}
		}
		for(i=1;i<=nparam2;i++){
			 if(param_tags2[i]==rtag){
				temp=rtag+"_z";
				temp=temp+numToString(exp_z2);
				 param_tags2[i]=temp;
			}
		}
	}

	//tote up shared parameters
	for(i=1;i<=nparam1;i++){
		for(j=1;j<=nparam2;j++){
			if(param_tags1[i]==param_tags2[j]){
			 	if(param_values1[i]!=param_values2[j]) cout<<"warning: mismatch in fiducial values: "<<param_tags1[i]<<endl;
				 nparamTot++;
			}
		}
	}
	cout<<nparamTot<<"\t shared parameters"<<endl;
	nexcess1=nparam1-nparamTot;
	nexcess2=nparam2-nparamTot;

	nparamTot+=nexcess1+nexcess2;  //total number of parameters
	fisherTot=dmatrix(1,nparamTot,1,nparamTot);

	//now determine indx for combining matrices
	indx1=ivector(1,nparam1);
	indx2=ivector(1,nparam2);

	for(i=1;i<=nparam1;i++){
		 indx1[i]=i;  //put in matrix 1 in same order
		 param_tags_new.push_back(param_tags1[i]);
		 param_values_new.push_back(param_values1[i]);
	}

	//now set up indx2
	indxMark=nparam1+1;
	for(i=1;i<=nparam2;i++) indx2[i]=-1; //negate index 2
	for(j=1;j<=nparam2;j++){
		for(i=1;i<=nparam1;i++){
			if(param_tags1[i]==param_tags2[j]) indx2[j]=i;
		}
		if(indx2[j]<0){  //if don't find parameter put at end
			indx2[j]=indxMark;
			indxMark++;
		}
	}

	//update parameter tag list
	for(j=1;j<=nparam2;j++){
		if(indx2[j]>nparam1){
		 param_tags_new.push_back(param_tags2[j]);
		 param_values_new.push_back(param_values2[j]);
		}
	}

	//now combine fisher matrices
	//zero fisher matrix
	for(i=1;i<=nparamTot;i++){
		for(j=1;j<=nparamTot;j++){
			fisherTot[i][j]=0.0;
		}
	}

	//add in information from fisher1
	for(i=1;i<=nparam1;i++){
		for(j=1;j<=nparam1;j++){
			fisherTot[indx1[i]][indx1[j]]+=fisher1[i][j];
		}
	}

	//add in information from fisher2
	for(i=1;i<=nparam2;i++){
		for(j=1;j<=nparam2;j++){
			fisherTot[indx2[i]][indx2[j]]+=fisher2[i][j];
		}
	}

	//printMatrix(fisher1,nparam1);
	//printMatrix(fisher2,nparam2);
	//printMatrix(fisherTot,nparamTot);
	//now update class values

	setFisherMatrix(fisherTot, nparamTot, fiducial);
	param_tags.clear();
	param_values.clear();



	for(i=0;i<=nparamTot;i++){
		param_tags.push_back(param_tags_new[i]);
		param_values.push_back(param_values_new[i]);
	}

	free_dmatrix(fisher1,1,nparam1,1,nparam1);
	free_dmatrix(fisher2,1,nparam2,1,nparam2);
	free_dmatrix(fisherTot,1,nparamTot,1,nparamTot);

	saveFisher(".dat");

}

//////////////////////////////////////////////////////////////
// Power spectrum from CAMB transfer function
//////////////////////////////////////////////////////////////

//calculate the CDM power spectrum at the present day
//I assume that the transfer function is constant between
//z=0 (where CAMB calculates it) and z=0.3, z=3.0 (where I apply it).
//
//Choose normalisation of power spectrum to agree with CAMB
double Fisher::powerSpectrum(double k, CosmParam *c, Spline *spline)
{
  double ns,ps(0.0);
  double kfid(0.05);

  if(transfer_spline_flag==1){
	  ns=c->nscal+c->alpha*log(k/kfid)/2.0;
  	  ps=c->ascal2*1.0e-10*pow(k/kfid,ns-1.0);  //primordial power spectrum

	  //note that CAMB transfer function is tabulated [k/h, T(k)]
  	  if(spline!=NULL){
     	 //matter power spectrum
     	ps*=k*2.0*PI*PI;
     	ps*= pow(spline->returnValue(k/c->h),2.0);  //h^{-3}Mpc^{3}
//	ps*=pow(c->h,3.0+ns)
//     	ps*=pow(c->h,-3.0);  //Mpc^{3}
      	return ps;
  	}
  }else{
     	ps=spline->returnValue(k/c->h);  //h^{-3}Mpc^{3}
     	ps*=pow(c->h,-3.0);  //Mpc^{3}
	ps*=c->ascal2/fiducialFC.ascal2;
      	return ps;
  }

  //cout <<"Primordial power spectrum"<<endl;
  return ps;

}

//fiducial power spectrum
double Fisher::powerSpectrumFid(double k, double z)
{
	static int iflag(0);
	static Spline splineFID;
	string tagm, base;
	base=dirbase;

	if(iflag==0){
	iflag++;
	base=dirbase;
	cout<<base<<endl;
	fiducialFC=loadCosmParam(base+"_fiducial"+"_param_m.dat");
	assignSpline(&splineFID,base,"_fiducial_m.dat");
	}

	return powerSpectrum(k,&fiducialFC,&splineFID);
}

/////////////////////////////////////////////////////////////////////////////
// LaTex Table
/////////////////////////////////////////////////////////////////////////////
void Fisher::latexTable()
{
	int i,j,ntag;
	int found_flag(0);
	string file;
	vector<string> tags;
	string line;
	double err;
	ofstream fout;

	//put this here so it gets run by everything
	fisherContour();
	normaliseFisher();

	//parameter list and ordering
	tags.push_back("_ommhh");      //1
	tags.push_back("_ombhh");      //2
	tags.push_back("_oml");        //3
	tags.push_back("_w0");         //6
	tags.push_back("_nscal");      //4
	tags.push_back("_ascal2");     //5
	tags.push_back("_tau");        //7
//	tags.push_back("_ts");         //8
//	tags.push_back("_omk");        //9
	tags.push_back("_yHe");	       //10
	tags.push_back("_alpha");      //11
	tags.push_back("_Mnu");	//12
	tags.push_back("_bias");       //13

//  	tags.push_back("_lda");
//  	tags.push_back("_lh");
//  	tags.push_back("_lg");
  	tags.push_back("_epsilon");
  	tags.push_back("_eta");
 	tags.push_back("_xi");


	file=exp_name+"_latex.txt";
	fout.open(file.c_str());

	//header line

	for(j=0;j<=tags.size()-1;j++) line=line+"& $"+getLatexFromTag(tags[j])+"$ ";

	ntag=tags.size()-1;
	for(i=1;i<=param_tags.size()-1;i++){
		for(j=0;j<=ntag;j++){
			if(param_tags[i]==tags[j]) found_flag=1;
		}
		if(found_flag==0){
			tags.push_back(param_tags[i]);
			line=line+"& "+param_tags[i]+" ";
		}
		found_flag=0;
	}
	line=line+" \\\\";

	fout<<line<<endl;
	line="";

	//param line
	for(j=0;j<=tags.size()-1;j++){
		for(i=1;i<=param_tags.size()-1;i++){
			if(param_tags[i]==tags[j]){
				 line=line+"& "+numToStringLong(param_values[i],4)+" ";
				found_flag=1;
			}
		}
		if(found_flag==0) line=line+"& - ";
		found_flag=0;
	}
	line=line+" \\\\";
	fout<<line<<endl;
	line="";

	//error line
	for(j=0;j<=tags.size()-1;j++){
		for(i=1;i<=param_tags.size()-1;i++){
			if(param_tags[i]==tags[j]){
				err=sqrt(ifisherFC[i][i]);
				 line=line+"& "+numToStringLong(err,3)+" ";
				found_flag=1;
			}
		}
		if(found_flag==0) line=line+"& - ";
		found_flag=0;
	}
	line=line+" \\\\";

	fout<<line<<endl;
	fout.close();

}


// Use string tag to assign parameter value from CosmParam structure
//
string Fisher::getLatexFromTag(string tag)
{
  if(tag=="_ommhh") return "\\Omega_mh^2";
  if(tag=="_ombhh") return "\\Omega_bh^2";
  if(tag=="_oml") return "\\Omega_\\Lambda";
  if(tag=="_nscal") return "n_s";
  if(tag=="_ascal2") return "A^2_s";
  if(tag=="_w0") return "w";
  if(tag=="_tau") return "\\tau";
  if(tag=="_ts") return "T/S";
  if(tag=="_omk") return "\\Omega_k";
  if(tag=="_h") return "h";
  if(tag=="_bias") return "b";
  if(tag=="_omb") return "\\Omega_b";
  if(tag=="_omm") return "\\Omega_m";
  if(tag=="_omnu") return "\\Omega_\\nu";
  if(tag=="_alpha") return "\\alpha";
  if(tag=="_omnuhh") return "\\Omega_\\nu h^2";
  if(tag=="_yHe") return "Y_{He}";
  if(tag=="_Mnu") return "M_\\nu";
  if(tag=="_lda") return "ln D_A";
  if(tag=="_lh") return "ln H";
  if(tag=="_lg") return "ln G";
  if(tag=="_epsilon") return "\\epsilon";    //double use of these parameters
  if(tag=="_eta") return "\\eta";
  if(tag=="_xi") return "\\xi";
  if(tag=="_ntens") return "n_T";
  if(tag=="_thetas") return "\\theta_s";
  if(tag=="_w1") return "w_a";
  if(tag=="_omchh") return "\\Omega_ch^2";
  if(tag=="_omdmhh") return "\\omega_{dm}";
  if(tag=="_fnu") return "f_\\nu";


//	cout<<"Error latex for tag not known: "<<tag<<endl;
	tag.erase(0,1);
	return tag;
}

////////////////////////////////////////////////////////////////////////
// Normalised Fisher matrix
////////////////////////////////////////////////////////////////////////
// calculate normalised Fisher matrix - this gives the
// matrix of cross-correlations
void Fisher::normaliseFisher()
{
	double **fisher_norm;
	int i,j,k;
	double pivot;
	string file, tag;

	fisher_norm=dmatrix(1,nparamFC,1,nparamFC);

	//copy fisher matrix
	for(i=1;i<=nparamFC;i++){
		for(j=1;j<=nparamFC;j++){
			fisher_norm[i][j]=fisherFC[i][j];
		}
	}

	for(k=1;k<=nparamFC;k++){
		pivot=sqrt(fisherFC[k][k]);
		for(i=1;i<=nparamFC;i++){
			fisher_norm[k][i]/=pivot;
			fisher_norm[i][k]/=pivot;
		}
	}


	file=exp_name+"_fisher_norm.txt";
	writeMatrixS(fisher_norm,nparamFC,file);


	free_dmatrix(fisher_norm,1,nparamFC,1,nparamFC);
}


////////////////////////////////////////////////////////////////////////
// Project onto dark energy parameters
////////////////////////////////////////////////////////////////////////

// Project from logH, logD_A onto dark energy parameters
//
// de_flag=1 oml,  2: oml, w0  3: oml, w0, w1 4: oml, w0, w1, omk
//         5 oml, w0, omk 6 oml, omk
//
void Fisher::projectDarkEnergy(int de_flag)
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

//erase references to distances from new parameter list
   for(j=0;j<=param_tags.size()-1;j++){
	for(i=0;i<tags_new.size();i++){
  	 if(tags_new[i].find("_lh",0)!=string::npos){
		 tags_new.erase(tags_new.begin()+i,tags_new.begin()+i+1);
		 values_new.erase(values_new.begin()+i,values_new.begin()+i+1);
	}
  	 if(tags_new[i].find("_lda",0)!=string::npos){
		 tags_new.erase(tags_new.begin()+i,tags_new.begin()+i+1);
		 values_new.erase(values_new.begin()+i,values_new.begin()+i+1);
  	}
	}
   }

	//add in ommhh if not already present
   j=0;
   for(i=0;i<=tags_new.size()-1;i++){
  	if(tags_new[i].find("_ommhh",0)!=string::npos) j++;
   }
   if(j!=1){
	tags_new.push_back("_ommhh");
       	values_new.push_back(fiducialFC.ommhh);
   }

	//add in other dark energy parameters
	tags_new.push_back("_oml");
	values_new.push_back(fiducialFC.oml);
	if(de_flag>1 && de_flag<6){
		 tags_new.push_back("_w0");
		values_new.push_back(fiducialFC.w0);
	}
	if(de_flag>2  && de_flag<5){
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
					fnm=deDerivatives(param_tags[n],tags_new[i]);
					fnm*=deDerivatives(param_tags[m],tags_new[j]);
					fnm*=fisherFC[n][m];
					fisherNew[i][j]+=fnm;
				//	cout<<param_tags[n]<<"\t"<<tags_new[i]<<"\t"<<param_tags[m]<<"\t"<<tags_new[j]<<"\t"<<fisherFC[n][m]<<"\t"<<fnm<<endl;
				}
			}

		}
	}

	//now update class values
	for(i=0;i<=tags_new.size()-1;i++){
//		values_new.push_back(getParamFromTag(tags_new[i],&fiducialFC));
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

double Fisher::deDerivatives(string tagn1, string tagi)
{
   double dhdoml, ddadoml, dhdw0, ddadw0, dhdw1, ddadw1;
   CosmParam cm, cp;
   double pstep;
   double zuse;
   string tagn;

   tagn=tagn1;

   zuse=exp_z;
   if(tagn.find("_lh_z",0)!=string::npos){
	tagn.erase(0,5);
	zuse=strtod(tagn.c_str(),NULL);
	tagn="_lh";
	}

   if(tagn.find("_lda_z",0)!=string::npos){
	tagn.erase(0,6);
	zuse=strtod(tagn.c_str(),NULL);
	tagn="_lda";
	}

   // derivatives w.r.t to dark energy parameters
   if(tagn=="_lh" && tagi=="_oml"){

	cm=fiducialFC;
	cp=fiducialFC;
	pstep=0.01*(1.0-cp.oml);
	cm.oml-=pstep;
	cp.oml+=pstep;
	dhdoml=(hubbleZ(zuse,&cp)-hubbleZ(zuse,&cm))/(2.0*pstep);
	return dhdoml/hubbleZ(zuse,&fiducialFC);
   }else if(tagn=="_lda" && tagi=="_oml"){

	cm=fiducialFC;
	cp=fiducialFC;
	pstep=0.01*(1.0-cp.oml);
	cm.oml-=pstep;
	cp.oml+=pstep;
	ddadoml=(angDiamDistance(zuse,&cp)-angDiamDistance(zuse,&cm))/(2.0*pstep);
	return ddadoml/angDiamDistance(zuse,&fiducialFC);
   }

   if(tagn=="_lh" && tagi=="_w0"){

	cm=fiducialFC;
	cp=fiducialFC;
	pstep=0.01*cp.w0;
	cm.w0-=pstep;
	cp.w0+=pstep;
	dhdw0=(hubbleZ(zuse,&cp)-hubbleZ(zuse,&cm))/(2.0*pstep);
       	return dhdw0/hubbleZ(zuse,&fiducialFC);
  }else if(tagn=="_lda" && tagi=="_w0"){

	cm=fiducialFC;
	cp=fiducialFC;
	pstep=0.01*cp.w0;
	cm.w0-=pstep;
	cp.w0+=pstep;
	ddadw0=(angDiamDistance(zuse,&cp)-angDiamDistance(zuse,&cm))/(2.0*pstep);
       	return ddadw0/angDiamDistance(zuse,&fiducialFC);
  }

   if(tagn=="_lh" && tagi=="_w1"){

	cm=fiducialFC;
	cp=fiducialFC;
	pstep=0.001;
	cm.w1-=pstep;
	cp.w1+=pstep;
	dhdw1=(hubbleZ(zuse,&cp)-hubbleZ(zuse,&cm))/(2.0*pstep);
       	return dhdw1/hubbleZ(zuse,&fiducialFC);
  }else if(tagn=="_lda" && tagi=="_w1"){

	cm=fiducialFC;
	cp=fiducialFC;
	pstep=0.001;
	cm.w1-=pstep;
	cp.w1+=pstep;
	ddadw1=(angDiamDistance(zuse,&cp)-angDiamDistance(zuse,&cm))/(2.0*pstep);
       	return ddadw1/angDiamDistance(zuse,&fiducialFC);
  }

   if(tagn=="_lh" && tagi=="_ommhh"){

	cm=fiducialFC;
	cp=fiducialFC;
	pstep=0.01*cp.ommhh;
	cm.ommhh-=pstep;
	cp.ommhh+=pstep;
	dhdw0=(hubbleZ(zuse,&cp)-hubbleZ(zuse,&cm))/(2.0*pstep);
       	return dhdw0/hubbleZ(zuse,&fiducialFC);
  }else if(tagn=="_lda" && tagi=="_ommhh"){

	cm=fiducialFC;
	cp=fiducialFC;
	pstep=0.01*cp.ommhh;
	cm.ommhh-=pstep;
	cp.ommhh+=pstep;
	ddadw0=(angDiamDistance(zuse,&cp)-angDiamDistance(zuse,&cm))/(2.0*pstep);
       	return ddadw0/angDiamDistance(zuse,&fiducialFC);
  }

   if(tagn=="_lh" && tagi=="_omk"){

	cm=fiducialFC;
	cp=fiducialFC;
	pstep=0.01;
	cm.omk-=pstep;
	cp.omk+=pstep;
	dhdw0=(hubbleZ(zuse,&cp)-hubbleZ(zuse,&cm))/(2.0*pstep);
       	return dhdw0/hubbleZ(zuse,&fiducialFC);
  }else if(tagn=="_lda" && tagi=="_omk"){

	cm=fiducialFC;
	cp=fiducialFC;
	pstep=0.01;
	cm.omk-=pstep;
	cp.omk+=pstep;
	ddadw0=(angDiamDistance(zuse,&cp)-angDiamDistance(zuse,&cm))/(2.0*pstep);
       	return ddadw0/angDiamDistance(zuse,&fiducialFC);
  }

  if(tagn==tagi) return 1.0;
  if(tagn!=tagi) return 0.0;

	cout<<"Error: in deDerivative for tags: "<<tagn<<"\t"<<tagi<<"\t"<<endl;
	return -1.0;
}

////////////////////////////////////////////////////////////////////////
// Project onto dark energy parameters - w_p and w_a
////////////////////////////////////////////////////////////////////////

// Project from logH, logD_A onto dark energy parameters
//
void Fisher::projectDarkEnergyWP()
{
   int nparamNew;
   double **fisherNew;
   int i,j, n, m;
   vector<string> tags_new;
   vector<double> values_new;
   double fnm;
   double wp, ap;

//create new fisher matrix
  tags_new=param_tags;
  values_new=param_values;

   //calculate pivot scale factor and value of w(a_p)
   m=-1;
   n=-1;
   for(i=0;i<=param_tags.size()-1;i++){
      if(param_tags[i].find("_w0",0)!=string::npos) m=i;
      if(param_tags[i].find("_w1",0)!=string::npos) n=i;
   }
   if(m<0 || n<0){
      cout<<"Error: must have w0 and w1 parameters to project onto w_p and w_a"<<endl;
      return;
   }

   ap=1.0+ifisherFC[m][n]/ifisherFC[n][n];
   wp=fiducialFC.w0+(1.0-ap)*fiducialFC.w1;
   cout<<"(ap,wp)="<<ap<<"\t"<<wp<<"\t zp: "<<1.0/ap-1.0<<endl;

//substitute wp for w0
   for(j=0;j<=param_tags.size()-1;j++){
	for(i=0;i<tags_new.size();i++){
  	 if(tags_new[i].find("_w0",0)!=string::npos){
		 tags_new[i]="_wp";
		 values_new[i]=wp;
         }

	}
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
                                   fnm=deDerivativesWP(param_tags[n],tags_new[i],ap);
                                   fnm*=deDerivativesWP(param_tags[m],tags_new[j],ap);
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

//  p_n    -> q_i
// (w0,wa) -> (wp,wa)
//
//  calculate dp_n/dq_i
//
double Fisher::deDerivativesWP(string tagn1, string tagi, double ap1)
{
   string tagn;
   double ap;

   tagn=tagn1;
   ap=ap1;

   // derivatives w.r.t to dark energy parameters
   if(tagn=="_w0" && tagi=="_wp"){
	return 1.0;
   }else if(tagn=="_w1" && tagi=="_wp"){
      return 0.0;
   }

   if(tagn=="_w0" && tagi=="_w1"){
	return -(1.0-ap);
   }else if(tagn=="_w1" && tagi=="_w1"){
      return 1.0;
   }

  if(tagn==tagi) return 1.0;
  if(tagn!=tagi) return 0.0;

	cout<<"Error: in deDerivativeWP for tags: "<<tagn<<"\t"<<tagi<<"\t"<<endl;
	return -1.0;
}



////////////////////////////////////////////////////////////////////////
// Project onto inflationary parameters
////////////////////////////////////////////////////////////////////////

//THIS NEEDS MODIFICATION AS SHOULD BE PROJECTING FROM
//(ns, alpha,r) to (epsilon, eta,xi)  but thinks its doing the otherway round

// Project from logH, logD_A onto dark energy parameters
//
// de_flag=1 oml,  2: oml, w0  3: oml, w0, w1
//
void Fisher::projectInflation(int inf_flag)
{
   int nparamNew;
   double **fisherNew;
   int i,j, n, m;
   vector<string> tags_new;
   vector<double> values_new;
   double fnm;
   Inflation sr;

//create new fisher matrix
  tags_new=param_tags;
  values_new=param_values;

//erase references to distances from new parameter list
   for(j=0;j<=param_tags.size()-1;j++){
	for(i=0;i<tags_new.size();i++){
  	 if(tags_new[i].find("_epsilon",0)!=string::npos){
		 tags_new.erase(tags_new.begin()+i,tags_new.begin()+i+1);
		 values_new.erase(values_new.begin()+i,values_new.begin()+i+1);
	}
  	 if(tags_new[i].find("_eta",0)!=string::npos){
		 tags_new.erase(tags_new.begin()+i,tags_new.begin()+i+1);
		 values_new.erase(values_new.begin()+i,values_new.begin()+i+1);
  	}
  	 if(tags_new[i].find("_xi",0)!=string::npos){
		 tags_new.erase(tags_new.begin()+i,tags_new.begin()+i+1);
		 values_new.erase(values_new.begin()+i,values_new.begin()+i+1);
  	}
	}
   }

	tags_new.push_back("_nscal");
	values_new.push_back(sr.getNscalar(fiducialFC.ts,fiducialFC.nscal,fiducialFC.alpha,0.0));
	if(inf_flag>1){
		 tags_new.push_back("_alpha");
		values_new.push_back(sr.getRunning(fiducialFC.ts,fiducialFC.nscal,fiducialFC.alpha,0.0));
	}
	if(inf_flag>2){
		 tags_new.push_back("_ts");
		values_new.push_back(sr.getTS(fiducialFC.ts,fiducialFC.nscal,fiducialFC.alpha,0.0));
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
					fnm=infDerivatives(param_tags[n],tags_new[i]);
					fnm*=infDerivatives(param_tags[m],tags_new[j]);
					fnm*=fisherFC[n][m];
					fisherNew[i][j]+=fnm;
				//	cout<<param_tags[n]<<"\t"<<tags_new[i]<<"\t"<<param_tags[m]<<"\t"<<tags_new[j]<<"\t"<<fisherFC[n][m]<<"\t"<<fnm<<endl;
				}
			}

		}
	}

	//now update class values
	for(i=0;i<=tags_new.size()-1;i++){
//		values_new.push_back(getParamFromTag(tags_new[i],&fiducialFC));
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

double Fisher::infDerivatives(string tagn1, string tagi)
{
   double dndi;
   CosmParam cm, cp;
   double pstep;
   string tagn;
   Inflation sr;

   tagn=tagn1;

   // derivatives w.r.t to dark energy parameters
   if(tagi=="_nscal" && tagn=="_epsilon"){

	cm=fiducialFC;
	cp=fiducialFC;
	pstep=0.01*fiducialFC.ts;
	if(pstep<1.0e-6) pstep=0.0001;
	cm.ts-=pstep;
	cp.ts+=pstep;
	dndi=(sr.getNscalar(cp.ts,cp.nscal,cp.alpha,0.0)-sr.getNscalar(cm.ts,cm.nscal,cm.alpha,0.0))/(2.0*pstep);
	return dndi;
   }
   if(tagi=="_nscal" && tagn=="_eta"){

	cm=fiducialFC;
	cp=fiducialFC;
	pstep=0.01*fiducialFC.nscal;
	if(pstep<1.0e-6) pstep=0.0001;
	cm.nscal-=pstep;
	cp.nscal+=pstep;
	dndi=(sr.getNscalar(cp.ts,cp.nscal,cp.alpha,0.0)-sr.getNscalar(cm.ts,cm.nscal,cm.alpha,0.0))/(2.0*pstep);
	return dndi;
   }
   if(tagi=="_nscal" && tagn=="_xi"){

	cm=fiducialFC;
	cp=fiducialFC;
	pstep=0.01*fiducialFC.alpha;
	if(pstep<1.0e-6) pstep=0.0001;
	cm.alpha-=pstep;
	cp.alpha+=pstep;
	dndi=(sr.getNscalar(cp.ts,cp.nscal,cp.alpha,0.0)-sr.getNscalar(cm.ts,cm.nscal,cm.alpha,0.0))/(2.0*pstep);
	return dndi;
   }
   if(tagi=="_alpha" && tagn=="_epsilon"){

	cm=fiducialFC;
	cp=fiducialFC;
	pstep=0.01*fiducialFC.ts;
	if(pstep<1.0e-6) pstep=0.0001;
	cm.ts-=pstep;
	cp.ts+=pstep;
	dndi=(sr.getRunning(cp.ts,cp.nscal,cp.alpha,0.0)-sr.getRunning(cm.ts,cm.nscal,cm.alpha,0.0))/(2.0*pstep);
	return dndi;
   }
   if(tagi=="_alpha" && tagn=="_eta"){

	cm=fiducialFC;
	cp=fiducialFC;
	pstep=0.01*fiducialFC.nscal;
	if(pstep<1.0e-6) pstep=0.0001;
	cm.nscal-=pstep;
	cp.nscal+=pstep;
	dndi=(sr.getRunning(cp.ts,cp.nscal,cp.alpha,0.0)-sr.getRunning(cm.ts,cm.nscal,cm.alpha,0.0))/(2.0*pstep);
	return dndi;
   }
   if(tagi=="_alpha" && tagn=="_xi"){

	cm=fiducialFC;
	cp=fiducialFC;
	pstep=0.01*fiducialFC.alpha;
	if(pstep<1.0e-6) pstep=0.0001;
	cm.alpha-=pstep;
	cp.alpha+=pstep;
	dndi=(sr.getRunning(cp.ts,cp.nscal,cp.alpha,0.0)-sr.getRunning(cm.ts,cm.nscal,cm.alpha,0.0))/(2.0*pstep);
	return dndi;
   }
   if(tagi=="_ts" && tagn=="_epsilon"){

	cm=fiducialFC;
	cp=fiducialFC;
	pstep=0.01*fiducialFC.ts;
	if(pstep<1.0e-6) pstep=0.0001;
	cm.ts-=pstep;
	cp.ts+=pstep;
	dndi=(sr.getTS(cp.ts,cp.nscal,cp.alpha,0.0)-sr.getTS(cm.ts,cm.nscal,cm.alpha,0.0))/(2.0*pstep);
	return dndi;
   }
   if(tagi=="_ts" && tagn=="_eta"){

	cm=fiducialFC;
	cp=fiducialFC;
	pstep=0.01*fiducialFC.nscal;
	if(pstep<1.0e-6) pstep=0.0001;
	cm.nscal-=pstep;
	cp.nscal+=pstep;
	dndi=(sr.getTS(cp.ts,cp.nscal,cp.alpha,0.0)-sr.getTS(cm.ts,cm.nscal,cm.alpha,0.0))/(2.0*pstep);
	return dndi;
   }
   if(tagi=="_ts" && tagn=="_xi"){

	cm=fiducialFC;
	cp=fiducialFC;
	pstep=0.01*fiducialFC.alpha;
	if(pstep<1.0e-6) pstep=0.0001;
	cm.alpha-=pstep;
	cp.alpha+=pstep;
	dndi=(sr.getTS(cp.ts,cp.nscal,cp.alpha,0.0)-sr.getTS(cm.ts,cm.nscal,cm.alpha,0.0))/(2.0*pstep);
	return dndi;
   }

  if(tagn==tagi) return 1.0;
  if(tagn!=tagi) return 0.0;

	cout<<"Error: in deDerivative for tags: "<<tagn<<"\t"<<tagi<<"\t"<<endl;
	return -1.0;
}

////////////////////////////////////////////////////////////////////////
// Project onto neutrino parameters
////////////////////////////////////////////////////////////////////////

// Project from m1, m2, m3 onto M_nu
//
// de_flag=0 Mnu  1  m1+m2, m3
//
void Fisher::projectNeutrinoMasses(int mnu_flag)
{
   int nparamNew;
   double **fisherNew;
   int i,j, n, m;
   vector<string> tags_new;
   vector<double> values_new;
   double fnm;
   Inflation sr;

//create new fisher matrix
  tags_new=param_tags;
  values_new=param_values;

//erase references to distances from new parameter list
   for(j=0;j<=param_tags.size()-1;j++){
	for(i=0;i<tags_new.size();i++){
  	 if(tags_new[i].find("_numass1",0)!=string::npos){
		 tags_new.erase(tags_new.begin()+i,tags_new.begin()+i+1);
		 values_new.erase(values_new.begin()+i,values_new.begin()+i+1);
	}
  	 if(tags_new[i].find("_numass2",0)!=string::npos){
		 tags_new.erase(tags_new.begin()+i,tags_new.begin()+i+1);
		 values_new.erase(values_new.begin()+i,values_new.begin()+i+1);
  	}
  	 if(mnu_flag!=1 && tags_new[i].find("_numass3",0)!=string::npos){
		 tags_new.erase(tags_new.begin()+i,tags_new.begin()+i+1);
		 values_new.erase(values_new.begin()+i,values_new.begin()+i+1);
  	}
	}
   }

	if(mnu_flag==1){
		tags_new.push_back("_Msum12");
	values_new.push_back(fiducialFC.numass1+fiducialFC.numass2);
	}else{
		tags_new.push_back("_Mnu");
	values_new.push_back(fiducialFC.numass1+fiducialFC.numass2+fiducialFC.numass3);
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
					fnm=mnuDerivatives(param_tags[n],tags_new[i]);
					fnm*=mnuDerivatives(param_tags[m],tags_new[j]);
					fnm*=fisherFC[n][m];
					fisherNew[i][j]+=fnm;
				//	cout<<param_tags[n]<<"\t"<<tags_new[i]<<"\t"<<param_tags[m]<<"\t"<<tags_new[j]<<"\t"<<fisherFC[n][m]<<"\t"<<fnm<<"\t"<<fisherNew[i][j]<<endl;
				}
			}

		}
	}

	//now update class values
	for(i=0;i<=tags_new.size()-1;i++){
//		values_new.push_back(getParamFromTag(tags_new[i],&fiducialFC));
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

//partial derivatives of m1 w.r.t M_nu given by fraction of Mnu
//that each m1 takes up - i.e. from hierarchy
double Fisher::mnuDerivatives(string tagn1, string tagi)
{
   double dndi;
   CosmParam cm, cp;
   double mnu;
   string tagn;

   tagn=tagn1;

   // derivatives w.r.t to dark energy parameters
   if(tagi=="_Mnu" && tagn=="_numass1"){
	mnu=fiducialFC.numass1+fiducialFC.numass2+fiducialFC.numass3;
	dndi=fiducialFC.numass1/mnu;
	return dndi;
   }
   if(tagi=="_Mnu" && tagn=="_numass2"){
	mnu=fiducialFC.numass1+fiducialFC.numass2+fiducialFC.numass3;
	dndi=fiducialFC.numass2/mnu;
	return dndi;
   }
   if(tagi=="_Mnu" && tagn=="_numass3"){
	mnu=fiducialFC.numass1+fiducialFC.numass2+fiducialFC.numass3;
	dndi=fiducialFC.numass3/mnu;
	return dndi;
   }

   // derivatives w.r.t to dark energy parameters
   if(tagi=="_Msum12" && tagn=="_numass1"){
	mnu=fiducialFC.numass1+fiducialFC.numass2;
	dndi=fiducialFC.numass1/mnu;
	return dndi;
   }
   if(tagi=="_Msum12" && tagn=="_numass2"){
	mnu=fiducialFC.numass1+fiducialFC.numass2;
	dndi=fiducialFC.numass2/mnu;
	return dndi;
   }

  if(tagn==tagi) return 1.0;
  if(tagn!=tagi) return 0.0;

	cout<<"Error: in deDerivative for tags: "<<tagn<<"\t"<<tagi<<"\t"<<endl;
	return -1.0;
}


////////////////////////////////////////////////////////////////////////
// Project from omegam h^2 to omegac h ^2
////////////////////////////////////////////////////////////////////////

// Project onto cold dark matter parameter
//
// Recall: omm=omc+omb+omnu  -> omc=omm-omb-omnu
//
void Fisher::projectOmegaC()
{
   int nparamNew;
   double **fisherNew;
   int i,j, n, m;
   vector<string> tags_new;
   vector<double> values_new;
   double fnm;
   double omchh, omnuhh;

//create new fisher matrix
  tags_new=param_tags;
  values_new=param_values;

  omchh=fiducialFC.ommhh-fiducialFC.ombhh-fiducialFC.omnuhh;
  if(omchh<0.0) cout<<"omchh negative in projectOmegaC()"<<endl;

//substitute omegachh for omegamhh
   for(j=0;j<=param_tags.size()-1;j++){
	for(i=0;i<tags_new.size();i++){
  	 if(tags_new[i].find("_ommhh",0)!=string::npos){
		 tags_new[i]="_omchh";
		 values_new[i]=omchh;
         }

	}
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
                                   fnm=omegacDerivatives(param_tags[n],tags_new[i]);
                                   fnm*=omegacDerivatives(param_tags[m],tags_new[j]);
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

//derivatives for transformation
// (ommhh,ombhh,omnuhh) -> (omchh,ombhh,omnuhh)
//omchh=ommhh-ombhh-omnuhh
//
//ommhh=omchh+ombhh+omnuhh
//ombhh=ombhh
//omnuhh=omnuhh
//
//  p_n    -> q_i
//
//  calculate dp_n/dq_i
double Fisher::omegacDerivatives(string tagn1, string tagi)
{
   string tagn;

   tagn=tagn1;

   // derivatives of p parameters  w.r.t to q parameter set
   if(tagn=="_ommhh" && tagi=="_omchh"){
	return 1.0;
   }else if(tagn=="_ombhh" && tagi=="_omchh"){
      return 0.0;
   }else if(tagn=="_Mnu" && tagi=="_omchh"){
      return 0.0;
   }

   if(tagn=="_ommhh" && tagi=="_ombhh"){
	return 1.0;
   }else if(tagn=="_ombhh" && tagi=="_ombhh"){
      return 1.0;
   }else if(tagn=="_Mnu" && tagi=="_ombhh"){
      return 0.0;
   }

   if(tagn=="_ommhh" && tagi=="_Mnu"){
      return 1.0*omnuhhFromMnu(1.0);  //conversion from omnuhh to Mnu
   }else if(tagn=="_ombhh" && tagi=="_Mnu"){
      return 0.0;
   }else if(tagn=="_Mnu" && tagi=="_Mnu"){
      return 1.0;
   }


  if(tagn==tagi) return 1.0;
  if(tagn!=tagi) return 0.0;

	cout<<"Error: in omegacDerivative for tags: "<<tagn<<"\t"<<tagi<<"\t"<<endl;
	return -1.0;
}

////////////////////////////////////////////////////////////////////////
// Project from omegam h^2 to omegadm h ^2
////////////////////////////////////////////////////////////////////////

// Project onto cold dark matter parameter
//
// Recall: omm=omc+omnu  -> omc=omm-omnu
//
void Fisher::projectOmegaDM()
{
   int nparamNew;
   double **fisherNew;
   int i,j, n, m;
   vector<string> tags_new;
   vector<double> values_new;
   double fnm;
   double omdmhh, omnuhh;

//create new fisher matrix
  tags_new=param_tags;
  values_new=param_values;

  omdmhh=fiducialFC.ommhh-fiducialFC.ombhh;
  if(omdmhh<0.0) cout<<"omdmhh negative in projectOmegaC()"<<endl;

//substitute omegachh for omegamhh
   for(j=0;j<=param_tags.size()-1;j++){
	for(i=0;i<tags_new.size();i++){
  	 if(tags_new[i].find("_ommhh",0)!=string::npos){
		 tags_new[i]="_omdmhh";
		 values_new[i]=omdmhh;
         }

	}
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
                                   fnm=omegadmDerivatives(param_tags[n],tags_new[i]);
                                   fnm*=omegadmDerivatives(param_tags[m],tags_new[j]);
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

//derivatives for transformation
// (ommhh,ombhh,omnuhh) -> (omdmhh,ombhh,omnuhh)
//omchh=ommhh-ombhh-omnuhh
//
//ommhh=omdmhh+ombhh
//ombhh=ombhh
//omnuhh=omnuhh
//
//  p_n    -> q_i
//
//  calculate dp_n/dq_i
double Fisher::omegadmDerivatives(string tagn1, string tagi)
{
   string tagn;

   tagn=tagn1;

   // derivatives of p parameters  w.r.t to q parameter set
   if(tagn=="_ommhh" && tagi=="_omdmhh"){
	return 1.0;
   }else if(tagn=="_ombhh" && tagi=="_omdmhh"){
      return 0.0;
   }else if(tagn=="_Mnu" && tagi=="_omdmhh"){
      return 0.0;
   }

   if(tagn=="_ommhh" && tagi=="_ombhh"){
	return 1.0;
   }else if(tagn=="_ombhh" && tagi=="_ombhh"){
      return 1.0;
   }else if(tagn=="_Mnu" && tagi=="_ombhh"){
      return 0.0;
   }

   if(tagn=="_ommhh" && tagi=="_Mnu"){
      return 0.0;
   }else if(tagn=="_ombhh" && tagi=="_Mnu"){
      return 0.0;
   }else if(tagn=="_Mnu" && tagi=="_Mnu"){
      return 1.0;
   }


  if(tagn==tagi) return 1.0;
  if(tagn!=tagi) return 0.0;

	cout<<"Error: in omegadmDerivative for tags: "<<tagn<<"\t"<<tagi<<"\t"<<endl;
	return -1.0;
}


////////////////////////////////////////////////////////////////////////
// Project from Mnu to f_nu
////////////////////////////////////////////////////////////////////////

// Project onto cold dark matter parameter
//
// Recall: omnuhh=omdmhh*fnu
//
void Fisher::projectFNu()
{
   int nparamNew;
   double **fisherNew;
   int i,j, n, m;
   vector<string> tags_new;
   vector<double> values_new;
   double fnm;
   double fnu;
   int found_flag(0);


//check have cdm parameter
  for(i=0;i<tags_new.size();i++){
     if(tags_new[i].find("_omdmhh",0)!=string::npos) found_flag=1;
  }
  if(found_flag==0){
     cout<<"First projecting onto omegac parameter set"<<endl;
     projectOmegaDM();
  }

//create new fisher matrix
  tags_new=param_tags;
  values_new=param_values;

  fnu=fiducialFC.omnuhh/getParamFromTag("_omchh",&fiducialFC);

//substitute fnu for Mnu
   for(j=0;j<=param_tags.size()-1;j++){
	for(i=0;i<tags_new.size();i++){
  	 if(tags_new[i].find("_Mnu",0)!=string::npos){
		 tags_new[i]="_fnu";
		 values_new[i]=fnu;
         }

	}
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
                                   fnm=fnuDerivatives(param_tags[n],tags_new[i]);
                                   fnm*=fnuDerivatives(param_tags[m],tags_new[j]);
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

//derivatives for transformation
// (omdmhh,Mnu) -> (omdmhh,fnu)
// f_nu=omnuhh/omchh
//  => omnuhh=f_nu * omdmhh
//     omdmhh = omdmhh
//
//  p_n    -> q_i
//
//  calculate dp_n/dq_i
double Fisher::fnuDerivatives(string tagn1, string tagi)
{
   string tagn;
   double omdmhh,fnu;

   tagn=tagn1;
   omdmhh=fiducialFC.ommhh-fiducialFC.ombhh;
   fnu=fiducialFC.omnuhh/omdmhh;

   // derivatives of p parameters  w.r.t to q parameter set
   if(tagn=="_ommhh" && tagi=="_omdmhh"){
	return 1.0;
   }else if(tagn=="_Mnu" && tagi=="_omdmhh"){
      return fnu/omnuhhFromMnu(1.0);
   }

   if(tagn=="_ommhh" && tagi=="_fnu"){
	return 0.0;
   }else if(tagn=="_Mnu" && tagi=="_fnu"){
      return omdmhh/omnuhhFromMnu(1.0);
   }

  if(tagn==tagi) return 1.0;
  if(tagn!=tagi) return 0.0;

	cout<<"Error: in fnuDerivative for tags: "<<tagn<<"\t"<<tagi<<"\t"<<endl;
	return -1.0;
}

///////////////////////////////////////////////////////////////////////////////
// Figure of merit
///////////////////////////////////////////////////////////////////////////////
double Fisher::figureOfMerit(vector<string> params)
{
	int ndim;
	double volume,fom;
	double **M, **iM;
	int i,j, *indx;

	ndim=params.size();
	M=dmatrix(1,ndim,1,ndim);
	iM=dmatrix(1,ndim,1,ndim);
	indx=ivector(1,ndim);

	for(i=1;i<=ndim;i++){
		indx[i]=0;
		for(j=1;j<=nparamFC;j++){
			if(params[i-1]==param_tags[j]) indx[i]=j;
		}
		if(indx[i]==0) cout<<"Error did not find tag in FoM:"<<params[i-1]<<endl;
	}

	//form sub matrix of inverse fisher matrix
	for(i=1;i<=ndim;i++){
		for(j=1;j<=ndim;j++){
			//cout<<i<<"\t"<<j<<"\t"<<param_tags[indx[i]]<<"\t"<<param_tags[indx[j]]<<endl;
			M[i][j]=ifisherFC[indx[i]][indx[j]];
		}
	}

	//invert reduced fisher matrix
	invertMatrix(M,ndim,iM);

	//get volume of ellipse
	volume=1.0/sqrt(detMatrix(iM,ndim));
	fom=1.0/volume;

	free_dmatrix(M,1,ndim,1,ndim);
	free_dmatrix(iM,1,ndim,1,ndim);
	free_ivector(indx,1,ndim);
	return fom;
}

double Fisher::marginalisedErrors(vector<string> params_marg)
{
	int ndim;
	double **M, **iM, **A;
	int i,j, *indx;
        int nparam;
        vector<string> params, tags_new;
        vector<double> values_new;

        tags_new=param_tags;
        values_new=param_values;

        //form list of parameters less those being marginalised

        params_marg.push_back("_fiducial");

        for(j=0;j<=params_marg.size()-1;j++){
           for(i=0;i<tags_new.size();i++){
              if(tags_new[i].find(params_marg[j],0)!=string::npos){
		 tags_new.erase(tags_new.begin()+i,tags_new.begin()+i+1);
		 values_new.erase(values_new.begin()+i,values_new.begin()+i+1);
              }
           }
        }

        params=tags_new;


	ndim=params.size();
	M=dmatrix(1,ndim,1,ndim);
	iM=dmatrix(1,ndim,1,ndim);
	A=dmatrix(1,ndim,1,ndim);
	indx=ivector(1,ndim);

        //locate tags in fisher matrix
	for(i=1;i<=ndim;i++){
		indx[i]=0;
		for(j=1;j<=nparamFC;j++){
			if(params[i-1]==param_tags[j]) indx[i]=j;
		}
		if(indx[i]==0) cout<<"Error did not find tag in marginalisedErrors:"<<params[i-1]<<endl;
	}

	//form sub matrix of inverse fisher matrix
	for(i=1;i<=ndim;i++){
		for(j=1;j<=ndim;j++){
			//cout<<i<<"\t"<<j<<"\t"<<param_tags[indx[i]]<<"\t"<<param_tags[indx[j]]<<endl;
			iM[i][j]=ifisherFC[indx[i]][indx[j]];
		}
	}

	invertMatrix(iM,ndim,M);
	//print errors
	invertMatrix(M,ndim,A);
//	for(i=1;i<=ndim;i++){
//		cout<<param_tags[indx[i]]<<"\t"<<param_values[indx[i]]<<"\t"<<sqrt(A[i][i])<<endl;
//	}

	nparam=tags_new.size();
	setFisherMatrix(M, nparam, &fiducialFC);
        param_tags=tags_new;
        param_values=values_new;

        //vector<int>::iterator it;

        //it=param_tags.begin();
        param_tags.insert(param_tags.begin(),"_fiducial");
        //it=param_tags.begin();
        param_values.insert(param_values.begin(),0.0);

	free_dmatrix(M,1,ndim,1,ndim);
	free_dmatrix(iM,1,ndim,1,ndim);
	free_dmatrix(A,1,ndim,1,ndim);
	free_ivector(indx,1,ndim);
	return 0.0;
}

////////////////////////////////////////////////////////////////////////
// Add systematic errors to distance measurements
////////////////////////////////////////////////////////////////////////
void Fisher::addSystematicErrors(double sigmaD, double sigmaH)
{
	int i,j;
	string tag;
	double **ifishernew, **fishernew;

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

void Fisher::addSystematicErrorFromTag(string param, double sigma)
{
	int i,j;
	string tag;
	double **ifishernew, **fishernew;

	ifishernew=dmatrix(1,nparamFC,1,nparamFC);
	fishernew=dmatrix(1,nparamFC,1,nparamFC);

	//copy fisherFC into ifishernew
	copyMatrix(ifisherFC,nparamFC,ifishernew);

	//add systematic error in quadrature to inverse fisher matrix
	for(i=1;i<=param_tags.size()-1;i++){
		 tag=param_tags[i];
		   if(tag.find(param,0)!=string::npos){
			 ifishernew[i][i]=sigma*sigma+ifisherFC[i][i];
		   }
	}

	//now invert inverse fisher matrix to get new fisher matrix w/sys err
	invertMatrix(ifishernew,nparamFC,fishernew);

	setFisherMatrix(fishernew,nparamFC,&fiducialFC);

	free_dmatrix(fishernew,1,nparamFC,1,nparamFC);
	free_dmatrix(ifishernew,1,nparamFC,1,nparamFC);
}

void Fisher::addPriorFromTag(string param, double sigma)
{
	int i,j;
	string tag;
	double **fishernew;

	fishernew=dmatrix(1,nparamFC,1,nparamFC);

	//copy fisherFC into ifishernew
	copyMatrix(fisherFC,nparamFC,fishernew);

	//add prior to fisher matrix
	for(i=1;i<=param_tags.size()-1;i++){
		 tag=param_tags[i];
		   if(tag.find(param,0)!=string::npos){
			 fishernew[i][i]=1.0/sigma/sigma+fisherFC[i][i];
		   }
	}

	//now copy updated fisher matrix into class

	setFisherMatrix(fishernew,nparamFC,&fiducialFC);

	free_dmatrix(fishernew,1,nparamFC,1,nparamFC);
}


///////////////////////////////////////////////////////////////////
//return inverse of matrix leaving original untouched
//
// Uses gaussj - somewhat unreliable
//
void invertMatrix2(double **fisher,int nparam, double **ifisher)
{
  int i,j,m(1);
  double **b;

  b=dmatrix(1,nparam,1,m);

  for(i=1;i<=nparam;i++){
    for(j=1;j<=nparam;j++){
      ifisher[i][j]=fisher[i][j];
    }
  }

  for(i=1;i<=nparam;i++){
    b[i][1]=1.0;
  }
  //need this routine from NR
  gaussj(ifisher,nparam,b,m);

  free_dmatrix(b,1,nparam,1,m);
}

//save a matrix to file
void writeMatrix(double **M,int nparam, char *file)
{
  ofstream fout;
  int i,j;

  fout.open(file);
  fout.precision(NPRECISION);
  for(i=1;i<=nparam;i++){
    for(j=1;j<=nparam;j++){
      if(j==nparam){
	fout<<M[i][j]<<endl;
      }else{
	fout<<M[i][j]<<"\t";
      }
    }
  }

  fout.close();
}

void writeMatrixS(double **M,int nparam, string file)
{
  ofstream fout;
  int i,j;

  fout.open(file.c_str());
  fout.precision(NPRECISION);
  for(i=1;i<=nparam;i++){
    for(j=1;j<=nparam;j++){
      if(j==nparam){
	fout<<M[i][j]<<endl;
      }else{
	fout<<M[i][j]<<"\t";
      }
    }
  }

  fout.close();
}

double detMatrix(double **M, int nparam)
{
	int i,j, *indx;
	double d;

	if(nparam>3){
		cout<<"fine if non-singular"<<endl;
		ludcmp(M,nparam,indx,&d);
		for(j=1;j<=nparam;j++) d*=M[j][j];
	}else if(nparam==2){
		//2 by 2 matrix
		d=M[1][1]*M[2][2]-M[1][2]*M[2][1];
	}else{
		//3 by 3 matrix
		d=M[1][1]*(M[2][2]*M[3][3]-M[3][2]*M[2][3]);
		d+= -M[1][2]*(M[2][1]*M[3][3]-M[2][3]*M[3][1]);
		d+=M[1][3]*(M[2][1]*M[3][2]-M[2][2]*M[3][1]);
	}

	return d;
}

// Multiply the square N*N matrices A and B to obtain C
//
void multMatrix(double **A, double **B, double **C, int N)
{
	int i,j,k;

	for(i=1;i<=N;i++){
		for(j=1;j<=N;j++){
			C[i][j]=0.0;
			for(k=1;k<=N;k++){
				C[i][j]+=A[i][k]*B[k][j];
			}
		}
	}

}

//return inverse of matrix leaving original untouched
//
// Uses LU decomposition - fairly reliable even for nearly singular matrices
//
void invertMatrix(double **fisher,int nparam, double **ifisher)
{
  int i,j;
  double **a, d, *col;
  int *indx;

  a=dmatrix(1,nparam,1,nparam);
  col=dvector(1,nparam);
  indx=ivector(1,nparam);

  //write fisher into a since inversion destroys original matrix
  for(i=1;i<=nparam;i++){
    for(j=1;j<=nparam;j++){
      a[i][j]=fisher[i][j];
    }
  }

  //perform inversion

  ludcmp(a,nparam,indx,&d);

  for(j=1;j<=nparam;j++){
	for(i=1;i<=nparam;i++) col[i]=0.0;
	col[j]=1.0;
	lubksb(a,nparam,indx,col);
	for(i=1;i<=nparam;i++) ifisher[i][j]=col[i];
  }

  free_dmatrix(a,1,nparam,1,nparam);
  free_dvector(col,1,nparam);
  free_ivector(indx,1,nparam);
}

//return inverse of matrix leaving original untouched
//
// Uses SVD decomposition - note clear that its an improvement on LU decomp
//
void invertMatrixSVD(double **fisher,int nparam, double **ifisher)
{
  int i,j;
  double wmax, wmin, **a, **u, *w, **v, *x, *col;

  a=dmatrix(1,nparam,1,nparam);
  u=dmatrix(1,nparam,1,nparam);
  v=dmatrix(1,nparam,1,nparam);
  x=dvector(1,nparam);
  w=dvector(1,nparam);
  col=dvector(1,nparam);

  //write fisher into a since inversion destroys original matrix
  for(i=1;i<=nparam;i++){
    for(j=1;j<=nparam;j++){
      a[i][j]=fisher[i][j];
    }
  }

  //perform inversion

  svdcmp(a,nparam,nparam,w,v);
  wmax=0.0;
  for(j=1;j<=nparam;j++) if(w[j]>wmax) wmax=w[j];
  //set threshold
  wmin=1.0e-50*wmax;
  for(j=1;j<=nparam;j++){
	if(w[j]<wmin){
		cout<<"zeroing\t"<<j<<"\t"<<w[j]<<endl;
		w[j]=0.0;
	}
  }

  // calculate inverse matrix
  for(j=1;j<=nparam;j++){
	for(i=1;i<=nparam;i++) col[i]=0.0;
	col[j]=1.0;
	svbksb(a,w,v,nparam,nparam,col,x);
	for(i=1;i<=nparam;i++) ifisher[i][j]=x[i];
  }

  free_dmatrix(a,1,nparam,1,nparam);
  free_dmatrix(u,1,nparam,1,nparam);
  free_dmatrix(v,1,nparam,1,nparam);
  free_dvector(x,1,nparam);
  free_dvector(col,1,nparam);
  free_dvector(w,1,nparam);
}


////////////////////////////////////////////////////////////////////////
// Numerical recipes gaussj routine ///////////////////////////////////
///////////////////////////////////////////////////////////////////////

#define NRANSI
#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}

void gaussj(double **a, int n, double **b, int m)
{
	int *indxc,*indxr,*ipiv;
	int i,icol,irow,j,k,l,ll;
	double big,dum,pivinv,temp;

	indxc=ivector(1,n);
	indxr=ivector(1,n);
	ipiv=ivector(1,n);
	for (j=1;j<=n;j++) ipiv[j]=0;
	for (i=1;i<=n;i++) {
		big=0.0;
		for (j=1;j<=n;j++)
			if (ipiv[j] != 1)
				for (k=1;k<=n;k++) {
					if (ipiv[k] == 0) {
						if (fabs(a[j][k]) >= big) {
							big=fabs(a[j][k]);
							irow=j;
							icol=k;
						}
					}
				}
		++(ipiv[icol]);
		if (irow != icol) {
			for (l=1;l<=n;l++) SWAP(a[irow][l],a[icol][l])
			for (l=1;l<=m;l++) SWAP(b[irow][l],b[icol][l])
		}
		indxr[i]=irow;
		indxc[i]=icol;
		if (a[icol][icol] == 0.0) nrerror("gaussj: Singular Matrix");
		pivinv=1.0/a[icol][icol];
		a[icol][icol]=1.0;
		for (l=1;l<=n;l++) a[icol][l] *= pivinv;
		for (l=1;l<=m;l++) b[icol][l] *= pivinv;
		for (ll=1;ll<=n;ll++)
			if (ll != icol) {
				dum=a[ll][icol];
				a[ll][icol]=0.0;
				for (l=1;l<=n;l++) a[ll][l] -= a[icol][l]*dum;
				for (l=1;l<=m;l++) b[ll][l] -= b[icol][l]*dum;
			}
	}
	for (l=n;l>=1;l--) {
		if (indxr[l] != indxc[l])
			for (k=1;k<=n;k++)
				SWAP(a[k][indxr[l]],a[k][indxc[l]]);
	}
	free_ivector(ipiv,1,n);
	free_ivector(indxr,1,n);
	free_ivector(indxc,1,n);
}
#undef SWAP
#undef NRANSI
