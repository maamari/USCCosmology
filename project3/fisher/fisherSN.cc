/* fisherSN.cc
 * Contains functions for calculating fisher matrix for a Galaxy experiment
  */

#include <math.h>
#include <iostream>
#include <fstream>
#include "fisherSN.h"
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

 FisherSN::FisherSN()
{
	//cout<<"FisherSN constructor called"<<endl;
	fisher_flag=0;
       
}

 FisherSN::~FisherSN()
{
	//cout<<"FisherSN destructor called"<<endl;

}

//////////////////////////////////////////////////////////////////////
// Initialisation routine
//////////////////////////////////////////////////////////////////////

//initialise Galaxy survey parameters
void FisherSN::initSurvey(string name1, double nbin1, double zmin1, double zmax1, double sigmaz1, double sigmaD1)
{
	exp_name=name1;
        nbin=nbin1;
	zmin=zmin1;	
	zmax=zmax1;
	exp_z=(zmin+zmax)/2.0;
	sigmaz=sigmaz1;
 	sigmaD=sigmaD1;

	dirbase=dirbase+"_z"+numToString(exp_z);


	cout<<"SN 1a survey: "<<exp_name<<" initialised"<<endl;
}

//////////////////////////////////////////////////////////////////////
// Member functions
//////////////////////////////////////////////////////////////////////

string FisherSN::getSurveyName()
{
	return exp_name;
}

double FisherSN::getSurveyZ()
{
	return survey_z;
}

double FisherSN::getZMin()
{
	return zmin;
}

double FisherSN::getZMax()
{
	return zmax;
}

double FisherSN::getNbin()
{
	return nbin;
}

void FisherSN::setCosmology(CosmParam fiducial)
{
	fiducialFC=fiducial;
}

////////////////////////////////////////////////////////////////////////
// magnitude of supernova
////////////////////////////////////////////////////////////////////////
double FisherSN::magnitudeSNfull(double z, CosmParam *c, double magoff, double magL, double magQ, double zbias)
{
	double mag;	
	double zuse;

	zuse=z+zbias;

	mag=magnitudeSN(zuse,c);

	mag+=magoff;
	mag+=magL*zuse;
  	mag+=magQ*zuse*zuse;

	return mag;
}

double FisherSN::magnitudeSN(double z, CosmParam *c)
{
	double mag, dL;	

	dL=lumDistance(z,c)*SPEEDOFLIGHT_CGS/H0_CMSMPC;
	mag=5.0*log10(dL)+25.0;

	return mag;
}

////////////////////////////////////////////////////////////////////////
// Calculate Fisher matrix
////////////////////////////////////////////////////////////////////////

void FisherSN::fisherMatrixSN()
{
  double **fisher, **ifisher;
  vector<double> element;
  int i,j,nmat;
  string tag;
  double sigmaMag;
  double dmudmuoff, pstep;
  double dmudmuQ, dmudmuL;
  double fij;
  double sigmaL, sigmaQ;


  //establish parameters
  param_tags.clear();
  param_values.clear();
  param_tags.push_back("_fiducial");
  param_values.push_back(0.0);

  nparamFC=0;
  addParamWithValue("_mag",magnitudeSN(exp_z,&fiducialFC));
  addParamWithValue("_muoff",0.0);
  addParamWithValue("_muL",0.0);
  addParamWithValue("_muQ",0.0);

  nmat=nparamFC;
  fisher=dmatrix(1,nmat,1,nmat);
  ifisher=dmatrix(1,nmat,1,nmat);

  saveCosmParam(&fiducialFC,fiducialFC.name);

  //calculate fisher matrix information
  sigmaMag=sqrt((sigmaD*sigmaD+sigmaz*sigmaz)/nbin);

  pstep=0.01;
  dmudmuoff=(magnitudeSNfull(exp_z,&fiducialFC,pstep,0.0,0.0,0.0)-magnitudeSNfull(exp_z,&fiducialFC,-pstep,0.0,0.0,0.0))/(2.0*pstep);
  dmudmuL=(magnitudeSNfull(exp_z,&fiducialFC,0.0,pstep,0.0,0.0)-magnitudeSNfull(exp_z,&fiducialFC,0.0,-pstep,0.0,0.0))/(2.0*pstep);
  dmudmuQ=(magnitudeSNfull(exp_z,&fiducialFC,0.0,0.0,pstep,0.0)-magnitudeSNfull(exp_z,&fiducialFC,0.0,0.0,-pstep,0.0))/(2.0*pstep);

  element.push_back(0.0);
  element.push_back(1.0);
  element.push_back(dmudmuoff);
  element.push_back(dmudmuL);
  element.push_back(dmudmuQ);
      
  for(i=1;i<=nparamFC;i++){
	for(j=1;j<=nparamFC;j++){
		fij=0.0;
		fisher[i][j]=element[i]*element[j];
	}
  }

  //multiply by overall error factor
  for(i=1;i<=nparamFC;i++){
	for(j=1;j<=nparamFC;j++){
		fisher[i][j]/=pow(sigmaMag,2.0);
	}
  }

  //priors on quadratic mu offset
  sigmaL=0.01/sqrt(2.0);
  sigmaQ=0.01/sqrt(2.0);
  fisher[3][3]+=1.0/pow(sigmaL,2.0);
  fisher[4][4]+=1.0/pow(sigmaQ,2.0);

  //save information about fisher matrix and inverse into class
  setFisherMatrix(fisher,nmat,&fiducialFC);

  //write result to file
  tag="_z"+numToString(exp_z)+".dat";
  saveFisher(tag);
  free_dmatrix(fisher,1,nmat,1,nmat);
  free_dmatrix(ifisher,1,nmat,1,nmat);

  return;
}


////////////////////////////////////////////////////////////////////////
// Add redshift bin
////////////////////////////////////////////////////////////////////////
void FisherSN::addRedshiftBin(double zmin1, double zmax1, double nbin1, double sigmaz1)
{
	fisherData oldData, newData;

	//save old fisherData somewhere safe
	oldData=returnFisherData();
	saveFisher("_outs");
	oldData.fisher_file=fisher_fileFC;

	//update redshift information
	zmin=zmin1;	
	zmax=zmax1;
	nbin=nbin1;
	sigmaz=sigmaz1;
	exp_z=(zmin+zmax)/2.0;

	//calculate new redshift bin data
	fisherMatrixSN();
	newData=returnFisherData();
	
	combineFisherMatrices(oldData,newData,&fiducialFC);

}

////////////////////////////////////////////////////////////////////////
// Project onto dark energy parameters
////////////////////////////////////////////////////////////////////////

// Project from logH, logD_A onto dark energy parameters
//
// de_flag=1 oml,  2: oml, w0  3: oml, w0, w1 
//
void FisherSN::projectDarkEnergy(int de_flag)
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
  	 if(tags_new[i].find("_mag",0)!=string::npos){
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
	if(de_flag>1){
		 tags_new.push_back("_w0");
		values_new.push_back(fiducialFC.w0);
	}
	if(de_flag>2){
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

double FisherSN::deDerivatives(string tagn1, string tagi)
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

   if(tagn.find("_mag_z",0)!=string::npos){
	tagn.erase(0,6);
	zuse=strtod(tagn.c_str(),NULL);
	tagn="_mag";
	}

 if(tagn=="_mag" && tagi=="_ommhh"){

	cm=fiducialFC;
	cp=fiducialFC;
	pstep=0.01*cp.ommhh;
	cm.ommhh-=pstep;
	cp.ommhh+=pstep;
	ddadw0=(magnitudeSN(zuse,&cp)-magnitudeSN(zuse,&cm))/(2.0*pstep);
       	return ddadw0;
  }
 if(tagn=="_mag" && tagi=="_oml"){

	cm=fiducialFC;
	cp=fiducialFC;
	pstep=0.01*cp.oml;
	cm.oml-=pstep;
	cp.oml+=pstep;
	ddadw0=(magnitudeSN(zuse,&cp)-magnitudeSN(zuse,&cm))/(2.0*pstep);
       	return ddadw0;
  }
 if(tagn=="_mag" && tagi=="_w0"){

	cm=fiducialFC;
	cp=fiducialFC;
	pstep=0.01*cp.w0;
	cm.w0-=pstep;
	cp.w0+=pstep;
	ddadw0=(magnitudeSN(zuse,&cp)-magnitudeSN(zuse,&cm))/(2.0*pstep);
       	return ddadw0;
  }
 if(tagn=="_mag" && tagi=="_w1"){

	cm=fiducialFC;
	cp=fiducialFC;
	pstep=0.01;
	cm.w1-=pstep;
	cp.w1+=pstep;
	ddadw0=(magnitudeSN(zuse,&cp)-magnitudeSN(zuse,&cm))/(2.0*pstep);
       	return ddadw0;
  }
 if(tagn=="_mag" && tagi=="_omk"){

	cm=fiducialFC;
	cp=fiducialFC;
	pstep=0.01;
	cm.omk-=pstep;
	cp.omk+=pstep;
	ddadw0=(magnitudeSN(zuse,&cp)-magnitudeSN(zuse,&cm))/(2.0*pstep);
       	return ddadw0;
  }

  if(tagn==tagi) return 1.0;
  if(tagn!=tagi) return 0.0;

	cout<<"Error: in deDerivative for tags: "<<tagn<<"\t"<<tagi<<"\t"<<endl;
	return -1.0;
}
