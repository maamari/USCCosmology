/* fisherCMB.cc
 * Contains functions for calculating fisher matrix for a Galaxy experiment
  */

#include <math.h>
#include <iostream>
#include <fstream>
#include "fisherGAL.h"
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

 FisherGAL::FisherGAL()
{
	//cout<<"FisherGAL constructor called"<<endl;
	fisher_flag=0;
	//nl_flag=0;
        nparamFC=nparamCORE+nparamGALX;
       
}

 FisherGAL::~FisherGAL()
{
	//cout<<"FisherGAL destructor called"<<endl;

}

//////////////////////////////////////////////////////////////////////
// Initialisation routine
//////////////////////////////////////////////////////////////////////

//initialise Galaxy survey parameters
void FisherGAL::initSurvey(string name1, double z1, double volume1, double density1, double kmin1, double kmax1, double sigma8gal1, double sigmaz1)
{
	exp_name=name1;
	exp_z=z1;
	volume=volume1;
	density=density1;
	kmin=kmin1;	
	kmax=kmax1;
	sigma8gal=sigma8gal1;
	sigmaz=sigmaz1;

	sigmar=0.0;
 	sigpara2=0.0;
	sigperp2=0.0;

//	dirbase="./data/cosmology";
	dirbase=dirbase+"_z"+numToString(exp_z);

	cout<<"survey: "<<exp_name<<" initialised"<<endl;
}

//////////////////////////////////////////////////////////////////////
// Member functions
//////////////////////////////////////////////////////////////////////

string FisherGAL::getSurveyName()
{
	return exp_name;
}

double FisherGAL::getSurveyZ()
{
	return survey_z;
}

double FisherGAL::getKMin()
{
	return kmin;
}

double FisherGAL::getKMax()
{
	return kmax;
}

double FisherGAL::getKNL(double z)
{
   //return c->nonlinearScale(z);
}


double FisherGAL::getKMaxUse()
{
	double kmaxuse,kcut;
	double cutoff(1.0e-10);

	kmaxuse=getKMax();

	//reduce range if photo z render high k ps zero
	if(sigmar>0.0){
		kcut=sqrt(-log(cutoff))/sigmar;
		cout<<kmaxuse<<"\t"<<kcut<<"\t"<<sigmar<<endl;
		kmaxuse=fmin(kmaxuse,kcut);
	}

	return kmaxuse;
}

double FisherGAL::getVolume()
{
	return volume;
}

double FisherGAL::getDensity()
{
	return density;
}

double FisherGAL::getSigperp2()
{
	return sigperp2;
}

double FisherGAL::getSigpara2()
{
	return sigpara2;
}

//////////////////////////////////////////////////////////////////////
// Supplementary galaxy and bias routines
//////////////////////////////////////////////////////////////////////

//calculate the redshift distortion parameter beta
double FisherGAL::getBeta(double z, CosmParam *c)
{
  double gB;

  gB=pow(omegaMZ(z,c),0.6)/c->bias;
 
  return gB;
}

double FisherGAL::getBias(CosmParam *c)
{
	double bias(0.0);

	if(sigma8gal<0.0){
		return fabs(sigma8gal);
	}else{
		return sigma8gal/c->sigma8;
	}

	return bias;
}



////////////////////////////////////////////////////////////////////////////
//  Set galaxy survey parameter files
////////////////////////////////////////////////////////////////////////////
double FisherGAL::setGalaxyFiles()
{
  int i;
  //parameters for calculating derivatives
  for(i=0;i<=nparamFC;i++) setGalaxyPair(i);
  return 1.0;
}


double FisherGAL::setGalaxyPair(int iflag)
{
  double p_step;
  CosmParam c_p;
  CosmParam c_m;
  string base,name,file,tagp,tagm,tag;
  double vary(0.05);
  double bias;
  int i;
  int noncosm_flag(0);  //indicates noncosmological parameter

  base=dirbase;
  tagp="_param_p.dat";
  tagm="_param_m.dat";

  	tag=param_tags[iflag];
	name=tag;

	//for non-cosmological parameters use fiducial
	if(name=="_bias") noncosm_flag=1;
	if(name=="_lda") noncosm_flag=1;
	if(name=="_lh") noncosm_flag=1;
	if(name=="_lg") noncosm_flag=1;

	if(noncosm_flag==1){
		name="_fiducial";
	}

	file=base+name+tagp;
	//cout<<file<<endl;
	c_p=loadCosmParam(file);
	bias=fiducialFC.bias;
	c_p.bias=bias;

	file=base+name+tagm;
	//cout<<file<<endl;
	c_m=loadCosmParam(file);
	c_m.bias=bias;

	//handle distance parameters
	if(distance_flag==0){
		c_m.lda=log(angDiamDistance(exp_z,&c_m));
		c_p.lda=log(angDiamDistance(exp_z,&c_p));
		c_p.lh=log(hubbleZ(exp_z,&c_p));
		c_m.lh=log(hubbleZ(exp_z,&c_m));
		c_m.lg=log(growthZ(exp_z,&c_m));
		c_p.lg=log(growthZ(exp_z,&c_p));
        }else if(distance_flag==1 && tag=="_fiducial"){
		c_m.lda=log(angDiamDistance(exp_z,&c_m));
		c_p.lda=log(angDiamDistance(exp_z,&c_p));
		c_p.lh=log(hubbleZ(exp_z,&c_p));
		c_m.lh=log(hubbleZ(exp_z,&c_m));
		c_m.lg=log(growthZ(exp_z,&c_m));
		c_p.lg=log(growthZ(exp_z,&c_p));

	}else{
		c_m.lda=fiducialFC.lda;
		c_p.lda=fiducialFC.lda;	
		c_m.lh=fiducialFC.lh;
		c_p.lh=fiducialFC.lh;
		c_m.lg=fiducialFC.lg;
		c_p.lg=fiducialFC.lg;
	}
		
	//handle derivatives with respect to non-cosmological parameters
	if(tag=="_fiducial"){
		bias=getBias(&c_p);
		c_p.bias=bias;
		c_m.bias=bias;
		fiducialFC=c_p;
	}else if(tag=="_bias"){
		param_values[iflag]=bias;  //insert fiducial value of bias
		p_step=vary*bias;
		c_p.bias+=p_step;
		c_m.bias-=p_step;
	}else if(tag=="_lda"){
		param_values[iflag]=c_p.lda;  //insert fiducial value of lda
		p_step=vary*c_p.lda;
		c_p.lda+=p_step;
		c_m.lda-=p_step;
	}else if(tag=="_lg"){
		param_values[iflag]=c_p.lg;  //insert fiducial value of lg
		p_step=vary*c_p.lg;
		c_p.lg+=p_step;
		c_m.lg-=p_step;
	}else if(tag=="_lh"){
		param_values[iflag]=c_p.lh;  //insert fiducial value of lh
		p_step=vary*c_p.lh;
		c_p.lh+=p_step;
		c_m.lh-=p_step;
	}

	//handle file related details for non-cosmological parameters
	if(noncosm_flag==1){
		c_p.file=base+tag+"_p.dat";
		c_m.file=base+tag+"_m.dat";

		c_p.name=base+tag+"_param_p.dat";
		c_m.name=base+tag+"_param_m.dat";

		//replicate fiducial transfer and matter spectrum - wasteful
		name="cp "+dirbase+"_fiducial_p.dat"+" "+dirbase+tag+"_p.dat";
		system(name.c_str()); 
		name="cp "+dirbase+"_fiducial_m.dat"+" "+dirbase+tag+"_m.dat";
		system(name.c_str()); 

		base=base+"_transfer";
		name="cp "+base+"_fiducial_p.dat"+" "+base+tag+"_p.dat";
		system(name.c_str()); 
		name="cp "+base+"_fiducial_m.dat"+" "+base+tag+"_m.dat";
		system(name.c_str()); 
	}

  saveCosmParam(&c_m,c_m.name);
  saveCosmParam(&c_p,c_p.name);

  return 1.0;
}

/////////////////////////////////////////////////////////////////////////
// Calculate Galaxy Fisher Matrix
/////////////////////////////////////////////////////////////////////////

//global variables for fisher matrix integration
double lkint;
//double lkminint,lkmaxint;

//Evaluate the Fisher Matrix via the Tegmark approximation
// Generates the fisher matrix in fisher[nmat][nmat] given
// the survey volume (vsurvey), the galaxy number density (density)
// and a fiducial cosmology (c)
void FisherGAL::fisherMatrixGAL()
{
  double tol(1.0e-5);
  double fM;
  double lkminint,lkmaxint;
  string file_dlPi;
  string file_dlPj;

  int i,j;
  ofstream fout;
  string file;
  string file_fid;
  string tag;
  int nmat;
  double **fisher, **ifisher;

  nparamFC=param_tags.size()-1;
  nmat=nparamFC;
//  for(i=0;i<=nparamFC;i++) cout<<param_tags[i]<<endl;
  string file_deriv[nmat+1];

  fisher=dmatrix(1,nmat,1,nmat);
  ifisher=dmatrix(1,nmat,1,nmat);

  cout<<"setting galaxy files"<<endl;
  setGalaxyFiles();
  cout<<"done"<<endl;

  //distance errors
  sigmar=SPEEDOFLIGHT_CGS*sigmaz;
  sigmar/=UNH*hubbleZ(exp_z,&fiducialFC);
  sigmar/=MPC;

  cout<<sigmar<<endl;

  //non-linear growth factors
  double sigma0, sigmapara, sigmaperp, Dz;
  Dz=growthZ(exp_z,&fiducialFC);
  sigma0=(fiducialFC.sigma8/growthZ(exp_z,&fiducialFC)/0.9);
  sigma0*=12.4/fiducialFC.h;
  sigma0*=0.758*growthZ(exp_z,&fiducialFC);

	cout<<sigma0<<endl;

  sigmaperp=sigma0;
  sigperp2=sigmaperp*sigmaperp;
  sigmapara=sigma0*(1.0+pow(omegaMZ(exp_z,&fiducialFC),0.6));
  sigpara2=sigmapara*sigmapara;

  file_fid=param_tags[0];

  for(i=0;i<=nparamFC-1;i++) file_deriv[i]=param_tags[i+1];

  cout<<"calculaing fisher matrix"<<endl;
  //everything is done using files to create the fisher matrix
  //I'll exploit this to automate the proceedure

    lkminint=log(getKMin());
    lkmaxint=log(getKMaxUse());

  for(i=0;i<nmat;i++){
    for(j=0;j<=i;j++){
      file_dlPi=file_deriv[i];
      file_dlPj=file_deriv[j];
       setFisherKernel(0.0,file_fid,this,file_dlPi,file_dlPj,0);
       fM=qromb(getFisherKernelK,lkminint,lkmaxint,tol);
      setFisherKernel(0.0,file_fid,this,file_dlPi,file_dlPj,2);
      fisher[i+1][j+1]=fM;
      if(i!=j) fisher[j+1][i+1]=fM;
      cout <<file_deriv[i]<<"\t"<<file_deriv[j]<<"\t"<<fM<<endl;
    }
  }

  //save information about fisher matrix and inverse into class
  setFisherMatrix(fisher,nmat,&fiducialFC);

  //write result to file
  tag="_z"+numToString(exp_z)+".dat";
  saveFisher(tag);

  //output error information
  cout <<"Galaxy survey errors: "<<exp_name<<endl;

  for(i=1;i<=nparamFC;i++){
	cout<<"delta"<<param_tags[i]<<"\t"<<param_values[i]<<"\t"<<sqrt(ifisherFC[i][i])<<"\t"<<sqrt(ifisherFC[i][i])/param_values[i]<<endl;
  }

  free_dmatrix(fisher,1,nmat,1,nmat);
  free_dmatrix(ifisher,1,nmat,1,nmat);

}

double getFisherKernelK(double lk)
{
  double gFKM;
  double tol(1.0e-5);

  lkint=lk;
  gFKM=qromb(getFisherKernel,-1.0,1.0,tol);

  return gFKM;
}

double getFisherKernel(double mu)
{
  string blank;
  return setFisherKernel(mu,blank,NULL,blank,blank,1);
}

double setFisherKernel(double mu, string file_fid, FisherGAL *survey1,string file_dlPi1in, string file_dlPj1in,int iflag)
{
  static CosmParam c;
  static FisherGAL *survey;
  double sFK;
  double k,lk;
  static CosmParam ci_m, ci_p;
  static CosmParam cj_m, cj_p;
  double dlPi,dlPj;
  string file,base,tagm,tagp;
  static string file_dlPi1, file_dlPj1;

  tagm="_param_m.dat";
  tagp="_param_p.dat";


  lk=lkint;
  k=exp(lk);

  if(iflag==0){
    survey=survey1;
    base=survey->returnDirbase();

    file_dlPi1=file_dlPi1in;
    file_dlPj1=file_dlPj1in;


    file=base+file_fid+tagm;
    c=survey->loadCosmParam(file);
    file=base+file_dlPi1+tagm;
    ci_m=survey->loadCosmParam(file);
    file=base+file_dlPi1+tagp;
    ci_p=survey->loadCosmParam(file);
    file=base+file_dlPj1+tagm;
    cj_m=survey->loadCosmParam(file);
    file=base+file_dlPj1+tagp;
    cj_p=survey->loadCosmParam(file);

	ci_m.daratio=survey1->daratio(&ci_m);
	ci_p.daratio=survey1->daratio(&ci_p);
	ci_m.hratio=survey1->hratio(&ci_m);
	ci_p.hratio=survey1->hratio(&ci_p);
	cj_m.daratio=survey1->daratio(&cj_m);
	cj_p.daratio=survey1->daratio(&cj_p);
	cj_m.hratio=survey1->hratio(&cj_m);
	cj_p.hratio=survey1->hratio(&cj_p);

    survey->setDerivSplines(file_dlPi1,file_dlPj1,1);
    return 0.0;
  }

  if(iflag==2){
    //clean up spline routines
	    survey->setDerivSplines(file_dlPi1,file_dlPj1,2);
    return 0.0;
  }

  //calculate power spectrum derivatives
  dlPi=survey->derivePowerGAL(k,mu,ci_p,ci_m,file_dlPi1,1); 
  dlPj=survey->derivePowerGAL(k,mu,cj_p,cj_m,file_dlPj1,2);

  //if difference is very small the above method becomes numerically 
  //unstable
  if((dlPi!=dlPi)||(dlPj!=dlPj)) cout<<k<<"\t"<<mu<<"\t"<<"Error in setFisherKernel"<<file_dlPi1<<"\t"<<file_dlPj1<<"\t"<<dlPi<<"\t"<<dlPj<<"\t"<<survey->returnDirbase()<<endl;

  sFK=dlPi*dlPj*survey->volumeEffective(k,mu);
  sFK*=k*k/(4.0*PI*PI);
  sFK/=2.0;   //not sure where this comes from
  sFK*=k;  //integral over log(k)
  
  if(nl_flag==1){
  sFK*=exp(-0.5*k*k*mu*mu*survey->getSigpara2());   //non-linear growth
  sFK*=exp(-0.5*k*k*(1.0-mu*mu)*survey->getSigperp2()); //non-linear growth
  }

  return sFK;

}

//Calculate effective volume of survey for fiducial cosmology
double FisherGAL::volumeEffective(double k, double mu)
{
  double nsur,vsur,ps,veff;

  ps=getPowerSpectrumObs(k,mu,&fiducialFC,&splineFid);
  nsur=density;
  vsur=volume;

  veff=nsur*ps/(nsur*ps+1.0);
  veff=vsur*veff*veff;

  return veff;
}


///////////////////////////////////////////////////////////////////////
//  Code to calculate the galaxy power spectrum
///////////////////////////////////////////////////////////////////////

//calculate the derivative of power sepctrum
double FisherGAL::derivePowerGAL(double k, double mu, CosmParam c_p, CosmParam c_m, string tag, int iflag)
{
	double ps;
	double dlP, p_step, dP;
	Spline *splineP, *splineM;

	if(iflag==1){
		splineP=&splinePow1p;
		splineM=&splinePow1m;
	}else if(iflag==2){
		splineP=&splinePow2p;
		splineM=&splinePow2m;
	}	

    p_step=getParamFromTag(tag,&c_p)-getParamFromTag(tag,&c_m);

//     dlP=(log(getPowerSpectrumObs(k,mu,&c_p,splineP))-log(getPowerSpectrumObs(k,mu,&c_m,splineM)))/(p_step);


  //possibly get better cancellation from direct differencing
//	dP=dlP;
  dlP=((getPowerSpectrumObs(k,mu,&c_p,splineP))-(getPowerSpectrumObs(k,mu,&c_m,splineM)))/(p_step);

  ps=getPowerSpectrumObs(k,mu,&fiducialFC,&splineFid);
  dlP/=ps;

   if(dlP!=dlP) return 0.0;

	return dlP;
} 

//set splines for power spectrum derivatives
void FisherGAL::setDerivSplines(string tag1, string tag2, int iflag)
{
 	string base,tagm,tagp;
	string transfertag;
	transfertag="_transfer";
	tagm="_m.dat";
  	tagp="_p.dat";

	//set splines
	if(iflag==1){
	    base=dirbase;
	assignSpline(&splineFid,base,"_fiducial"+tagp);
	assignSpline(&splinePow1p,base,tag1+tagp);
	assignSpline(&splinePow1m,base,tag1+tagm);
	assignSpline(&splinePow2p,base,tag2+tagp);
	assignSpline(&splinePow2m,base,tag2+tagm);

	return;
	}

	//clean splines
	if(iflag==2){
	splineFid.cleanSpline();	
	splinePow1p.cleanSpline();
	splinePow1m.cleanSpline();
	splinePow2p.cleanSpline();
	splinePow2m.cleanSpline();

	return;
	}

	cout<<"Must enter 1 or 2 for iflag"<<endl;
	return;

}


//Calculate the galaxy power spectrum for the parameters specified in c
//at the observed values of (kref,muref) using a spline transfer function 
double FisherGAL::getPowerSpectrumObs(double kref, double muref,CosmParam *c, Spline *spline)
{
  double mu;
  double ps;
  double beta;
  double pshot;
  double alcockp;
  double kpara2_ref,kperp2_ref;
  double kpara2,kperp2;
  double k,kref2;
  double hratio, daratio;
  int ap_flag(1);
  int rss_flag(1);

  pshot=0.0;
  beta=getBeta(exp_z,c);

  //relate reference k,mu to the true k,mu
  if(ap_flag==1){
	  daratio=exp(c->lda)/exp(fiducialFC.lda);
  	hratio=exp(c->lh)/exp(fiducialFC.lh);
  	kref2=kref*kref;
  	kpara2_ref=kref2*muref*muref;
  	kperp2_ref=kref2*(1.0-muref*muref);
 	 kperp2=kperp2_ref/pow(daratio,2.0);
  	kpara2=kpara2_ref*pow(hratio,2.0);
  	k=sqrt(kperp2+kpara2);
  	alcockp=hratio/pow(daratio,2.0);  //Alcock-Paczynski factor
  }else{
 	 k=kref;
 	mu=muref;
 	 kpara2=k*k*mu*mu;
 	 kperp2=k*k*(1.0-mu*mu);
  }

  if(k!=k)  cout<<k<<"\t"<<kref<<"\t"<<muref<<"\t"<<daratio<<"\t"<<hratio<<"\t"<<alcockp<<endl;

  ps=powerSpectrum(k,c,spline);		      // raw power spectrum
  ps*=pow(c->bias,2.0); 		      //bias
  if(rss_flag==1) ps*=pow(1.0+beta*kpara2/k/k,2.0);//linear redshift distortion
  if(ap_flag==1) ps*=alcockp;                              //AP-effect
  ps*=exp(-kpara2*sigmar*sigmar); 	    //redshift errors

  return ps+pshot;
}

///////////////////////////////////////////////////////////////////////
// Miscellaneous extra quantities that are interesting to output
///////////////////////////////////////////////////////////////////////

double FisherGAL::powerSpectrumError(double k)
{
	double tol(1.0e-5);
	double mumin(-1.0), mumax(1.0);
	double errInt;
	double z;
	double bias;
	string base,file;

	z=exp_z;

    file=dirbase+"_fiducial_param_p.dat";;
    fiducialFC=loadCosmParam(file);
	bias=fiducialFC.bias;

	fiducialFC.lda=log(angDiamDistance(exp_z,&fiducialFC));
	fiducialFC.lh=log(hubbleZ(exp_z,&fiducialFC));

	    base=dirbase;  // use transfer function
	assignSpline(&splineFid,base,"_fiducial_p.dat");

	setKErrorKernelGAL(k,0.0,this,1);
	errInt=qromb(getKErrorKernelGAL,mumin,mumax,tol);
//	errInt/=pow(bias,4.0);  //account for bias that was missed in integrand
	return 1.0/sqrt(fabs(errInt));
}

double setKErrorKernelGAL(double k1, double mu, FisherGAL *survey1, int iflag)
{
	static double k;
	static FisherGAL *survey;
	double sFK;
	double ps;
	double eta(0.5); //k step size dk=eta*k

	if(iflag==1){
		k=k1;
		survey=survey1;
		return 0.0;
	}

  sFK=survey->volumeEffective(k,mu);
  sFK*=k*k*k*eta;
  sFK/=8.0*PI*PI;
  //need a factor of 1/ps/ps here to get power spectrum error
//  ps=survey->powerSpectrumFid(k);
//  sFK/=ps*ps;

	return sFK;
}

double getKErrorKernelGAL(double theta)
{
	return setKErrorKernelGAL(0.0,theta,NULL,0);
}

void FisherGAL::outPowerErrors()
{
	double bias;
	double k, ps, err, veff;
	string file;
	ofstream fout;

	file=exp_name+"_k_errors.dat";

	bias=fiducialFC.bias;

	fout.open(file.c_str());
	k=kmin;
	while(k<2.0){
		ps=powerSpectrumFid(k,exp_z)*bias*bias;
		err=powerSpectrumError(k);   //fractional error
		veff=averagedEffectiveVolume(k);
		fout<<k<<"\t"<<ps<<"\t"<<err<<"\t"<<bias<<"\t"<<veff<<endl;
		k*=1.2;
	}
	fout.close();
	
}


double FisherGAL::averagedEffectiveVolume(double k)
{
	double tol(1.0e-5);
	double mumin(-1.0), mumax(1.0);
	double errInt;
	double z;
	double bias;
	string base,file;

	z=exp_z;

    file=dirbase+"_fiducial_param_p.dat";;
    fiducialFC=loadCosmParam(file);
	bias=fiducialFC.bias;

	fiducialFC.lda=log(angDiamDistance(exp_z,&fiducialFC));
	fiducialFC.lh=log(hubbleZ(exp_z,&fiducialFC));

	    base=dirbase;  // use transfer function
	assignSpline(&splineFid,base,"_fiducial_p.dat");

	setKVolumeKernelGAL(k,0.0,this,1);
	errInt=qromb(getKVolumeKernelGAL,mumin,mumax,tol);
	return errInt;
}

double setKVolumeKernelGAL(double k1, double mu, FisherGAL *survey1, int iflag)
{
	static double k;
	static FisherGAL *survey;
	double sFK;
	double ps;
	double eta(0.5); //k step size dk=eta*k

	if(iflag==1){
		k=k1;
		survey=survey1;
		return 0.0;
	}

  sFK=survey->volumeEffective(k,mu);
  sFK/=2.0;  //averaging over mu
  return sFK;
}

double getKVolumeKernelGAL(double theta)
{
	return setKVolumeKernelGAL(0.0,theta,NULL,0);
}


///////////////////////////////////////////////////////////////////////
//  Code to calculate the galaxy power spectrum derivatives for output
///////////////////////////////////////////////////////////////////////
void FisherGAL::derivGAL(string tag,string file)
{
	CosmParam cp, cm, c;
	ofstream fout;
	double dkmin, dkmax;
 	string base, tagm,tagp;
	double p_step;
	int i,nstep;
	double k,kstep;
	double Pp, Pm, dP, dlP, Ps, dlP2;
	double dlP3, dlP4;
	double mu(1.0);

	Spline splineP, splineM, splineFID;

	tagm="_m.dat";
  	tagp="_p.dat";

	base=dirbase;

	cout<<base<<endl;
	cout<<base+tag+"_param"+tagp<<endl;

	cp=loadCosmParam(base+tag+"_param"+tagp);
	cm=loadCosmParam(base+tag+"_param"+tagm);
	c=loadCosmParam(base+"_fiducial"+"_param"+tagm);

        p_step=getParamFromTag(tag,&cp)-getParamFromTag(tag,&cm);
    
	assignSpline(&splineFID,base,"_fiducial"+tagm);
	assignSpline(&splineP,base,tag+tagp);
	assignSpline(&splineM,base,tag+tagm);

	fiducialFC=c;

	c.lda=log(angDiamDistance(exp_z,&c));
	c.lh=log(hubbleZ(exp_z,&c));

	fiducialFC=c;

	cp.lda=log(angDiamDistance(exp_z,&cp));
	cp.lh=log(hubbleZ(exp_z,&cp));
	cm.lda=log(angDiamDistance(exp_z,&cm));
	cm.lh=log(hubbleZ(exp_z,&cm));

	cout<<c.lh<<"\t"<<c.lda<<endl;
	cout<<cm.lh<<"\t"<<cm.lda<<endl;
	cout<<cp.lh<<"\t"<<cp.lda<<endl;	

        setDerivSplines(tag,tag,1);

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
		dlP2=derivePowerGAL(k,1.0,cp,cm,tag,1);
		dlP3=derivePowerGAL(k,0.5,cp,cm,tag,1);
		dlP4=derivePowerGAL(k,0.0,cp,cm,tag,1);
		
		fout<<k<<"\t"<<Pp<<"\t"<<Pm<<"\t"<<dP<<"\t"<<dlP<<"\t"<<Ps<<"\t"<<dlP2<<"\t"<<dlP3<<"\t"<<dlP4<<endl;
//		cout<<k<<"\t"<<Pp<<"\t"<<Pm<<"\t"<<dP<<"\t"<<dlP<<endl;
	}
	fout.close();
	
	cout<<c.lh<<"\t"<<c.lda<<endl;
	cout<<cm.lh<<"\t"<<cm.lda<<endl;
	cout<<cp.lh<<"\t"<<cp.lda<<endl;
	cout<<base<<endl;

        setDerivSplines(tag,tag,2);

}

