/* fisherVPAIR.cc
 * Contains functions for calculating fisher matrix for a Galaxy Cluster experiment
  */
#include "spline.h"
#include "fisherVPAIR.h"

/*********************************************************************
 **************** Constructor/Destructor *****************************
 ********************************************************************/

 FisherVPAIR::FisherVPAIR()
{
	fisher_flag=0;

}

 FisherVPAIR::~FisherVPAIR()
{
  	//cout<<"FisherCMB destructor called"<<endl;
}



void FisherVPAIR::initExperiment(int sur, string label,double SkyCov, double zmax, double intzzz, double rmin,double rmax, double intrrr, int der_flag)
{
  dFlag=der_flag;
  zmaxx=zmax;
  intzz=intzzz;
  survey=sur;
  Rmin= rmin; Rmax=rmax; deltar=intrrr; // in units of Mpc/h

  redshift_list.clear();
  //generate redshift list
  double z=zmax;
  while (z >= 0.0) {
      redshift_list.push_back(z);
      z -= intzzz;
  }
  redshift_list.push_back(0.0);


  CLL_name="./data/Aresults6/CLDer_"+label;
  nameSave="./data/Aresults6/Fisher_"+label;
  COV_name="./data/Aresults6/COV_"+label;
  iCOV_name="./data/Aresults6/COV_"+label+"_inverse";
  SCCC=SkyCov;
}


/////////////////////////////////////////////////////////////////////////////////////////////
/****************************************Cosmology*******************************************/
/////////////////////////////////////////////////////////////////////////////////////////////
double FisherVPAIR::Eofz(double z, Cosmo cosmo)
{
  double E=sqrt(cosmo.wm*pow(1+z,3.) + cosmo.wr*pow(1+z, 4.) + cosmo.wl*pow(1+z, 3*(1+cosmo.w)) + (1-cosmo.wr-cosmo.wl-cosmo.wm));
  return E;
}

double FisherVPAIR::Growth(double z, void *params)
{
  Cosmo cosmo = *(Cosmo *) params;
  Cosmology my_cosmo(cosmo.wm, cosmo.wl, cosmo.wb,cosmo.h, cosmo.sigma8, cosmo.nspec, cosmo.omnu, cosmo.w, cosmo.wa);

  return my_cosmo.growthFac(z);
}

double FisherVPAIR::Dv(double z, Cosmo cosmo)
{
  gsl_function d;
  d.function=&FisherVPAIR::Growth;
  d.params=&cosmo;
  double result, error;

  gsl_deriv_central(&d, z, 1e-4, &result, &error);
  double Hofz=100*cosmo.h*Eofz(z,cosmo);
  return Hofz*result;
}

/////////////////////////////////////////////////////////////////////////////////////////////
/*****************************Functions for theoretical Vpair (Ferreira)*********************/
/////////////////////////////////////////////////////////////////////////////////////////////
double FisherVPAIR::theory_vpair_ferreira(double r, double z, double fr0, double minmass,double w1,double w0,double wm,double wl,double wb,double sigma8, double nscal)
{
  Cosmo cosmo;
  cosmo.wm=wm; cosmo.wl=wl; cosmo.wr=0; cosmo.sigma8=sigma8, cosmo.nspec=nscal; cosmo.omnu=0.; cosmo.wb=wb;
  cosmo.w=w0; cosmo.wa=w1;
  My_params my_params;
  my_params.z=z;
  my_params.cosmo=cosmo;
  my_params.r=r;
  my_params.fr0=fr0;
  my_params.M=minmass;///1e15;
  Cosmology c(cosmo.wm, cosmo.wl, cosmo.wb,cosmo.h, cosmo.sigma8, cosmo.nspec, cosmo.omnu, cosmo.w, cosmo.wa);

  my_params.moment=1;
  double b=effective_b(my_params, minmass);
  double xi=xi_halo(my_params);
//  my_params.moment=1;
  double xi_bar=xi_halo_bar(my_params);
  double dv=abs(Dv(z, cosmo)) ;
  double da=c.growthFac(z);

//  cout << r << "\t" << xi << "\t" << xi_bar << "\t" << dv << "\t" << da << "\t" << b <<endl;

  return -2./3. * (dv/da) * (r*xi_bar*b) /(1.+xi*b*b);

}


double FisherVPAIR::xi_halo_bar(My_params my_params)
{
  Cosmo cosmo=my_params.cosmo;
  Cosmology c(cosmo.wm, cosmo.wl, cosmo.wb,cosmo.h, cosmo.sigma8, cosmo.nspec, cosmo.omnu, cosmo.w, cosmo.wa);
  double z=my_params.z;
  double R=my_params.r;
  double Da=c.growthFac(z);
  double result, error;
  gsl_integration_workspace *aa=gsl_integration_workspace_alloc(1000);
  gsl_function F;
  F.function =&FisherVPAIR::xi_halo_bar_kernelr;
  F.params=&my_params;
  gsl_integration_qags(&F,0,R,0,1e-4,1000, aa, &result, &error);
  gsl_integration_workspace_free(aa);

  return 3.*Da*Da/(2.*PI*PI*R*R*R) * result;
}

double FisherVPAIR::xi_halo_bar_kernelr(double r, void *params)
{
  My_params my_params= *(My_params *) params;
  My_params new_params; new_params.moment=my_params.moment; new_params.M=my_params.M;new_params.cosmo=my_params.cosmo; new_params.z=my_params.z; new_params.r=r;//, new_params.fr0=my_params.fr0, new_params.sfr0=my_params.sfr0;
  double result, error;
  gsl_integration_workspace *aa=gsl_integration_workspace_alloc(1000);
  gsl_function F;
  F.function =&FisherVPAIR::xi_halo_kernel;
  F.params=&new_params;
  gsl_integration_qags(&F,1e-4,10,0,1e-4,1000, aa, &result, &error);
  gsl_integration_workspace_free(aa);

  return result*r;
}

double FisherVPAIR::xi_halo(My_params my_params)
{
  Cosmo cosmo=my_params.cosmo;
  Cosmology c(cosmo.wm, cosmo.wl, cosmo.wb,cosmo.h, cosmo.sigma8, cosmo.nspec, cosmo.omnu, cosmo.w, cosmo.wa);
  double z=my_params.z;
  double Da=c.growthFac(z);
  double result, error;
  gsl_integration_workspace *aa=gsl_integration_workspace_alloc(1000);
  gsl_function F;
  F.function =&xi_halo_kernel;
  F.params=&my_params;
  gsl_integration_qags(&F,1e-5,10,0,1e-4,1000, aa, &result, &error);
  gsl_integration_workspace_free(aa);
  return Da*Da/(2.*PI*PI*my_params.r) * result;
}


double FisherVPAIR::xi_halo_kernel(double k, void *params)
{
  My_params my_params= *(My_params *) params;
  Cosmo cosmo=my_params.cosmo;
  Cosmology c(cosmo.wm, cosmo.wl, cosmo.wb,cosmo.h, cosmo.sigma8, cosmo.nspec, cosmo.omnu, cosmo.w, cosmo.wa);
  double pk=c.powerSpectrum(k);
  return k*pk*sin(k*my_params.r);//*j1;
}


//Input mass in units of Msun/h
double FisherVPAIR::effective_b(My_params my_params, double minmass)
{
  Cosmo cosmo=my_params.cosmo;
  Cosmology c(cosmo.wm, cosmo.wl, cosmo.wb,cosmo.h, cosmo.sigma8, cosmo.nspec, cosmo.omnu, cosmo.w, cosmo.wa);
  double z=my_params.z;//, k=my_params.k,th, R ;

  double Md=minmass, deltaM= pow(10,0.1),NNN=0., NNN2=0.;

  for (int i=0; i<50; i++) {
    double b=c.biasTinker(z,Md);
    double MassF=c.dndlMTinker(z,Md)/Md;
    NNN += MassF*Md*(deltaM-1.);
    NNN2 += MassF*Md*(deltaM-1.0)*pow(b,my_params.moment);
    Md *= deltaM;

  }

  return NNN2/NNN;
}

/////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////covariance////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
void FisherVPAIR::calccov() {

  CosmParam *c;
  c=&fiducialFC;

  Cosmo cosmo;
  cosmo.wm=c->omm; cosmo.wl=c->oml; cosmo.wr=0; cosmo.sigma8=c->sigma8, cosmo.nspec=c->nscal; cosmo.omnu=0.; cosmo.wb=c->omb;cosmo.w=c->w0;cosmo.wa=c->w1;

  for(int i=0;i<=redshift_list.size()-1;i++) {
    cosmic_variance(cosmo,redshift_list[i],1e13);
  }

}



void FisherVPAIR::cosmic_variance(Cosmo cosmo, double z, double Mmin)
{
  double **cov, **icov;
  //initialize the cov matrix
  int nbin= int((Rmax-Rmin)/deltar);

  cov=dmatrix(1,nbin,1,nbin) ;
  icov=dmatrix(1,nbin,1,nbin) ;

  Cosmology c(cosmo.wm, cosmo.wl, cosmo.wb,cosmo.h, cosmo.sigma8, cosmo.nspec, cosmo.omnu, cosmo.w, cosmo.wa);
  double vol=pow(c.angDiamDistance(z)*pow(1+z,3),3)*4.*PI/3.*SCCC/(4*PI);
  double result, error;
  gsl_integration_workspace *aa=gsl_integration_workspace_alloc(1000);
  gsl_function F;

  My_params_cov my_params_cov;
  my_params_cov.z1=z, my_params_cov.z2=z;
  my_params_cov.deltar=deltar;
  my_params_cov.Mmin=Mmin;
  my_params_cov.cosmo=cosmo;

  My_params my_params1, my_params2;
  my_params1.z=z,my_params1.cosmo=cosmo, my_params1.M=Mmin, my_params1.fr0=0.;
  my_params2.z=z, my_params2.cosmo=cosmo, my_params2.M=Mmin, my_params2.fr0=0;

  // start computing the cov matrix .....
  double dv1=abs(Dv(z,cosmo)), D1=c.growthFac(z);
  double dv2=abs(Dv(z,cosmo)), D2=c.growthFac(z);
  for (int i=0; i<nbin; i++) {
     double r1 = Rmin + deltar*i;
     for (int j=i; j<nbin; j++) {
       double r2 = Rmin + deltar*j;
       my_params1.r=r1, my_params2.r=r2; my_params1.moment=1, my_params2.moment=1;
       double b1=effective_b(my_params1,Mmin);
       double b2=effective_b(my_params2,Mmin);
       double fac1=(dv1/D1)/(1.+b1*b1*xi_halo(my_params1));
       double fac2=(dv2/D2)/(1.+b2*b2*xi_halo(my_params2));

        my_params_cov.r1=r1, my_params_cov.r2=r2;
        F.function =&FisherVPAIR::cv_kernel;
        F.params=&my_params_cov;
        gsl_integration_qagiu(&F,0,0,1e-4,1000, aa, &result, &error);

        cov[i][j]=result*4/(PI*PI*vol) * fac1*fac2;
        if (i!=j) {
            cov[j][i]=cov[i][j];
        }

  //      cout << "Computed cov at "<< "\t" << i << "\t"  << j << "\t" << cov[i][j] <<  endl;
     }
  }
  gsl_integration_workspace_free(aa);
  invertMatrix(cov,nbin,icov);

  std::ostringstream ostrs;
  ostrs << z;
  string strr = ostrs.str();
  writeMatrixS(cov,nbin,COV_name+"_z"+strr);
  writeMatrixS(icov,nbin,iCOV_name+"_z"+strr);

  //return cov;
}

double FisherVPAIR::Wx(double x)
{
  return (2*cos(x) + x*sin(x) )/(x*x*x);
}


double FisherVPAIR::cv_kernel(double k, void *params)
{
  FisherVPAIR fisher;
  My_params_cov my_params_cov= *(My_params_cov *) params;
  Cosmo cosmo=my_params_cov.cosmo;
  double deltar=my_params_cov.deltar;
  double r1=my_params_cov.r1, kr1min=k*r1, kr1max=k*(r1+deltar);
  double r2=my_params_cov.r2, kr2min=k*r2, kr2max=k*(r2+deltar);
  double z1=my_params_cov.z1, z2=my_params_cov.z2;
  Cosmology c(cosmo.wm, cosmo.wl, cosmo.wb,cosmo.h, cosmo.sigma8, cosmo.nspec, cosmo.omnu, cosmo.w, cosmo.wa);
  double D1=c.growthFac(z1), D2=c.growthFac(z2);
  double pk=c.powerSpectrum(k);

  double nz1=0.0, nz2=0.0 ;
  double Md=my_params_cov.Mmin, deltaM= pow(10,0.1);
  for (int i=0; i<50; i++) {
    double MassF=c.dndlMTinker(z1,Md)/Md;
    nz1 += MassF* Md*(deltaM-1.);
    MassF=c.dndlMTinker(z2,Md)/Md;
    nz2 += MassF* Md*(deltaM-1.);
    Md *= deltaM;
  }

  double kernel1 = 3.*(r1*r1*r1 *fisher.Wx(kr1min) - pow(r1+deltar,3.)*fisher.Wx(kr1max))/(pow(r1+deltar,3.)-r1*r1*r1);
  double kernel2 = 3.*(r2*r2*r2 *fisher.Wx(kr2min) - pow(r2+deltar,3.)*fisher.Wx(kr2max))/(pow(r2+deltar,3.)-r2*r2*r2);
  //return (1./nz1   )*(1./nz2 ) * kernel1 * kernel2;
  //return (pk*D1*D1   )*(pk*D2*D2 ) * kernel1 * kernel2;
  return (pk*D1*D1 + 1./nz1 )*(pk*D2*D2 + 1./nz2) * kernel1 * kernel2;
}




/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////Derivatives//////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////

void FisherVPAIR::derivatives()
{
  if (dFlag==1)
    {
      double zz=intzz, intz=intzz, mm=1e13;
      double rmin=20, rmax=100, rr=rmin;
      double pstep;//, p_step;
      double cpw0, cmw0,cpwa, cmwa, cmomm,cpomm,cmomb,cpomb,cmoml,cpoml,cmsigma8,cpsigma8, cmnscal,cpnscal;
      nparamFC=7;
      ofstream fout,fout2;
      double *dNN;
      CosmParam *c;
      c=&fiducialFC;
      double h=c->h;

      pstep=0.001;

      dNN=dvector(1,nparamFC);

      cmw0= c->w0- c->w0*pstep/10.0;
      cpw0= c->w0+ c->w0*pstep/10.0;

      cmwa= c->w1- pstep;
      cpwa= c->w1+ pstep;

      cmomm=c->omm*h*h- c->omm*h*h*pstep;
      cpomm=c->omm*h*h+ c->omm*h*h*pstep;

      cmoml=c->oml- c->oml*pstep;
      cpoml=c->oml+ c->oml*pstep;

      cmsigma8=c->sigma8 - c->sigma8*pstep;;
      cpsigma8=c->sigma8 + c->sigma8*pstep;

      cmnscal=c->nscal- c->nscal*pstep;
      cpnscal=c->nscal+ c->nscal*pstep;

      cmomb=c->omb*h*h- c->omb*h*h*pstep;
      cpomb=c->omb*h*h+ c->omb*h*h*pstep;


      while ( zz < zmaxx ) {
        std::ostringstream ostrs;
        ostrs << zz;
        string strr = ostrs.str();
        string fzname=CLL_name+"_z"+strr;
        fout.open(fzname.c_str());
        rr = rmin;
        while (  rr< rmax ) {
            dNN[1]=(theory_vpair_ferreira(rr, zz, 0, 1e13,c->w1,c->w0,cpomm,c->oml,c->omb,c->sigma8, c->nscal)
                  -theory_vpair_ferreira(rr, zz, 0, 1e13,c->w1,c->w0,cmomm,c->oml,c->omb,c->sigma8, c->nscal))/(2*c->omm*pstep);
            dNN[2]=(theory_vpair_ferreira(rr, zz, 0, 1e13,c->w1,c->w0,c->omm,cpoml,c->omb,c->sigma8, c->nscal)
                  -theory_vpair_ferreira(rr, zz, 0, 1e13,c->w1,c->w0,c->omm,cmoml,c->omb,c->sigma8, c->nscal))/(2*c->oml*pstep);
            dNN[3]=(theory_vpair_ferreira(rr, zz, 0, 1e13,c->w1,c->w0,c->omm,c->oml,c->omb,cpsigma8, c->nscal)
                  -theory_vpair_ferreira(rr, zz, 0, 1e13,c->w1,c->w0,c->omm,c->oml,c->omb,cmsigma8, c->nscal))/(2*c->sigma8*pstep);
            dNN[4]=(theory_vpair_ferreira(rr, zz, 0, 1e13,c->w1,c->w0,c->omm,c->oml,cpomb,c->sigma8, c->nscal)
                  -theory_vpair_ferreira(rr, zz, 0, 1e13,c->w1,c->w0,c->omm,c->oml,cmomb,c->sigma8, c->nscal))/(2*c->omb*pstep);
            dNN[5]=(theory_vpair_ferreira(rr, zz, 0, 1e13,c->w1,c->w0,c->omm,c->oml,c->omb,c->sigma8, cpnscal)
                  -theory_vpair_ferreira(rr, zz, 0, 1e13,c->w1,c->w0,c->omm,c->oml,c->omb,c->sigma8, cmnscal))/(2*c->nscal*pstep);
            dNN[6]=(theory_vpair_ferreira(rr, zz, 0, 1e13,c->w1,cpw0,c->omm,c->oml,c->omb,c->sigma8, c->nscal)
                  -theory_vpair_ferreira(rr, zz, 0, 1e13,c->w1,cmw0,c->omm,c->oml,c->omb,c->sigma8, c->nscal))/(2*c->w0*pstep/10.0);
            dNN[7]=(theory_vpair_ferreira(rr, zz, 0, 1e13,cpwa,c->w0,c->omm,c->oml,c->omb,c->sigma8, c->nscal)
                  -theory_vpair_ferreira(rr, zz, 0, 1e13,cmwa,c->w0,c->omm,c->oml,c->omb,c->sigma8, c->nscal))/(2*pstep);

            fout << zz <<"\t"<< rr << "\t" <<  dNN[1] <<"\t"<< dNN[2] <<"\t"<<dNN[3]  <<"\t"<<dNN[4]<<"\t"<<dNN[5]<<"\t"<<dNN[6]<<"\t"<<dNN[7]<< endl;
        rr += deltar;
      }
      fout.close();
      zz += intz;
    }
    free_dvector(dNN,1,nparamFC);
    }
}

void FisherVPAIR::fisherMatrix()
{
  double **fisher, **ifisher;
  ifstream fin, finc;
  double cov, fM, zzz, rrr, zzz2=redshift_list[1];
  double *dnn;
  CosmParam *c;

  c=&fiducialFC;
  nparamFC=7;
  dnn=dvector(1,nparamFC);
  fisher=dmatrix(1,nparamFC,1,nparamFC);
  ifisher=dmatrix(1,nparamFC,1,nparamFC);



  for(int i=1;i<nparamFC+1;i++){
    for(int j=1;j<=i;j++){
      while (zzz2 <= zmaxx) {
        std::ostringstream ostrs;
        ostrs << zzz2;
        string strr = ostrs.str();
        string CLLz=CLL_name+"_z"+strr, COVz=iCOV_name+"_z"+strr ;
        fin.open(CLLz.c_str());
        finc.open(COVz.c_str());
        fM=0.0;
        fin.clear() ;
        fin.seekg(0, ios::beg) ;
        finc.clear() ;
        finc.seekg(0, ios::beg) ;
        fin>> zzz>>rrr>>dnn[1]>>dnn[2]>>dnn[3]>>dnn[4]>>dnn[5]>>dnn[6]>>dnn[7];
        finc>> cov;
        while (!fin.eof()){
          fin>>zzz>>rrr>>dnn[1]>>dnn[2]>>dnn[3]>>dnn[4]>>dnn[5]>>dnn[6]>>dnn[7];
            if (cov !=0) {
              fM += cov*dnn[i]*dnn[j]*intzz*SCCC;
	    //fM+=(1/NNN)*dnn[i]*dnn[j]*intzz*SCCC;    //I multiply by the sky coverage and the z bin here, then I can change the parameter and don't have to calculate derivatives all the times
            }
        }
    }
        fisher[i][j]=fM;
        if(i!=j) fisher[j][i]=fM;
      }
  }
  free_dvector(dnn,1,nparamFC);

  invertMatrix(fisher,nparamFC,ifisher);
  cout<<"Errors:"<<endl;
  cout<<"OmegaM= "<<c->omm<<"\t"<<"OmegaL= "<<c->oml<<"\t"<<"s8= "<<c->sigma8<<"\t"<<"Omegab= "<<c->omb<<"\t"<<"ns= "<<c->nscal<<"\t"<<"Bm= "<<-0.15<<"\t"<<"Sm= "<<0.25<<"\t"<<"alpha= "<<0<<"\t"<<  "beta= "<<0<<"\t"<<"fR= "<<0<<endl;
  for(int i=1;i<nparamFC+1;i++){
    cout<<sqrt(ifisher[i][i])<<"\t\t";
  }
  cout<<endl;
  writeMatrixS(fisher,nparamFC,nameSave);
  writeMatrixS(ifisher,nparamFC,nameSave+"_inverse");

  free_dmatrix(fisher,1,nparamFC,1,nparamFC);
  free_dmatrix(ifisher,1,nparamFC,1,nparamFC);
  fin.close();

}
