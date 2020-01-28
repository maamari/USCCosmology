/* fisherGAL_CLUSTER.cc
 * Contains functions for calculating fisher matrix for a Galaxy Cluster experiment
  */
#include "spline.h"
#include "fisherGAL_CLUSTER.h"

/*********************************************************************
 **************** Constructor/Destructor *****************************
 ********************************************************************/

 FisherGAL_CLUSTER::FisherGAL_CLUSTER()
{
	//cout<<"FisherCMB constructor called"<<endl;
	//nchannel=0;
	fisher_flag=0;
        //nparamFC=nparamCORE+nparamCMBX;
        //lminPol=2;

}

 FisherGAL_CLUSTER::~FisherGAL_CLUSTER()
{
  	//cout<<"FisherCMB destructor called"<<endl;
}



void FisherGAL_CLUSTER::initExperiment(int sur, string label, string Mdynname, string FRname,string FRbias, string FRPk ,double SkyCov, double zmax, double intzzz, double FRR, int der_flag)
{
  std::ostringstream ostrs;
  ostrs << FRR;
  string strr = ostrs.str();

  dFlag=der_flag;
  zmaxx=zmax;
  intzz=intzzz;
  survey=sur;

  CLL_name="./data2/Aresults/CLDer_"+label+"_"+strr;
  CLL_namePS="./data2/Aresults/CLDerPS_Tinker_"+label+"_"+strr;
  nameSave="./data2/Aresults/Fisher_"+label+"_"+strr;
  nameSavePS="./data2/Aresults/FisherPS_Tinker_"+label+"_"+strr;
  FR_name=FRname;
  FR_bias=FRbias;
  FR_Pk=FRPk;
  Mdyn_name=Mdynname;
  nameSaveFR="./data/Afr";
  SCCC=SkyCov;
  fr_0=FRR;
}


double FisherGAL_CLUSTER::Omegamz(double z, double wm)
{
  return wm * pow(1+z,3);
}


double FisherGAL_CLUSTER::OmegaLambdaz(double z, double wl,double w, double wa)
{
  double lamZ;

  if (wa !=0) {
      lamZ = wl*pow(1+z,3.0*(1.0+w+wa)) * exp(-3.0*wa*z/(1+z));
  } else {
    lamZ= wl *pow(1+z,3*(1+w));
  }

  return lamZ;
}

double FisherGAL_CLUSTER::wzKernel(double z, void * params)
{
  Cosmo cosmo = *(Cosmo *) params;

  return (1+ ( cosmo.w + cosmo.wa*z/(1+z)) ) / (1+z);
}

double FisherGAL_CLUSTER::comovingdistance(double z, Cosmo cosmo) // h^-1 Mpc
{
  double Dh=cl*1e-5; // h^-1 Mpc
  double result, error;
  gsl_integration_workspace *aa=gsl_integration_workspace_alloc(1000);
  gsl_function F;
  F.function =&FisherGAL_CLUSTER::k;
  F.params=&cosmo;
  gsl_integration_qags(&F,0,z,0,1e-7,1000, aa, &result, &error);
  gsl_integration_workspace_free(aa);
  return Dh*result/cosmo.h;
}

double FisherGAL_CLUSTER::ComovingVolume(double z, double ommhh, double wr, double wl,double w, double wa)
{
  double radian2deg = 57.296;
  double h=sqrt(ommhh/(1-wl));
  double wm=ommhh/(h*h);
  Cosmo cosmo;
  cosmo.wm=wm;
  cosmo.wl=wl;
  cosmo.wr=wr;
  cosmo.w=w;
  cosmo.wa=wa;
  cosmo.h=h;

  /*double mat = Omegamz(z,wm);
  double rad = wr*pow(1+z,4);
  double lam = OmegaLambdaz(z,wl,w, wa);
  double Hofz = 1000 * (100*h) * sqrt(mat+rad+lam); // in units of ms-1 Mpc-1
  CosmParam *c = &fiducialFC;
  c->omm = wm; c->oml = wl; c->w0=w; c->w1=wa; c->h=h;*/
  double Hofz = h*Eofz(wm,wl,wr,w,wa,z); // in units of ms-1 Mpc-1
  double da=angdiamdist(ommhh, wl, wr,w,wa, z);
  double coord=comovingdistance(z,cosmo);

  return coord*coord*cl*1e-5 / (Hofz); // * radian2deg*radian2deg) ;
 // return pow( (1+z)*da,2 ) * cl*1e-5 / (Hofz*h); // * radian2deg*radian2deg) ;
}

double FisherGAL_CLUSTER::getVolume(double z, double ommhh, double wl, double wr, double w, double wa)
{
  double result, error;
  gsl_integration_workspace * ww= gsl_integration_workspace_alloc(1000);
  double h=sqrt(ommhh/(1-wl));
  My_params my_params;
  my_params.wm=ommhh; // wm *h*h
  my_params.wl=wl;
  my_params.wr=wr;
  my_params.w=w;
  my_params.wa=wa;

  gsl_function F;
  F.function =&FisherGAL_CLUSTER::setVol;
  F.params =&my_params;

  gsl_integration_qags (&F, 0, z, 0, 1e-7, 1000, ww, &result, &error);

  gsl_integration_workspace_free(ww);
  return result;

}

double FisherGAL_CLUSTER::setVol (double z, void * params)
{
  My_params my_params = *(My_params *) params;
  FisherGAL_CLUSTER fisher;

  return fisher.ComovingVolume(z,my_params.wm,my_params.wl,my_params.wr,my_params.w,my_params.wa);
}



double FisherGAL_CLUSTER::Eofz(double wm, double wl, double wr, double w, double wa, double z)
{
  //return sqrt(wm*pow(1+z,3.) + wr*pow(1+z,4.) + wl + (1-wl-wm-wr)*pow(1+z,2.));
  return sqrt(Omegamz(z,wm) + OmegaLambdaz(z,wl,w, wa)+ wr*pow(1+z,4));
}

double FisherGAL_CLUSTER::angdiamdist (double ommhh, double wl, double wr, double w, double wa, double z)
{
  gsl_integration_workspace * ww= gsl_integration_workspace_alloc(1000);

  double result, error;
  double h=sqrt(ommhh/(1-wl));
  double wm=ommhh/(h*h);

  Cosmo cosmo;
  cosmo.wm=wm;
  cosmo.wl=wl;
  cosmo.wr=wr;
  cosmo.w=w;
  cosmo.wa=wa;

  gsl_function F;
  F.function =&FisherGAL_CLUSTER::k;
  F.params =&cosmo;

  gsl_integration_qags (&F, 0, z, 0, 1e-7, 1000, ww, &result, &error);
  gsl_integration_workspace_free(ww);
  double Dh = cl*1e-5;
  return Dh*result/(h*(1+z));
}

double FisherGAL_CLUSTER::k (double z, void * params)
{
  Cosmo cosmo = *(Cosmo *) params;
  FisherGAL_CLUSTER fisher;

  return 1.0/fisher.Eofz(cosmo.wm,cosmo.wl,cosmo.wr,cosmo.w,cosmo.wa, z);
}

// search for scale factor for dndm and halo bias files
double FisherGAL_CLUSTER::findNnn_z (double z, double m,string filename, double wm)
{
  ifstream fin;
  double m0,m1, zz, ncon0, nloo0, ncon1, nloo1, mint=pow(10,0.05), intzzz=0.005;
  double m500crit=M200toM500(m,500./wm, z);
  fin.open(filename.c_str());
  if (fin.is_open()) {
    for (int i=0; i< floor(z/intzzz)*70; i++) fin >> m0 >> zz >> ncon0 >> nloo0; // jump to the corresponding z
    for (int i=0; i<70; i++) {
          fin >> m0 >> zz >> ncon0 >> nloo0;
          if ( (m0*mint > m500crit) && ( i!=69) ) {
              fin >> m1 >> zz >> ncon1 >> nloo1;
              return ncon0 + (ncon1-ncon0)/(log10(m1/m0)) * (log10(m500crit/m0));
              //return nloo0 + (nloo1-nloo0)/(log10(m1/m0)) * (log10(m/m0));
          }
    }
  } else {
    cerr << "Error. Cannot open dndm/bias file" << endl;
    exit(0);
  }

  fin.close();
  //return nloo1; //the mass range is beyond the tabulated values,scale according to the last pair
  return ncon1; //the mass range is beyond the tabulated values,scale according to the last pair
}

double FisherGAL_CLUSTER::findNnn (double z, double m,string filename, double wm)
{
  double intzzz=0.005;

  if ( z/intzzz-floor(z/intzzz) == 0 ) {
    return findNnn_z(z,m,filename, wm);
  } else {
      double z0= floor(z/intzzz)*intzzz;
      double z1= z0+intzzz;
      double y_z0= findNnn_z(z0,m,filename,wm);
      double y_z1= findNnn_z(z1,m,filename, wm);
      return y_z0 + (z-z1)/(z0-z1) * (y_z1- y_z0);
  }
}

double FisherGAL_CLUSTER::findlambdaC(double fr0)
{
  ifstream fin;
  fin.open("./data/lambdaC.dat");
  double f, lambc, f1, lambc1;

  if (fin.is_open()) {
    while ( !fin.eof()) {
    fin >> f >> lambc;
    if ( f == fr0 ) {
      return lambc;
    } else if ( f > fr0) {
      fin >> f1 >> lambc1;
      return lambc + (lambc1-lambc)/(f1-f) * (fr0-f);
    }
    }
  }
  else {
    cout << "File not found !" << endl;
  }
  fin.close();
  cout << "Not found !" << endl;
  return 0;
}


double FisherGAL_CLUSTER::findMyn_z (double z, double m, double wm)
{
  ifstream fin;
  double intzzz=0.005;
  double mtrue_0,mtrue_1, msz_0,msz_1, zz, mint=pow(10,0.05);
  double m500crit=M200toM500(m,500./wm, z);
  fin.open(Mdyn_name.c_str());
  cout << Mdyn_name.c_str() << endl;
  
  if (fin.is_open()) {
      for (int i=0; i< floor(z/intzzz)*50; i++) fin >> mtrue_0 >> zz >> msz_0; // jump to the corresponding z
      for (int i=0; i<50; i++) {
          fin >> mtrue_0 >> zz >> msz_0;
          if ( (msz_0*mint > m500crit) && (i != 49) ) {
              fin >> mtrue_1 >> zz >> msz_1;
              double resultm= pow(10,log10(mtrue_0)+ (log10(mtrue_1/mtrue_0))/(log10(msz_1/msz_0)) * (log10(m500crit/msz_0)));
              //if (resultm <= 5e13) return 5e13;
              //else {
                //cout << z << "\t" << m << "\t" << resultm << "\t" << m/resultm << endl;
                return resultm/m500crit *m;
              //}
          }
      }
  } else {
    cerr << "Error. Cannot open Mdyn file" << endl;
    exit(0);
  }

  fin.close();
  //cout << z << "\t" << m << "\t" << m*(mtrue_0/msz_0) << "\t" << m/(m*(mtrue_0/msz_0)) << endl;
  return m*(mtrue_0/msz_0); //the mass range is beyond the tabulated values,scale according to the last pair
}


double FisherGAL_CLUSTER::findMyn (double z, double m, double wm)
{
  double intzzz=0.005;

  if ( z/intzzz-floor(z/intzzz) == 0 ) return findMyn_z(z,m, wm);
  else {
      double z0= floor(z/intzzz)*intzzz;
      double z1= z0+intzz;
      double y_z0= findMyn_z(z0,m, wm);
      double y_z1= findMyn_z(z1,m, wm);
      return y_z0 + (z-z1)/(z0-z1) * (y_z1- y_z0);
  }
}

// Return the Mlim(z) in M200,mean
double FisherGAL_CLUSTER::cal_mlim (double ommhh, double wl, double wr, double w,double wa,double z)
{
  CosmParam *c = &fiducialFC;
  double h=sqrt(ommhh/(1-wl));
  double wm=ommhh/(h*h);
  double A=4.21, B=1.31, C=1.6, sn=3.6;
  c->omm = wm; c->oml = wl; c->w0=w; c->w1=wa; c->h=0.72;

  double da=angDiamDistance(z,c)*cl*1e-5, E=Eofz(wm,wl,wr,w,wa, z), E6=Eofz(wm, wl,wr,w,wa,0.6);
  //double da=angdiamdist(wm,wl,wr,w,wa, z), E=Eofz(wm,wl,wr,w,wa, z);
  double mlim=0, Slim_spt=5.0, // mJy
         Ylim_act=1.0e-3, Ylim_planck=2e-3, // arcmin^2
         masscut=0;

  if (survey == 0 ) { // Planck
        double M200crit;
        /*masscut=3.5e14;
	if (z<=0.11)
	  //mlim=((1e15)/h)*exp(-1.924+8.333*z);
          M200crit=((1e15)/h)*pow(10,-1.924+8.333*z);
        else
          //mlim=((1e15)/h)*exp(-1.2+1.469*atan(pow(z-0.1,0.44)));
          M200crit=((1e15)/h)*pow(10,-1.2+1.469*atan(pow(z-0.1,0.44)));*/
        masscut=1e14*(0.72/0.7);
        if ( z ==0 ) mlim=masscut;
        else M200crit=(1e15*0.72/0.7)*pow(da*da*pow(E,-2./3.)*(Ylim_planck*pow(d2r/60,2)/2.504e-4), 1./1.876);
        if (M200crit <= masscut) M200crit=masscut;
        //mlim=M200crit;
        mlim=MhtoM200(M200crit, 200./wm, z);
  } else if (survey ==1 ) { // SPT
        masscut=0.9e14/h;
        if ( z ==0 ) mlim=masscut;
        else mlim=(5.e14/h)*pow((sqrt(sn*sn-3.)/A)*pow((1+z)/1.6,-C),1.0/B);
        if (mlim <= masscut) mlim=masscut;
        //else mlim=(1e15)*pow(da*da*pow(E,-2./3.)*(Slim_spt/2.592e8), 1./1.876);
  } else if (survey ==2 ) { // ACTpol
        masscut=5.e14/h;
        //masscut=5.e14/h;
        mlim=masscut;
        //if ( z ==0 ) mlim=masscut;
        //else mlim=(1e15)*pow(da*da*pow(E,-2./3.)*(Ylim_act*pow(d2r/60,2)/2.504e-4), 1./1.876);
  } else if (survey ==3) {//SPT pol
        masscut=0.9e14/h;
        if ( z ==0 ) mlim=masscut;
        else mlim=(5.e14/h)*pow((sqrt(sn*sn-3.)/A)*pow((1+z)/1.6,-C),1.0/B)*(3.01/5.95);
        if (mlim <= masscut) mlim=masscut;
  } else if (survey ==4) {//CCAT
        double m500crit;
        masscut=0.9e13/h;
        sn=5.;
        if ( z ==0 ) mlim=masscut;
        else m500crit=0.54*(3.e14/h)*pow((sqrt(sn*sn-3.)/10.2)*pow(E6/E,0.73),1.0/1.33);
        if (m500crit <= masscut) m500crit=masscut;
        mlim=MhtoM200(m500crit, 500./wm, z);
  } else if (survey ==5) {//CCAT3g
        double m500crit;
        masscut=0.9e13/h;
        sn=8.74;
        if ( z ==0 ) mlim=masscut;
        else m500crit=0.54*(3.e14/h)*pow((sqrt(sn*sn-3.)/10.2)*pow(E6/E,0.73),1.0/1.33);
        if (m500crit <= masscut) m500crit=masscut;
        mlim=MhtoM200(m500crit, 500./wm, z);
  } else if (survey >=6) {//CCATwide
        double m500crit;
        masscut=0.9e13/h;
        if (survey ==6) sn=10.9; //33,000 sq deg
        else if (survey ==7) sn=8.83; //10,000 sq deg
        else if (survey ==8) sn=8.27; //5,000 sq deg
        else if (survey ==9) sn=7.83; //1,000 sq deg
        if ( z ==0 ) mlim=masscut;
        else m500crit=0.54*(3.e14/h)*pow((sqrt(sn*sn-3.)/10.2)*pow(E6/E,0.73),1.0/1.33);
        if (m500crit <= masscut) m500crit=masscut;
        mlim=MhtoM200(m500crit, 500./wm, z);
  }

  //double mmin=M200toM500(mlim);
  double mmin=mlim;
  return mmin;
}

double FisherGAL_CLUSTER::fff(double x)
{
  return pow(x,3)*(log(1.0+1.0/x)-1.0/(1.0+x));
}

double FisherGAL_CLUSTER::M200toM500(double M200, double overdensity, double z)
{
  valarray<double> a(4);
  a[0]=0.5116;
  a[1]=-0.4283;
  a[2]=-3.13e-3;
  a[3]=-3.52e-5;

  double c=9.0/(1.0+z) *pow(M200/0.8e13/0.7,-0.13) ;
  //double c=5.0;
  double f_h=overdensity/200. *fff(1.0/c);
  double p=a[1] + a[2]*log(f_h)+a[3]*log(f_h)*log(f_h);
  double xoff=1.0/sqrt(a[0]*pow(f_h,2.0*p)+9.0/16.0) +2.0*f_h;
  return M200*overdensity/200. / pow(c*xoff,3.0);
}

double FisherGAL_CLUSTER::MhtoM200(double Mh, double overdensity, double z)
{
  valarray<double> a(4);
  a[0]=0.5116;
  a[1]=-0.4283;
  a[2]=-3.13e-3;
  a[3]=-3.52e-5;

  double c=9.0/(1.0+z) *pow(Mh/0.8e13/0.7,-0.13) ;
  //double c=9.0/(1.0+z) *pow(M200/0.8e13/0.7,-0.13) ;
  //double c=5.0;
  double f_h=overdensity/200. *fff(1.0/c);
  double p=a[1] + a[2]*log(f_h)+a[3]*log(f_h)*log(f_h);
  double xoff=1.0/sqrt(a[0]*pow(f_h,2.0*p)+9.0/16.0) +2.0*f_h;
  //return M200*overdensity/200. / pow(c*xoff,3.0);
  return Mh/(overdensity/200.) * pow(c*xoff,3.0);
}
double FisherGAL_CLUSTER::SetdndM(double m, void * params)
{
  //FisherGAL_CLUSTER my_params =  *(FisherGAL_CLUSTER *) params;
  //CosmParam c;
  //c=my_params.fiducialFC;
  My_params my_params =  *(My_params *) params;
  //double omnu = 0.0;
  //Cosmology my_cosmo(c.omm, c.oml, c.omb, c.h, c.sigma8, c.nscal, c.omnu, c.w0, c.w1);
  Cosmology my_cosmo(my_params.wm, my_params.wl, my_params.wb, my_params.h, my_params.sigma8, my_params.ns, 0., my_params.w, my_params.wa);
  double z =my_params.z;
  double nn=1.;

  //FisherGAL_CLUSTER fisher;
  //int *choiceVec; choiceVec=ivector(1,22);
  //for (int i=1; i<=22; i++) choiceVec[i]=1;
  //fisher.specifyParams(choiceVec,my_params.c);
  //fisher.initExperiment(0,"",my_params.mdyfile, my_params.filename, "", "", 100 , 1 , 0.05 , 1e-4, 1);
  //fisher.initExperiment(survey,"", "data/Mdyn-Planck-"+FR_name+".dat","data/dndm-Planck-"+FR_name+".dat", "data/haloBias-Planck-"+FR_name+".dat", "data/Pk-Planck-"+FR_name+".dat", SCCC , 1 , 0.05 , fr_0, 1);
  /*if (c.frflag) nn=my_params.findNnn(z,m,my_params.FR_name);
  if (nn<1) nn=1;

  double Bm=c.Bm0*pow((1+z),c.alpha_cluster),
        sigmaMM=c.sigma0*pow((1+z),c.beta),
        xm1=(log(c.m1)-Bm-log(m))/sqrt(2*sigmaMM),
        xm2=(log(c.m2)-Bm-log(m))/sqrt(2*sigmaMM);*/
  return my_cosmo.dndlMTinker(z,m)/m;
  //double MassF=nn*fisher.dndlMTinker(z,m,my_params)/m;
  //return MassF*(erfc(xm2)-erfc(xm1));
}

/*double getNofzFR(double m)
{
  return setNofzFR(m,NULL,0);
}

double setNofzFR(double m, My_params *c1, int iflag)
{
  static My_params *c;
  if (iflag==1) {
    c=c1;
    return 0.0;
  }

  double omnu = 0.0;
  double nn=1;
  //Cosmology my_cosmo(c->wm+c->wr, c->wl, c->wb, h, c->sigma8, c->ns, omnu);
  Cosmology my_cosmo(c->wm+c->wr, c->wl, c->wb, h, c->sigma8, c->ns, omnu, c->w);
  double z = c->z;
  FisherGAL_CLUSTER fisher;
  if (c->frflag) nn=fisher.findNnn(z,m,c->filename);
  if (nn<1) nn=1;

  double Bm=c->Bm0*pow((1+z),c->alpha),
         sigmaMM=c->sigma0*pow((1+z),c->beta),
         xm1=(log(c->m1)-Bm-log(m))/sqrt(2*sigmaMM),
         xm2=(log(c->m2)-Bm-log(m))/sqrt(2*sigmaMM);

  double MassF=nn*fisher.dndlMTinker(z,m,*c)/m;
  //double MassF=nn*my_cosmo.dndlMTinker(z,m)/m;

  return MassF*(erfc(xm2)-erfc(xm1));
}*/

/*   Useless, it's just NoFZ * nn*/
double FisherGAL_CLUSTER::NofzFR(bool frflag, double wa,double w,double wmmmhh, double wlll, double wbbbhh, double wrrr,  double sigmaa8, double nspecc, double zzz, double Mmm1,double Mmm2, double Bm0, double sigma0,double alpha, double beta)
//double FisherGAL_CLUSTER::NofzFR(My_params my_params, double zzz)
{
  double MassF=0, deltaM=pow(10,0.05),xm1=0, xm2=0, Bm=0.1, sigmaMM=0.1, Md=0;
  double NNN=0,nn=1;
  double omnu = 0.0;
  double h=sqrt(wmmmhh/(1-wlll));
  Cosmology my_cosmo(wmmmhh/(h*h), wlll, wbbbhh/(h*h), h, sigmaa8, nspecc, omnu, w, wa);

  Bm=Bm0*pow((1+zzz),alpha);
  sigmaMM=sigma0*pow((1+zzz),beta);
  Md=1e13;
  for(int k=0; k< 100; k++){
    if (frflag) nn=findNnn(zzz,Md,FR_name, wmmmhh/(h*h));
    if (nn <1) nn=1;
    MassF=nn*my_cosmo.dndlMTinker(zzz,Md)/Md;
    //I have to define Bm and sigmaMM
    xm1=(log(Mmm1)-Bm-log(Md))/(sqrt(2)*sigmaMM);
    xm2=(log(Mmm2)-Bm-log(Md))/(sqrt(2)*sigmaMM);
    NNN +=MassF*Md *(deltaM-1) *(erfc(xm2)-erfc(xm1));
    Md *=deltaM;
  }
  //return (NNN*getDVolumeDOmega(zzz,&myfiducial))/2.0;
  //return result*my_cosmo.volumeComoving(zzz);
   return NNN*ComovingVolume(zzz,wmmmhh,wrrr,wlll,w, wa)/2.0;
  //return ComovingVolume(zzz,wmmm,wrrr,wlll,w)/2 * result;

}


double FisherGAL_CLUSTER::SetN(double z, void *params)
{
  My_params my_params= *(My_params *) params;

  FisherGAL_CLUSTER fisher;
  fisher.initExperiment(my_params.survey,"","","","", "",0,1,0.05,0,1);
  cout << "okay " << endl;
  double m=fisher.cal_mlim(my_params.wm,my_params.wl,0,my_params.w,my_params.wa,z);
  cout << m ;
  double result=fisher.NofzFR(my_params.frflag,my_params.wa,my_params.w,my_params.wm, my_params.wl, my_params.wr, my_params.wb, my_params.sigma8, my_params.ns, z, my_params.m1, m, my_params.Bm0, my_params.sigma0, my_params.alpha, my_params.beta);
  cout << "\t" << result << endl;
  //double result=fisher.NofzFR(my_params,z);
  return result;

}

double setNFR(double z, My_params *c1, int iflag)
{
  static My_params *c;
  if (iflag==1) {
    c=c1;
    return 0.0;
  }

  /*double NNN=0,nn=1, result;
  double omnu = 0.0;

  My_params my_params=*c;
  my_params.z=z;

  setNofzFR(0.0,&my_params,1);
  result= qromb(getNofzFR,1e13,1e17,1e-5);
  FisherGAL_CLUSTER fisher;
  NNN=fisher.ComovingVolume(z,c->wm,0,c->wl,c->w)/2 * result;*/
  FisherGAL_CLUSTER fisher;
  double NNN=fisher.NofzFR(c->frflag,c->wa,c->w,c->wm, c->wl, c->wr, c->wb, c->sigma8, c->ns, z, c->m1, c->m2, c->Bm0, c->sigma0, c->alpha, c->beta);

  return NNN;

}

double getNFR(double z)
{
  return setNFR(z,NULL,0);
}


// Total number clusters in the redshift bin (zzz1,zzz2), and mass bin (Mmm1,Mmm2) in the f(R) model
double FisherGAL_CLUSTER::NFR(bool frflag,double h, double wa,double w,double wmmmhh, double wlll, double wbbbhh,double wrrr,  double sigmaa8, double nspecc, double zzz1, double zzz2,double Mmm1,double Mmm2, double Bm0, double sigma0,double alpha, double beta)
{
  gsl_integration_workspace *ww= gsl_integration_workspace_alloc(1000);
  double result, error;

  double zz= zzz1, intz= (zzz2-zzz1)/50., factor, ml;
  for (int k=0; k<50; k++) {
    ml=cal_mlim(wmmmhh,wlll,wrrr,w,wa,zz);
    result += NofzFR(frflag,wa,w,wmmmhh,wlll,wbbbhh,wrrr,sigmaa8,nspecc,zz,Mmm1,ml,Bm0,sigma0,alpha,beta) * intz;
    zz +=intz ;
  }

  return result;
}


void FisherGAL_CLUSTER::derivatives()
{
  if (dFlag==1)
    {
      double zz;
      if (survey == 0 || survey >=6) zz=intzz;
      else zz=0.15;
      double Mmax=5e16, Mmin=1e13;
      double pstep, p_step;
      double No;
      double cpw0, cmw0,cpwa, cmwa, cmomm,cpomm,cmomb,cpomb,cmoml,cpoml,cmsigma8,cpsigma8, cmnscal,cpnscal;
      double ML=0,MMM,MM,mmm=0,mm=Mmin,intz=intzz, intm=pow(10,0.3),Nm, Np, Nfr;
      double beta=0,alpha=0,sigma0=0.1,Bm0=0.;//Bm0=-0.15;
      double cmBm0,cpBm0,cmSigma0,cpSigma0,cmAlpha,cpAlpha,cmBeta,cpBeta;
      nparamFC=12;
      ifstream fin;
      ofstream fout,fout2;
      double *dNN;
      CosmParam *c;
      c=&fiducialFC;
      double h=c->h;

      dNN=dvector(1,nparamFC);

      pstep=0.001;

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

      cmBm0=Bm0- pstep;
      cpBm0=Bm0+ pstep;
      //cmBm0=Bm0- Bm0*pstep;
      //cpBm0=Bm0+ Bm0*pstep;

      cmSigma0=sigma0- sigma0*pstep;
      cpSigma0=sigma0+ sigma0*pstep;

      cmAlpha=alpha - pstep;   //since alpha and beta are zero.
      cpAlpha=alpha + pstep;

      cmBeta=beta - pstep;
      cpBeta=beta + pstep;

      fout.open(CLL_name.c_str());

      while ( zz < zmaxx ) {
        mm=cal_mlim(c->omm*h*h,c->oml,0,c->w0,c->w1,zz);
        while ( mm < Mmax ) {
	//MMM=mmm/h;
	//MM=mm/h;
        double ml=mm;
        //MM=Mmax;
        MM=mm*intm;
        double factor=findMyn(zz,ml, c->omm);
        //double factor= ml;
	    No= NofzFR(0,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,zz,MM,ml,Bm0,sigma0,alpha,beta);
	    Nfr=NofzFR(1,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,zz,factor*intm,factor,Bm0,sigma0,alpha,beta);
            cout<< zz <<"\t" << ml  << "\t" << MM << "\t" << Nfr<<"\t"<<No<<endl;
            dNN[1]=(NofzFR(0,c->w1,c->w0,cpomm,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,zz,MM,ml,Bm0,sigma0,alpha,beta)-NofzFR(0,c->w1,c->w0,cmomm,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,zz,MM,ml,Bm0,sigma0,alpha,beta))/(2*c->omm*h*h*(pstep));
            dNN[2]=(NofzFR(0,c->w1,c->w0,c->omm*h*h,cpoml,c->omb*h*h,0,c->sigma8, c->nscal,zz,MM,ml,Bm0,sigma0,alpha,beta)-NofzFR(0,c->w1,c->w0,c->omm*h*h,cmoml,c->omb*h*h,0,c->sigma8, c->nscal,zz,MM,ml,Bm0,sigma0,alpha,beta))/(2*c->oml*(pstep));
            dNN[3]=(NofzFR(0,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,cpsigma8, c->nscal,zz,MM,ml,Bm0,sigma0,alpha,beta)-NofzFR(0,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,cmsigma8, c->nscal,zz,MM,ml,Bm0,sigma0,alpha,beta))/(2*c->sigma8*(pstep));
            dNN[4]=(NofzFR(0,c->w1,c->w0,c->omm*h*h,c->oml,cpomb,0,c->sigma8, c->nscal,zz,MM,ml,Bm0,sigma0,alpha,beta)-NofzFR(0,c->w1,c->w0,c->omm*h*h,c->oml,cmomb,0,c->sigma8, c->nscal,zz,MM,ml,Bm0,sigma0,alpha,beta))/(2*c->omb*h*h*pstep);
            dNN[5]=(NofzFR(0,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, cpnscal,zz,MM,ml,Bm0,sigma0,alpha,beta)-NofzFR(0,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, cmnscal,zz,MM,ml,Bm0,sigma0,alpha,beta))/(2*c->nscal*pstep);
            dNN[6]=(NofzFR(0,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,zz,MM,ml,cpBm0,sigma0,alpha,beta)-NofzFR(0,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,zz,MM,ml,cmBm0,sigma0,alpha,beta))/(2*pstep);
            dNN[7]=(NofzFR(0,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,zz,MM,ml,Bm0,cpSigma0,alpha,beta)-NofzFR(0,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,zz,MM,ml,Bm0,cmSigma0,alpha,beta))/(2*sigma0*pstep);
            dNN[8]=(NofzFR(0,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,zz,MM,ml,Bm0+1e-4,sigma0,cpAlpha,beta)-NofzFR(0,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,zz,MM,ml,Bm0+1e-3,sigma0,cmAlpha,beta))/(2*pstep);
            dNN[9]=(NofzFR(0,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,zz,MM,ml,Bm0,sigma0,alpha,cpBeta)-NofzFR(0,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,zz,MM,ml,Bm0,sigma0,alpha,cmBeta))/(2*pstep);
            dNN[10]=(Nfr-No)/(findlambdaC(fr_0));
            dNN[11]=(NofzFR(0,c->w1,cpw0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,zz,MM,ml,Bm0,sigma0,alpha,beta)-NofzFR(0,c->w1,cmw0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,zz,MM,ml,Bm0,sigma0,alpha,beta))/(2*c->w0*pstep/10.0);
            dNN[12]=(NofzFR(0,cpwa,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,zz,MM,ml,Bm0,sigma0,alpha,beta)-NofzFR(0,cmwa,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,zz,MM,ml,Bm0,sigma0,alpha,beta))/(2*pstep);

	    /*No= NFR(0, h,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,zz,zmaxx, Mmax,ml,Bm0,sigma0,alpha,beta);
	    Nfr=NFR(1,h, c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,zz,zmaxx, Mmax,ml,Bm0,sigma0,alpha,beta);
            cout<< zz <<"\t" << ml  << "\t" << MM << "\t" << Nfr<<"\t"<<No<<endl;
            dNN[1]=(NFR(0, h,c->w1,c->w0,cpomm,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,zz,zmaxx,Mmax,ml,Bm0,sigma0,alpha,beta)-NFR(0, h,c->w1,c->w0,cmomm,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,zz,zmaxx,Mmax,ml,Bm0,sigma0,alpha,beta))/(2*c->omm*h*h*(pstep));
            dNN[2]=(NFR(0, h,c->w1,c->w0,c->omm*h*h,cpoml,c->omb*h*h,0,c->sigma8, c->nscal,zz,zmaxx,Mmax,ml,Bm0,sigma0,alpha,beta)-NFR(0, h,c->w1,c->w0,c->omm*h*h,cmoml,c->omb*h*h,0,c->sigma8, c->nscal,zz,zmaxx,Mmax,ml,Bm0,sigma0,alpha,beta))/(2*c->oml*(pstep));
            dNN[3]=(NFR(0, h,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,cpsigma8, c->nscal,zz,zmaxx,Mmax,ml,Bm0,sigma0,alpha,beta)-NFR(0, h,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,cmsigma8, c->nscal,zz,zmaxx,Mmax,ml,Bm0,sigma0,alpha,beta))/(2*c->sigma8*(pstep));
            dNN[4]=(NFR(0, h,c->w1,c->w0,c->omm*h*h,c->oml,cpomb,0,c->sigma8, c->nscal,zz,zmaxx,Mmax,ml,Bm0,sigma0,alpha,beta)-NFR(0, h,c->w1,c->w0,c->omm*h*h,c->oml,cmomb,0,c->sigma8, c->nscal,zz,zmaxx,Mmax,ml,Bm0,sigma0,alpha,beta))/(2*c->omb*h*h*pstep);
            dNN[5]=(NFR(0, h,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, cpnscal,zz,zmaxx,Mmax,ml,Bm0,sigma0,alpha,beta)-NFR(0, h,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, cmnscal,zz,zmaxx,Mmax,ml,Bm0,sigma0,alpha,beta))/(2*c->nscal*pstep);
            dNN[6]=(NFR(0, h,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,zz,zmaxx,Mmax,ml,cpBm0,sigma0,alpha,beta)-NFR(0, h,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,zz,zmaxx,Mmax,ml,cmBm0,sigma0,alpha,beta))/(2*pstep);
            dNN[7]=(NFR(0, h,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,zz,zmaxx,Mmax,ml,Bm0,cpSigma0,alpha,beta)-NFR(0, h,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,zz,zmaxx,Mmax,ml,Bm0,cmSigma0,alpha,beta))/(2*sigma0*pstep);
            dNN[8]=(NFR(0, h,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,zz,zmaxx,Mmax,ml,Bm0+1e-4,sigma0,cpAlpha,beta)-NFR(0, h,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,zz,zmaxx,Mmax,ml,Bm0+1e-3,sigma0,cmAlpha,beta))/(2*pstep);
            dNN[9]=(NFR(0, h,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,zz,zmaxx,Mmax,ml,Bm0,sigma0,alpha,cpBeta)-NFR(0, h,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,zz,zmaxx,Mmax,ml,Bm0,sigma0,alpha,cmBeta))/(2*pstep);
            dNN[10]=(Nfr-No)/(findlambdaC(fr_0));
            dNN[11]=(NFR(0, h,c->w1,cpw0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,zz,zmaxx,Mmax,ml,Bm0,sigma0,alpha,beta)-NFR(0, h,c->w1,cmw0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,zz,zmaxx,Mmax,ml,Bm0,sigma0,alpha,beta))/(2*c->w0*pstep/10.0);
            dNN[12]=(NFR(0, h,cpwa,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,zz,zmaxx,Mmax,ml,Bm0,sigma0,alpha,beta)-NFR(0, h,cmwa,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,zz,zmaxx,Mmax,ml,Bm0,sigma0,alpha,beta))/(2*pstep);
           */
            if (dNN[10] < 0) {dNN[10]=0;}
            //fout << zz <<"\t"<< mm << "\t" << No <<"\t"<< dNN[1] <<"\t"<< dNN[2] <<"\t"<<dNN[3]  <<"\t"<<dNN[4]<<"\t"<<dNN[5]<<"\t"<<dNN[6]<<"\t"<<dNN[7]<<"\t"<<dNN[8]<<"\t"<<dNN[9]<<"\t"<<dNN[10]  << endl;
            fout << zz <<"\t"<< MM << "\t" << No <<"\t"<< dNN[1] <<"\t"<< dNN[2] <<"\t"<<dNN[3]  <<"\t"<<dNN[4]<<"\t"<<dNN[5]
                 <<"\t"<<dNN[6]<<"\t"<<dNN[7]<<"\t"<<dNN[8]<<"\t"<<dNN[9]<<"\t"<<dNN[10]  <<"\t"<<dNN[11] << "\t" << dNN[12] << endl;
        //}
	////mmm=mm;
        mm *= intm;
        //mm = Mmax;
      }
      zz += intz;
      //mm = Mmin;
    }
    fout.close();
    free_dvector(dNN,1,nparamFC);
    }
}

void FisherGAL_CLUSTER::fisherMatrix()
{
  double **fisher, **ifisher;
  ifstream fin;
  double NNN, fM, zzz, MMM;
  double *dnn;
  CosmParam *c;

  c=&fiducialFC;
  nparamFC=12;
  dnn=dvector(1,nparamFC);
  fisher=dmatrix(1,nparamFC,1,nparamFC);
  ifisher=dmatrix(1,nparamFC,1,nparamFC);


  fin.open(CLL_name.c_str());

  for(int i=1;i<nparamFC+1;i++){
    for(int j=1;j<=i;j++){
      fM=0.0;
      fin.clear() ;
      fin.seekg(0, ios::beg) ;
      //fin>>zzz>>MMM>>NNN>>dnn[1]>>dnn[2]>>dnn[3]>>dnn[4]>>dnn[5]>>dnn[6]>>dnn[7]>>dnn[8]>>dnn[9]>>dnn[10];
      fin>>zzz>>MMM>>NNN>>dnn[1]>>dnn[2]>>dnn[3]>>dnn[4]>>dnn[5]>>dnn[6]>>dnn[7]>>dnn[8]>>dnn[9]>>dnn[10]>>dnn[11]>> dnn[12];
      while (!fin.eof()){
	//fin>>zzz>>MMM>>NNN>>dnn[1]>>dnn[2]>>dnn[3]>>dnn[4]>>dnn[5]>>dnn[6]>>dnn[7]>>dnn[8]>>dnn[9]>>dnn[10];
	fin>>zzz>>MMM>>NNN>>dnn[1]>>dnn[2]>>dnn[3]>>dnn[4]>>dnn[5]>>dnn[6]>>dnn[7]>>dnn[8]>>dnn[9]>>dnn[10]>>dnn[11]>>dnn[12];
	cout<< "\t" <<zzz<< "\t" <<MMM<< "\t" <<NNN<< "\t" <<dnn[1]<< "\t" <<dnn[2]<< "\t" <<dnn[3]<< "\t" <<dnn[4]<< "\t" <<dnn[5]<< "\t" <<dnn[6]<< "\t" <<dnn[7]<< "\t" <<dnn[8]<< "\t" <<dnn[9]<< "\t" <<dnn[10]<< "\t" <<dnn[11]<< "\t" <<dnn[12]<<endl;
	//if(zzz/intzz-floor(zzz/intzz)==0){  //per cambiare l'intervallo in z (attenzione pero' dovrei combiare anche la moltiplicazione per intz)

	  if (NNN !=0) {
	    //fM+=(1/NNN)*dnn[i]*dnn[j]*SCCC;    //I multiply by the sky coverage and the z bin here, then I can change the parameter and don't have to calculate derivatives all the times
	    cout << SCCC << endl;
            fM+=(1/NNN)*dnn[i]*dnn[j]*intzz*SCCC;    //I multiply by the sky coverage and the z bin here, then I can change the parameter and don't have to calculate derivatives all the times
	  }
	//}
      }
      fisher[i][j]=fM;
      if(i!=j) fisher[j][i]=fM;
      // cout << fisher[i][j] << endl;
    }
  }
  free_dvector(dnn,1,nparamFC);

  //invertM(fisher,nparamFC,ifisher);
  invertMatrix(fisher,nparamFC,ifisher);
  cout<<"Errors:"<<endl;
  cout<<"OmegaM= "<<c->omm<<"\t"<<"OmegaL= "<<c->oml<<"\t"<<"s8= "<<c->sigma8<<"\t"<<"Omegab= "<<c->omb<<"\t"<<"ns= "<<c->nscal<<"\t"<<"Bm= "<<-0.15<<"\t"<<"Sm= "<<0.25<<"\t"<<"alpha= "<<0<<"\t"<<  "beta= "<<0<<"\t"<<"fR= "<<0<<endl;
  //cout<<"OmegaM= "<<c->omm<<"\t"<<"OmegaL= "<<c->oml<<"\t"<<"s8= "<<c->sigma8<<"\t"<<"Omegab= "<<c->omb<<"\t"<<"ns= "<<c->nscal<<"\t"<<"Bm= "<<-0.15<<"\t"<<"Sm= "<<0.25<<"\t"<<"alpha= "<<0<<"\t"<<  "beta= "<<0<<"\t"<< " w= " << c->w0 << "\t"<<"fR= "<<0<<endl;
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

//*****************************************************************
//********* Here I calculate the POWER SPECTRUM derivatives *******
//*****************************************************************

double FisherGAL_CLUSTER::findPk (double z, double k)
{
  ifstream fin;
  double kk0,kk1, nn0,nn1, zz, kint=pow(10,0.017), intzzz=0.02;
  fin.open(FR_Pk.c_str());
  if (fin.is_open()) {
      for (int i=0; i< floor(z/intzzz)*118; i++) fin >> kk0 >> zz >> nn0; // jump to the corresponding z
      for (int i=0; i<118; i++) {
          fin >> kk0 >> zz >> nn0;
          if ( (kk0*kint > k)  && (i !=117)) {
              fin >> kk1 >> zz >> nn1;
              return nn0 + (nn1-nn0)/(log10(kk1/kk0)) * log10(k/kk0) ;
             /* double dx=kk*(kint-1);
              double deltx=k-kk;
              fin >> kk >> zz >> nn2;
              double dy=nn2-nn;
              return nn+deltx*dy/dx;*/
          }
      }
  } else {
    cerr << "Error. Cannot open Pk file" << endl;
    exit(0);
  }

  fin.close();
  return nn0; //the mass range is beyond the tabulated values,scale according to the last pair
}


void FisherGAL_CLUSTER::derivatives_PS()
{
  if (dFlag==1)
    {
      double zz=0.02;
      if (survey == 0 || survey>=4 ) zz=0.02;
      else zz=0.15;
      double pstep;
      double Noii, Noij, Nojj;
      ifstream fin,fin2;
      ofstream fout;//,fout2;
      double *dNN;
      CosmParam *c;

      double cmw1, cpw1, cmw0,cpw0,cmomb,cpomb,cmomm,cpomm,cmoml,cpoml,cmsigma8,cpsigma8, cmnscal,cpnscal,Nfr,ML=0,Vol,MMV1, MMV2, MMVfr,Veff, Vefffr;
      double kk2,BZ0,BZfr,BZ1p,BZ1m,BZ2p,BZ2m,BZ3p,BZ3m,BZ4p,BZ4m,BZ5p,BZ5m,BZ6p,BZ6m,BZ7p,BZ7m,BZ8p,BZ8m,BZ9p,BZ9m,BZ10p,BZ10m, BZ11p, BZ11m;
      double B2Z0,B2Zfr,B2Z1p,B2Z1m,B2Z2p,B2Z2m,B2Z3p,B2Z3m,B2Z4p,B2Z4m,B2Z5p,B2Z5m,B2Z6p,B2Z6m,B2Z7p,B2Z7m,B2Z8p,B2Z8m,B2Z9p,B2Z9m,B2Z10p,B2Z10m, B2Z11p, B2Z11m, Bmdym, Bmdym2, Ball, Ball2;
      double beta=0,alpha=0,sigma0=0.1,Bm0=0.;//Bm0=-0.15;
      double cmBm0,cpBm0,cmSigma0,cpSigma0,cmAlpha,cpAlpha,cmBeta,cpBeta;
      double KK=1e-2,Kmaxx=0.1,kint=pow(10,0.017),ztemp=0, mm1=1e13,mm2=1e13, Mmax=1e16, MM1, MM2,  intm= pow(10,6*0.05) ;
      //double KK=1e-2,Kmaxx=0.15,kint=pow(10,0.017),ztemp=0, intk=5e-3;

      c=&fiducialFC;
      double h=c->h;
      Cosmology mycosmo(c->omm, c->oml, c->omb, c->h, c->sigma8  , c->nscal, 0, c->w0, c->w1);
      //cp=&fiducialFC;
      //cm=&fiducialFC;
      int nparam=12;
      dNN=dvector(1,nparam);

      pstep=0.001;

      cmw1=c->w1-pstep;
      cpw1=c->w1+pstep;

      cmw0=c->w0- c->w0*pstep;
      cpw0=c->w0+ c->w0*pstep;

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

      cmBm0=Bm0- pstep;
      cpBm0=Bm0+ pstep;
      //cmBm0=Bm0- Bm0*pstep;
      //cpBm0=Bm0+ Bm0*pstep;

      cmSigma0=sigma0- sigma0*pstep;
      cpSigma0=sigma0+ sigma0*pstep;

      cmAlpha=alpha - pstep;   //since alpha and beta are zero.
      cpAlpha=alpha + pstep;

      cmBeta=beta - pstep;
      cpBeta=beta + pstep;

      cout<<"Calculating Power Spectrum derivatives"<<endl;
      fout.open(CLL_namePS.c_str());    //qui salvo le derivate del Power Spectrum.

      while (zz < zmaxx) {
        mm1=cal_mlim(c->omm*h*h,c->oml,0,c->w0,c->w1,zz);
        while ( mm1<Mmax) {
          double ml1=mm1;
          MM1=mm1*intm;
          //double factor1=ml1;
          double factor1=findMyn(zz,ml1, c->omm);
          double mdy1=findMyn(zz,ml1, c->omm);
          BZ0=ComBiasFR(0,MM1,ml1,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0,sigma0,alpha,beta,zz);
          Bmdym=ComBiasFR(0,mdy1*intm,mdy1,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0,sigma0,alpha,beta,zz);
          Ball=ComBiasFR(1,mdy1*intm,mdy1,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0,sigma0,alpha,beta,zz);
          BZfr=ComBiasFR(1,factor1*intm,factor1,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0,sigma0,alpha,beta,zz);
          BZ1p=ComBiasFR(0,MM1,ml1,c->w1,c->w0,cpomm,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0,sigma0,alpha,beta,zz);
          BZ1m=ComBiasFR(0,MM1,ml1,c->w1,c->w0,cmomm,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0,sigma0,alpha,beta,zz);
          BZ2p=ComBiasFR(0,MM1,ml1,c->w1,c->w0,c->omm*h*h,cpoml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0,sigma0,alpha,beta,zz);
          BZ2m=ComBiasFR(0,MM1,ml1,c->w1,c->w0,c->omm*h*h,cmoml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0,sigma0,alpha,beta,zz);
          BZ3p=ComBiasFR(0,MM1,ml1,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,cpsigma8, c->nscal,Bm0,sigma0,alpha,beta,zz);
          BZ3m=ComBiasFR(0,MM1,ml1,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,cmsigma8, c->nscal,Bm0,sigma0,alpha,beta,zz);
          BZ4p=ComBiasFR(0,MM1,ml1,c->w1,c->w0,c->omm*h*h,c->oml,cpomb,0,c->sigma8, c->nscal,Bm0,sigma0,alpha,beta,zz);
          BZ4m=ComBiasFR(0,MM1,ml1,c->w1,c->w0,c->omm*h*h,c->oml,cmomb,0,c->sigma8, c->nscal,Bm0,sigma0,alpha,beta,zz);
          BZ5p=ComBiasFR(0,MM1,ml1,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, cpnscal,Bm0,sigma0,alpha,beta,zz);
          BZ5m=ComBiasFR(0,MM1,ml1,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, cmnscal,Bm0,sigma0,alpha,beta,zz);
          BZ6p=ComBiasFR(0,MM1,ml1,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,cpBm0,sigma0,alpha,beta,zz);
          BZ6m=ComBiasFR(0,MM1,ml1,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,cmBm0,sigma0,alpha,beta,zz);
          BZ7p=ComBiasFR(0,MM1,ml1,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0,cpSigma0,alpha,beta,zz);
          BZ7m=ComBiasFR(0,MM1,ml1,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0,cmSigma0,alpha,beta,zz);
          BZ8p=ComBiasFR(0,MM1,ml1,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0,sigma0,cpAlpha,beta,zz);
          BZ8m=ComBiasFR(0,MM1,ml1,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0,sigma0,cmAlpha,beta,zz);
          BZ9p=ComBiasFR(0,MM1,ml1,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0,sigma0,alpha,cpBeta,zz);
          BZ9m=ComBiasFR(0,MM1,ml1,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0,sigma0,alpha,cmBeta,zz);
          BZ10p=ComBiasFR(0,MM1,ml1,c->w1,cpw0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0,sigma0,alpha,beta,zz);
          BZ10m=ComBiasFR(0,MM1,ml1,c->w1,cmw0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0,sigma0,alpha,beta,zz);
          BZ11p=ComBiasFR(0,MM1,ml1,cpw1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0,sigma0,alpha,beta,zz);
          BZ11m=ComBiasFR(0,MM1,ml1,cmw1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0,sigma0,alpha,beta,zz);
          //mm2=cal_mlim(c->omm*h*h,c->oml,0,c->w0,c->w1,zz);
          mm2=mm1;
          while ( mm2<Mmax) {
            double ml2=mm2;
            MM2=mm2*intm;
            //double factor2=ml2;
            double factor2=findMyn(zz,ml2, c->omm);
            double mdy2=findMyn(zz,ml2, c->omm);
            B2Z0=ComBiasFR(0,MM2,ml2,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0,sigma0,alpha,beta,zz);
            Bmdym2=ComBiasFR(0,mdy2*intm,mdy2,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0,sigma0,alpha,beta,zz);
            Ball2=ComBiasFR(1,mdy2*intm,mdy2,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0,sigma0,alpha,beta,zz);
            B2Zfr=ComBiasFR(1,factor2*intm,factor2,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0,sigma0,alpha,beta,zz);
            B2Z1p=ComBiasFR(0,MM2,ml2,c->w1,c->w0,cpomm,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0,sigma0,alpha,beta,zz);
            B2Z1m=ComBiasFR(0,MM2,ml2,c->w1,c->w0,cmomm,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0,sigma0,alpha,beta,zz);
            B2Z2p=ComBiasFR(0,MM2,ml2,c->w1,c->w0,c->omm*h*h,cpoml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0,sigma0,alpha,beta,zz);
            B2Z2m=ComBiasFR(0,MM2,ml2,c->w1,c->w0,c->omm*h*h,cmoml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0,sigma0,alpha,beta,zz);
            B2Z3p=ComBiasFR(0,MM2,ml2,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,cpsigma8, c->nscal,Bm0,sigma0,alpha,beta,zz);
            B2Z3m=ComBiasFR(0,MM2,ml2,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,cmsigma8, c->nscal,Bm0,sigma0,alpha,beta,zz);
            B2Z4p=ComBiasFR(0,MM2,ml2,c->w1,c->w0,c->omm*h*h,c->oml,cpomb,0,c->sigma8, c->nscal,Bm0,sigma0,alpha,beta,zz);
            B2Z4m=ComBiasFR(0,MM2,ml2,c->w1,c->w0,c->omm*h*h,c->oml,cmomb,0,c->sigma8, c->nscal,Bm0,sigma0,alpha,beta,zz);
            B2Z5p=ComBiasFR(0,MM2,ml2,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, cpnscal,Bm0,sigma0,alpha,beta,zz);
            B2Z5m=ComBiasFR(0,MM2,ml2,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, cmnscal,Bm0,sigma0,alpha,beta,zz);
            B2Z6p=ComBiasFR(0,MM2,ml2,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,cpBm0,sigma0,alpha,beta,zz);
            B2Z6m=ComBiasFR(0,MM2,ml2,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,cmBm0,sigma0,alpha,beta,zz);
            B2Z7p=ComBiasFR(0,MM2,ml2,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0,cpSigma0,alpha,beta,zz);
            B2Z7m=ComBiasFR(0,MM2,ml2,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0,cmSigma0,alpha,beta,zz);
            B2Z8p=ComBiasFR(0,MM2,ml2,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0+1e-4,sigma0,cpAlpha,beta,zz);
            B2Z8m=ComBiasFR(0,MM2,ml2,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0+1e-4,sigma0,cmAlpha,beta,zz);
            B2Z9p=ComBiasFR(0,MM2,ml2,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0,sigma0,alpha,cpBeta,zz);
            B2Z9m=ComBiasFR(0,MM2,ml2,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0,sigma0,alpha,cmBeta,zz);
            B2Z10p=ComBiasFR(0,MM2,ml2,c->w1,cpw0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0,sigma0,alpha,beta,zz);
            B2Z10m=ComBiasFR(0,MM2,ml2,c->w1,cmw0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0,sigma0,alpha,beta,zz);
            B2Z11p=ComBiasFR(0,MM2,ml2,cpw1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0,sigma0,alpha,beta,zz);
            B2Z11m=ComBiasFR(0,MM2,ml2,cmw1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0,sigma0,alpha,beta,zz);
      //while (zz < zmaxx) {
       //double factor=ML;
            Vol=ComovingVolume(zz,c->omm*h*h,0,c->oml,c->w0, c->w1)*SCCC*intzz;
            if ( zz == 0 ) {
              MMV1 = 0;
              MMV2 = 0;
              MMVfr = 0;
            } else {
              MMV1=NofzFR(0,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,zz,MM1,ml1,Bm0,sigma0,alpha,beta)/(Vol/(2*SCCC*intzz));
              MMV2=NofzFR(0,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,zz,MM2,ml2,Bm0,sigma0,alpha,beta)/(Vol/(2*SCCC*intzz));
              MMVfr=NofzFR(1,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,zz,1e16,factor1,Bm0,sigma0,alpha,beta)/(Vol/(2*SCCC*intzz));
            }
            while ( KK < Kmaxx ) {
              dNN[1]=KK*(log(BZ1p*B2Z1p*PowSFR(0,c->w1,c->w0, cpomm,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0,sigma0,alpha,beta,zz,KK))-log(BZ1m*B2Z1m*PowSFR(0,c->w1,c->w0,cmomm,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0,sigma0,alpha,beta,zz,KK)))/(2*c->omm*h*h*(pstep));
              dNN[2]=KK*(log(BZ2p*B2Z2p*PowSFR(0,c->w1,c->w0,c->omm*h*h,cpoml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0,sigma0,alpha,beta,zz,KK))-log(BZ2m*B2Z2m*PowSFR(0,c->w1,c->w0,c->omm*h*h,cmoml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0,sigma0,alpha,beta,zz,KK)))/(2*c->oml*(pstep));
              dNN[3]=KK*(log(BZ3p*B2Z3p*PowSFR(0,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,cpsigma8, c->nscal,Bm0,sigma0,alpha,beta,zz,KK))-log(BZ3m*B2Z3m*PowSFR(0,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,cmsigma8, c->nscal,Bm0,sigma0,alpha,beta,zz,KK)))/(2*c->sigma8*(pstep));
              dNN[4]=KK*(log(BZ4p*B2Z4p*PowSFR(0,c->w1,c->w0,c->omm*h*h,c->oml,cpomb,0,c->sigma8, c->nscal,Bm0,sigma0,alpha,beta,zz,KK))-log(BZ4m*B2Z4m*PowSFR(0,c->w1,c->w0,c->omm*h*h,c->oml,cmomb,0,c->sigma8, c->nscal,Bm0,sigma0,alpha,beta,zz,KK)))/(2*c->omb*h*h*(pstep));
              dNN[5]=KK*(log(BZ5p*B2Z5p*PowSFR(0,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, cpnscal,Bm0,sigma0,alpha,beta,zz,KK))-log(BZ5m*B2Z5m*PowSFR(0,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, cmnscal,Bm0,sigma0,alpha,beta,zz,KK)))/(2*c->nscal*(pstep));
              dNN[6]=KK*(log(BZ6p*B2Z6p*PowSFR(0,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,cpBm0,sigma0,alpha,beta,zz,KK))-log(BZ6m*B2Z6m*PowSFR(0,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,cmBm0,sigma0,alpha,beta,zz,KK)))/(2*pstep);
              dNN[7]=KK*(log(BZ7p*B2Z7p*PowSFR(0,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0,cpSigma0,alpha,beta,zz,KK))-log(BZ7m*B2Z7m*PowSFR(0,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0,cmSigma0,alpha,beta,zz,KK)))/(2*sigma0*pstep);
              dNN[8]=KK*(log(BZ8p*B2Z8p*PowSFR(0,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0+1e-4,sigma0,cpAlpha,beta,zz,KK))-log(BZ8m*B2Z8m*PowSFR(0,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0+1e-4,sigma0,cmAlpha,beta,zz,KK)))/(2*pstep);
              dNN[9]=KK*(log(BZ9p*B2Z9p*PowSFR(0,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0,sigma0,alpha,cpBeta,zz,KK))-log(BZ9m*B2Z9m*PowSFR(0,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0,sigma0,alpha,cmBeta,zz,KK)))/(2*pstep);
              dNN[11]=KK*(log(BZ10p*B2Z10p*PowSFR(0,c->w1,cpw0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0,sigma0,alpha,beta,zz,KK))-log(BZ10m*B2Z10m*PowSFR(0,c->w1,cmw0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0,sigma0,alpha,beta,zz,KK)))/(2*c->w0*(pstep));
              dNN[12]=KK*(log(BZ11p*B2Z11p*PowSFR(0,cpw1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0,sigma0,alpha,beta,zz,KK))-log(BZ11m*B2Z11m*PowSFR(0,cmw1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0,sigma0,alpha,beta,zz,KK)))/(2*pstep);
              Noii=(BZ0*BZ0*PowSFR(0,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0,sigma0,alpha,beta,zz,KK));
              Nojj=(B2Z0*B2Z0*PowSFR(0,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0,sigma0,alpha,beta,zz,KK));
              Noij=(BZ0*B2Z0*PowSFR(0,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0,sigma0,alpha,beta,zz,KK));
              double psfr=PowSFR(0,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0,sigma0,alpha,beta,zz,KK);
              double psfr2=PowSFR(1,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0,sigma0,alpha,beta,zz,KK);
              Nfr=log(BZfr*B2Zfr*psfr2);
              if (ml1 == ml2) {
                    Veff=Noij*Noij*MMV1*MMV2/((1+Noii*MMV1)*(1+Nojj*MMV2)+MMV2*MMV1*pow(Noij+1.0/MMV1,2));
              } else {
                    Veff=Noij*Noij*MMV1*MMV2/((1+Noii*MMV1)*(1+Nojj*MMV2)+MMV2*MMV1*pow(Noij,2));
              }
              //Veff=pow(((No)*MMV)/(1+(No)*MMV),2);
              Vefffr=pow(((Nfr)*MMVfr)/(1+(Nfr)*MMVfr),2);
              dNN[10]=KK*(Nfr-log(Noij))/(findlambdaC(fr_0));
              //cout<<zz << "\t" << KK << "\t" << ml1 << "\t" << MM1 << "\t" << ml2 << "\t" << MM2 << "\t" << Vol << "\t" << Veff << endl;
              cout<<zz << "\t" << KK << "\t" << ml1 << "\t" <<  ml2 <<  "\t" << Ball*Ball2*psfr2 << "\t" << BZfr*B2Zfr*psfr <<"\t" << BZ0*B2Z0*psfr2 << "\t" <<  Bmdym*Bmdym2*psfr<<"\t" <<Noij << endl;
              fout << KK <<"\t"<<zz<<"\t"<< Noii <<"\t"<<Vol<<"\t"<<Veff<<"\t"<<ml1<<"\t"<<ml2<<"\t"
                  << dNN[1] <<"\t"<< dNN[2] <<"\t"<<dNN[3]  <<"\t"<<dNN[4] <<"\t"<< dNN[5] <<"\t"<< dNN[6] <<"\t"<<dNN[7]  <<"\t"<<dNN[8] <<"\t"<<dNN[9] <<"\t"<<dNN[10] << "\t" << dNN[11]  <<"\t"<< dNN[12]<< endl;

              //KK += intk;   //logarithmic deltaK, the deltaK is constant, I change just Kmax
              KK *= kint;   //logarithmic deltaK, the deltaK is constant, I change just Kmax
            }
            KK =1e-2;
            mm2 *= intm;
          }
          mm1 *= intm;
        }
        zz += 0.02;
        }
      fout.close();
      free_dvector(dNN,1,nparam);

    }
}

/*void FisherGAL_CLUSTER::derivatives_PS()
{
  if (dFlag==1)
    {
      double zz=0.02;
      double pstep;
      double No;
      ifstream fin,fin2;
      ofstream fout;//,fout2;
      double *dNN;
      CosmParam *c;

      double cmw1, cpw1, cmw0,cpw0,cmomb,cpomb,cmomm,cpomm,cmoml,cpoml,cmsigma8,cpsigma8, cmnscal,cpnscal,Nfr,ML=0,Vol,MMV,MMVfr,Veff, Vefffr;
      double kk2,BZ0,BZfr,BZ1p,BZ1m,BZ2p,BZ2m,BZ3p,BZ3m,BZ4p,BZ4m,BZ5p,BZ5m,BZ6p,BZ6m,BZ7p,BZ7m,BZ8p,BZ8m,BZ9p,BZ9m,BZ10p,BZ10m, BZ11p, BZ11m;
      double beta=0,alpha=0,sigma0=0.1,Bm0=0.;//Bm0=-0.15;
      double cmBm0,cpBm0,cmSigma0,cpSigma0,cmAlpha,cpAlpha,cmBeta,cpBeta;
      double KK=1e-2,Kmaxx=0.1,kint=pow(10,0.017),ztemp=0;
      //double KK=1e-2,Kmaxx=0.15,kint=pow(10,0.017),ztemp=0, intk=5e-3;

      c=&fiducialFC;
      double h=c->h;
      Cosmology mycosmo(c->omm, c->oml, c->omb, c->h, c->sigma8  , c->nscal, 0, c->w0, c->w1);
      //cp=&fiducialFC;
      //cm=&fiducialFC;
      int nparam=12;
      dNN=dvector(1,nparam);

      pstep=0.001;

      cmw1=c->w1-pstep;
      cpw1=c->w1+pstep;

      cmw0=c->w0- c->w0*pstep;
      cpw0=c->w0+ c->w0*pstep;

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

      cmBm0=Bm0- pstep;
      cpBm0=Bm0+ pstep;
      //cmBm0=Bm0- Bm0*pstep;
      //cpBm0=Bm0+ Bm0*pstep;

      cmSigma0=sigma0- sigma0*pstep;
      cpSigma0=sigma0+ sigma0*pstep;

      cmAlpha=alpha - pstep;   //since alpha and beta are zero.
      cpAlpha=alpha + pstep;

      cmBeta=beta - pstep;
      cpBeta=beta + pstep;

      cout<<"Calculating Power Spectrum derivatives"<<endl;
      //fin2.open(BZfr_name.c_str());     //qui ci metto il file di b(z)_fr
      //fin.open(FR_Pk.c_str());
      fout.open(CLL_namePS.c_str());    //qui salvo le derivate del Power Spectrum.

      //for (int i=1;i<=zmaxx/intzz;i++){
      while (zz < zmaxx) {
	//fin2 >> zz  >> BZ_fr;
        ML=cal_mlim(c->omm*h*h,c->oml,0,c->w0,c->w1,zz);
	BZ0=ComBiasFR(0,1e16,ML,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0,sigma0,alpha,beta,zz);
	BZfr=ComBiasFR(1,1e16,ML,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0,sigma0,alpha,beta,zz);
	BZ1p=ComBiasFR(0,1e16,ML,c->w1,c->w0,cpomm,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0,sigma0,alpha,beta,zz);
	BZ1m=ComBiasFR(0,1e16,ML,c->w1,c->w0,cmomm,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0,sigma0,alpha,beta,zz);
	BZ2p=ComBiasFR(0,1e16,ML,c->w1,c->w0,c->omm*h*h,cpoml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0,sigma0,alpha,beta,zz);
	BZ2m=ComBiasFR(0,1e16,ML,c->w1,c->w0,c->omm*h*h,cmoml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0,sigma0,alpha,beta,zz);
	BZ3p=ComBiasFR(0,1e16,ML,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,cpsigma8, c->nscal,Bm0,sigma0,alpha,beta,zz);
	BZ3m=ComBiasFR(0,1e16,ML,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,cmsigma8, c->nscal,Bm0,sigma0,alpha,beta,zz);
	BZ4p=ComBiasFR(0,1e16,ML,c->w1,c->w0,c->omm*h*h,c->oml,cpomb,0,c->sigma8, c->nscal,Bm0,sigma0,alpha,beta,zz);
	BZ4m=ComBiasFR(0,1e16,ML,c->w1,c->w0,c->omm*h*h,c->oml,cmomb,0,c->sigma8, c->nscal,Bm0,sigma0,alpha,beta,zz);
	BZ5p=ComBiasFR(0,1e16,ML,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, cpnscal,Bm0,sigma0,alpha,beta,zz);
	BZ5m=ComBiasFR(0,1e16,ML,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, cmnscal,Bm0,sigma0,alpha,beta,zz);
	BZ6p=ComBiasFR(0,1e16,ML,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,cpBm0,sigma0,alpha,beta,zz);
	BZ6m=ComBiasFR(0,1e16,ML,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,cmBm0,sigma0,alpha,beta,zz);
	BZ7p=ComBiasFR(0,1e16,ML,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0,cpSigma0,alpha,beta,zz);
	BZ7m=ComBiasFR(0,1e16,ML,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0,cmSigma0,alpha,beta,zz);
	BZ8p=ComBiasFR(0,1e16,ML,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0+1e-4,sigma0,cpAlpha,beta,zz);
	BZ8m=ComBiasFR(0,1e16,ML,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0+1e-4,sigma0,cmAlpha,beta,zz);
	BZ9p=ComBiasFR(0,1e16,ML,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0,sigma0,alpha,cpBeta,zz);
	BZ9m=ComBiasFR(0,1e16,ML,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0,sigma0,alpha,cmBeta,zz);
	BZ10p=ComBiasFR(0,1e16,ML,c->w1,cpw0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0,sigma0,alpha,beta,zz);
	BZ10m=ComBiasFR(0,1e16,ML,c->w1,cmw0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0,sigma0,alpha,beta,zz);
	BZ11p=ComBiasFR(0,1e16,ML,cpw1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0,sigma0,alpha,beta,zz);
	BZ11m=ComBiasFR(0,1e16,ML,cmw1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0,sigma0,alpha,beta,zz);
      //while (zz < zmaxx) {
       //double factor=ML;
        double factor=findMyn(zz,ML, c->omm);
	Vol=ComovingVolume(zz,c->omm*h*h,0,c->oml,c->w0, c->w1)*SCCC*intzz;
	if ( zz == 0 ) {
          MMV = 0;
          MMVfr = 0;
        } else {
          MMV=NofzFR(0,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,zz,1e16,ML,Bm0,sigma0,alpha,beta)/(Vol/(2*SCCC*intzz));
          MMVfr=NofzFR(1,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,zz,1e16,factor,Bm0,sigma0,alpha,beta)/(Vol/(2*SCCC*intzz));
        }
        while ( KK < Kmaxx ) {
          dNN[1]=KK*(log(BZ1p*PowSFR(0,c->w1,c->w0, cpomm,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0,sigma0,alpha,beta,zz,KK))-log(BZ1m*PowSFR(0,c->w1,c->w0,cmomm,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0,sigma0,alpha,beta,zz,KK)))/(2*c->omm*h*h*(pstep));
	  dNN[2]=KK*(log(BZ2p*PowSFR(0,c->w1,c->w0,c->omm*h*h,cpoml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0,sigma0,alpha,beta,zz,KK))-log(BZ2m*PowSFR(0,c->w1,c->w0,c->omm*h*h,cmoml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0,sigma0,alpha,beta,zz,KK)))/(2*c->oml*(pstep));
	  dNN[3]=KK*(log(BZ3p*PowSFR(0,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,cpsigma8, c->nscal,Bm0,sigma0,alpha,beta,zz,KK))-log(BZ3m*PowSFR(0,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,cmsigma8, c->nscal,Bm0,sigma0,alpha,beta,zz,KK)))/(2*c->sigma8*(pstep));
	  dNN[4]=KK*(log(BZ4p*PowSFR(0,c->w1,c->w0,c->omm*h*h,c->oml,cpomb,0,c->sigma8, c->nscal,Bm0,sigma0,alpha,beta,zz,KK))-log(BZ4m*PowSFR(0,c->w1,c->w0,c->omm*h*h,c->oml,cmomb,0,c->sigma8, c->nscal,Bm0,sigma0,alpha,beta,zz,KK)))/(2*c->omb*h*h*(pstep));
	  dNN[5]=KK*(log(BZ5p*PowSFR(0,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, cpnscal,Bm0,sigma0,alpha,beta,zz,KK))-log(BZ5m*PowSFR(0,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, cmnscal,Bm0,sigma0,alpha,beta,zz,KK)))/(2*c->nscal*(pstep));
	  dNN[6]=KK*(log(BZ6p*PowSFR(0,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,cpBm0,sigma0,alpha,beta,zz,KK))-log(BZ6m*PowSFR(0,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,cmBm0,sigma0,alpha,beta,zz,KK)))/(2*pstep);
	  dNN[7]=KK*(log(BZ7p*PowSFR(0,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0,cpSigma0,alpha,beta,zz,KK))-log(BZ7m*PowSFR(0,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0,cmSigma0,alpha,beta,zz,KK)))/(2*sigma0*pstep);
	  dNN[8]=KK*(log(BZ8p*PowSFR(0,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0,sigma0,cpAlpha,beta,zz,KK))-log(BZ8m*PowSFR(0,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0,sigma0,cmAlpha,beta,zz,KK)))/(2*pstep);
	  dNN[9]=KK*(log(BZ9p*PowSFR(0,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0,sigma0,alpha,cpBeta,zz,KK))-log(BZ9m*PowSFR(0,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0,sigma0,alpha,cmBeta,zz,KK)))/(2*pstep);
	  dNN[11]=KK*(log(BZ10p*PowSFR(0,c->w1,cpw0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0,sigma0,alpha,beta,zz,KK))-log(BZ10m*PowSFR(0,c->w1,cmw0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0,sigma0,alpha,beta,zz,KK)))/(2*c->w0*(pstep));
	  dNN[12]=KK*(log(BZ11p*PowSFR(0,cpw1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0,sigma0,alpha,beta,zz,KK))-log(BZ11m*PowSFR(0,cmw1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0,sigma0,alpha,beta,zz,KK)))/(2*pstep);
          No=(BZ0*PowSFR(0,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0,sigma0,alpha,beta,zz,KK));
          double psfr=PowSFR(0,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0,sigma0,alpha,beta,zz,KK);
          double psfr2=PowSFR(1,c->w1,c->w0,c->omm*h*h,c->oml,c->omb*h*h,0,c->sigma8, c->nscal,Bm0,sigma0,alpha,beta,zz,KK);
          Nfr=log(BZfr*psfr2);
          Veff=pow(((No)*MMV)/(1+(No)*MMV),2);
          Vefffr=pow(((Nfr)*MMVfr)/(1+(Nfr)*MMVfr),2);
	  dNN[10]=KK*(Nfr-log(No))/(findlambdaC(fr_0));
          ////dNN[10]=KK*log(BZfr*findPk(zz,KK)/BZ0)/fr_0;
          //cout<<zz << "\t" << KK << "\t" << BZ0 <<"\t"<< BZfr << "\t"  << psfr << "\t" << psfr2 << endl;
	 // if (No-MMV < 1e-1 )
            cout << zz << "\t" << KK << "\t" <<No  << "\t" << Nfr << endl;
          fout << KK <<"\t"<<zz<<"\t"<< No <<"\t"<<Vol<<"\t"<<Veff<<"\t"<<Vefffr<<"\t"
               //<< dNN[1] <<"\t"<< dNN[2] <<"\t"<<dNN[3]  <<"\t"<<dNN[4] <<"\t"<< dNN[5] <<"\t"<< dNN[6] <<"\t"<<dNN[7]  <<"\t"<<dNN[8] <<"\t"<<dNN[9] <<"\t"<<dNN[10] <<endl;
               << dNN[1] <<"\t"<< dNN[2] <<"\t"<<dNN[3]  <<"\t"<<dNN[4] <<"\t"<< dNN[5] <<"\t"<< dNN[6] <<"\t"<<dNN[7]  <<"\t"<<dNN[8] <<"\t"<<dNN[9] <<"\t"<<dNN[10] << "\t" << dNN[11]  <<"\t"<< dNN[12]<< endl;

	  //KK += intk;   //logarithmic deltaK, the deltaK is constant, I change just Kmax
	  KK *= kint;   //logarithmic deltaK, the deltaK is constant, I change just Kmax
	}
        KK =1e-2;
        zz += 0.02;
        }

      fout.close();
      free_dvector(dNN,1,nparam);

    }
}*/

/*double FisherGAL_CLUSTER::PowS(double wmmm, double wlll, double wbbb,double wrrr,  double sigmaa8, double nspecc, double zzz, double KK)
{
  double D,PS;
  double omnu = 0.0;
  Cosmology my_cosmo(wmmm+wrrr, wlll, wbbb, h, sigmaa8, nspecc, omnu);
  D=my_cosmo.growthFac(zzz);
  PS=pow(D,2)*my_cosmo.powerSpectrum(KK);

  return PS;

}*/

double FisherGAL_CLUSTER::PowSFR(bool frflag, double wa, double w, double wmmmhh, double wlll, double wbbbhh, double wrrr,  double sigmaa8, double nspecc, double Bm0, double sigma0, double alpha, double beta, double zzz, double KK)
{
  double D,PS,nn=1,beff;
  double omnu = 0.0;
  double h=sqrt(wmmmhh/(1-wlll));

  if (frflag) nn=findPk(zzz,KK);
  Cosmology my_cosmo(wmmmhh/(h*h), wlll, wbbbhh/(h*h), h, sigmaa8, nspecc, omnu, w, wa);
  D=my_cosmo.growthFac(zzz);
  PS=nn*pow(D,2)*my_cosmo.powerSpectrum(KK);

  return PS;

}

// search for scale factor for dndm and halo bias files
double FisherGAL_CLUSTER::findbias_z (double z, double m,string filename, double wm )
{
  ifstream fin;
  double m0,m1, zz, ncon0, nloo0, ncon1, nloo1, mint=pow(10,0.05), intzzz=0.005;
  double m500crit=M200toM500(m,500./wm, z);
  fin.open(filename.c_str());
  if (fin.is_open()) {
    for (int i=0; i< floor(z/intzzz)*70; i++) fin >> m0 >> zz >> nloo0 >> ncon0; // jump to the corresponding z
    for (int i=0; i<70; i++) {
          fin >> m0 >> zz >> nloo0 >> ncon0;
          if ( (m0*mint > m500crit) && ( i!=69) ) {
              fin >> m1 >> zz >> nloo1 >> ncon1;
              return ncon0 + (ncon1-ncon0)/(log10(m1/m0)) * (log10(m500crit/m0));
              //return nloo0 + (nloo1-nloo0)/(log10(m1/m0)) * (log10(m/m0));
          }
    }
  } else {
    cerr << "Error. Cannot open dndm/bias file" << endl;
    exit(0);
  }

  fin.close();
  //return nloo1; //the mass range is beyond the tabulated values,scale according to the last pair
  return ncon1; //the mass range is beyond the tabulated values,scale according to the last pair
}

double FisherGAL_CLUSTER::findbias (double z, double m,string filename, double wm)
{
  double intzzz=0.005;

  if ( z/intzzz-floor(z/intzzz) == 0 ) {
    return findbias_z(z,m,filename, wm);
  } else {
      double z0= floor(z/intzzz)*intzzz;
      double z1= z0+intzzz;
      double y_z0= findbias_z(z0,m,filename, wm);
      double y_z1= findbias_z(z1,m,filename, wm);
      return y_z0 + (z-z1)/(z0-z1) * (y_z1- y_z0);
  }
}

double  FisherGAL_CLUSTER::ComBiasFR(bool frflag, double mm1, double mm2, double wa, double w,double wmmmhh, double wlll, double wbbbhh,double wrrr,  double sigmaa8, double nspecc, double Bm0,double sigma0,double alpha, double beta, double zzz)
{
  double Md=0,Md0=0,ML=0, MassF=0,deltaM=pow(10,0.3), zz, nn, nn2;
  double NNN=0,NNN2=0,xm1,xm2,BIAS;
  double Bm,sigmaMM;
  double omnu = 0.0;
  double h=sqrt(wmmmhh/(1-wlll));
  Cosmology mycosmo(wmmmhh/(h*h)+wrrr, wlll, wbbbhh/(h*h), h, sigmaa8, nspecc, omnu, w, wa);

  Bm=Bm0*pow((1+zzz),alpha);
  sigmaMM=sigma0*pow((1+zzz),beta);


  Md=1e13;
  for(int k=0;k<100;k++){
        if (frflag) {
          nn=findNnn(zzz,Md,FR_name, wmmmhh/(h*h));
          nn2=findbias(zzz,Md,FR_bias, wmmmhh/(h*h));
          if (nn<1) nn=1;
          if (nn2>1) nn2=1;
        }
        else {
          nn=1;
          nn2=1;
        }
        MassF=nn*mycosmo.dndlMTinker(zzz,Md)/Md;
        xm1=(log(mm1)-Bm-log(Md))/(sqrt(2)*sigmaMM);
        xm2=(log(mm2)-Bm-log(Md))/(sqrt(2)*sigmaMM);
        NNN +=MassF*Md*(deltaM-1)*(erfc(xm2)-erfc(xm1)) * nn2*mycosmo.biasTinker(zzz,Md);
        NNN2 +=MassF*Md*(deltaM-1)*(erfc(xm2)-erfc(xm1)) ;
         Md *=deltaM;
  }
    BIAS=NNN/NNN2;
    return BIAS;
}


void FisherGAL_CLUSTER::fisherMatrix_PS()
{
  double **fisherPS, **ifisherPS;
  ifstream fin;
  double NNN, fM, zzz,Veff,Vefffr,Vol,KK=0,KK2=0, ml1, ml2;
  double *dnn;
  CosmParam *c;
  //string CL_name, nameSave;
  //double zmaxx;
  //string file;
  //poi devo inserire il numero di parametri


  c=&fiducialFC;

  nparamFC=12;

  dnn=dvector(1,nparamFC);

  fisherPS=dmatrix(1,nparamFC,1,nparamFC);
  ifisherPS=dmatrix(1,nparamFC,1,nparamFC);


  fin.open(CLL_namePS.c_str());

  for(int i=1;i<=nparamFC;i++){
    for(int j=1;j<=i;j++){
      fM=0.0;
      fin.clear() ;
      fin.seekg(0, ios::beg) ;
      //fin>>KK>>zzz>>NNN>>Vol>>Veff>>Vefffr>>dnn[1]>>dnn[2]>>dnn[3]>>dnn[4]>>dnn[5]>>dnn[6]>>dnn[7]>>dnn[8]>>dnn[9]>>dnn[10];
      fin>>KK>>zzz>>NNN>>Vol>>Veff>>ml1>>ml2>>dnn[1]>>dnn[2]>>dnn[3]>>dnn[4]>>dnn[5]>>dnn[6]>>dnn[7]>>dnn[8]>>dnn[9]>>dnn[10]>>dnn[11]>>dnn[12] ;
      while (!fin.eof()){
	KK2=KK;
	fin>>KK>>zzz>>NNN>>Vol>>Veff>>ml1>>ml2>>dnn[1]>>dnn[2]>>dnn[3]>>dnn[4]>>dnn[5]>>dnn[6]>>dnn[7]>>dnn[8]>>dnn[9]>>dnn[10]>>dnn[11]>>dnn[12] ;
	//fin>>KK>>zzz>>NNN>>Vol>>Veff>>Vefffr>>dnn[1]>>dnn[2]>>dnn[3]>>dnn[4]>>dnn[5]>>dnn[6]>>dnn[7]>>dnn[8]>>dnn[9]>>dnn[10];
	//if(zzz >= 0.15){  //per cambiare l'intervallo in z (attenzione pero' dovrei combiare anche la moltiplicazione per intz) DOESN'T WORK
	  //cout<<NNN<<endl;
          //if (i == 10 && i==10) {
          //  if(KK-KK2>0){fM+=dnn[i]*dnn[j]*Vol*Vefffr*(KK-KK2);}    //following Sartoris et al 2010
          //} else {
        if (ml1 != ml2)
            if(KK-KK2>0){fM+=2*dnn[i]*dnn[j]*Vol*Veff*(KK-KK2);}    //following Sartoris et al 2010
        else
            if(KK-KK2>0){fM+=dnn[i]*dnn[j]*Vol*Veff*(KK-KK2);}    //following Sartoris et al 2010
          //}
	}
      //}
      fisherPS[i][j]=fM/pow(2*PI,2);
      if(i!=j) fisherPS[j][i]=fM/pow(2*PI,2);
      // cout << fisherPS[i][j] << endl;
    }
  }

  /*for(int i=1;i<=nparamFC;i++){
    for(int j=1;j<=i;j++){
      fM=0.0;
      fin.clear() ;
      fin.seekg(0, ios::beg) ;
      //fin>>KK>>zzz>>NNN>>Vol>>Veff>>Vefffr>>dnn[1]>>dnn[2]>>dnn[3]>>dnn[4]>>dnn[5]>>dnn[6]>>dnn[7]>>dnn[8]>>dnn[9]>>dnn[10];
      fin>>KK>>zzz>>NNN>>Vol>>Veff>>Vefffr>>dnn[1]>>dnn[2]>>dnn[3]>>dnn[4]>>dnn[5]>>dnn[6]>>dnn[7]>>dnn[8]>>dnn[9]>>dnn[10]>>dnn[11]>>dnn[12] ;
      while (!fin.eof()){
	KK2=KK;
	fin>>KK>>zzz>>NNN>>Vol>>Veff>>Vefffr>>dnn[1]>>dnn[2]>>dnn[3]>>dnn[4]>>dnn[5]>>dnn[6]>>dnn[7]>>dnn[8]>>dnn[9]>>dnn[10]>>dnn[11]>>dnn[12] ;
      //cout<< "\t" <<KK<< "\t" <<zzz<< "\t" <<NNN<< "\t" <<Vol<< "\t" <<Veff<< "\t" <<Vefffr<< "\t" <<dnn[1]<< "\t" <<dnn[2]<< "\t" <<dnn[3]<< "\t" <<dnn[4]<< "\t" <<dnn[5]<< "\t" <<dnn[6]<< "\t" <<dnn[7]<< "\t" <<dnn[8]<< "\t" <<dnn[9]<< "\t" <<dnn[10]<< "\t" <<dnn[11]<< "\t" <<dnn[12]<< "\t" <<endl ;
	//fin>>KK>>zzz>>NNN>>Vol>>Veff>>Vefffr>>dnn[1]>>dnn[2]>>dnn[3]>>dnn[4]>>dnn[5]>>dnn[6]>>dnn[7]>>dnn[8]>>dnn[9]>>dnn[10];
	//if(zzz/intzz-floor(zzz/intzz)==0){  //per cambiare l'intervallo in z (attenzione pero' dovrei combiare anche la moltiplicazione per intz) DOESN'T WORK
	  //cout<<NNN<<endl;
          //if (i == 10 && i==10) {
          //  if(KK-KK2>0){fM+=dnn[i]*dnn[j]*Vol*Vefffr*(KK-KK2);}    //following Sartoris et al 2010
          //} else {
            if(KK-KK2>0){fM+=dnn[i]*dnn[j]*Vol*Veff*(KK-KK2);}    //following Sartoris et al 2010
          //}
	}
      //}
      fisherPS[i][j]=fM/pow(2*PI,2);
      if(i!=j) fisherPS[j][i]=fM/pow(2*PI,2);
    }
  }*/

  fin.close();
  free_dvector(dnn,1,nparamFC);

  //invertM(fisherPS,nparamFC,ifisherPS);
  invertMatrix(fisherPS,nparamFC,ifisherPS);
  cout<<"Errors:"<<endl;
   cout<<"OmegaM= "<<c->omm<<"\t"<<"OmegaL= "<<c->oml<<"\t"<<"s8= "<<c->sigma8<<"\t"<<"Omegab= "<<c->omb<<"\t"<<"ns= "<<c->nscal<<"\t"<<"Bm= "<<-0.15<<"\t"<<"Sm= "<<0.25<<"\t"<<"alpha= "<<0<<"\t"<<  "beta= "<<0<<"\t"<<"fR= "<<0<<endl;
   //cout<<"OmegaM= "<<c->omm<<"\t"<<"OmegaL= "<<c->oml<<"\t"<<"s8= "<<c->sigma8<<"\t"<<"Omegab= "<<c->omb<<"\t"<<"ns= "<<c->nscal<<"\t"<<"Bm= "<<-0.15<<"\t"<<"Sm= "<<0.25<<"\t"<<"alpha= "<<0<<"\t"<<  "beta= "<<0<<"\t"<<"w0= "<<c->w0<<"\t"<<"fR= "<<0<<endl;
  for(int i=1;i<nparamFC+1;i++){
    cout<<sqrt(ifisherPS[i][i])<<"\t\t";
  }
  cout<<endl;
  writeMatrixS(fisherPS,nparamFC,nameSavePS);
  writeMatrixS(ifisherPS,nparamFC,nameSavePS+"_inverse");

  free_dmatrix(fisherPS,1,nparamFC,1,nparamFC);
  free_dmatrix(ifisherPS,1,nparamFC,1,nparamFC);

}


void FisherGAL_CLUSTER::invertM(double **fisher,int nparam, double **ifisher)
{
  matrix2 <double> fmatrix(nparam,nparam);
  matrix2 <double> ifmatrix(nparam,nparam);
  matrix2 <double> product(nparam,nparam);
  matrix2 <double> f1(5,5);
  matrix2 <double> if1(5,5);
  matrix2 <double> pf(5,5);

  for (int i=0; i< nparam; i++) {
    for (int j=0; j<i; j++) {
      fmatrix.setvalue(i,j, fisher[i+1][j+1]);
      if ( i!=j) fmatrix.setvalue(j,i, fisher[i+1][j+1]);
      if ( i < 5 && j < 5 ) {
        f1.setvalue(i,j,fisher[i+1][j+1]);
        if ( i!=j) f1.setvalue(j,i, fisher[i+1][j+1]);
        }
    }
  }

  ifmatrix.copymatrix2(fmatrix);
  ifmatrix.invert();
  if1.copymatrix2(f1);
  if1.invert();

  // Check
  product.settoproduct(fmatrix, ifmatrix);
  product.comparetoidentity();

  cout << "smaller one" << endl;
  pf.settoproduct(f1, if1);
  pf.comparetoidentity();

  bool xyz;
  double rv;
  for (int i=0; i < ifmatrix.getactualsize(); i++)
    {
    cout << "i=" << i << ": ";
    for (int j=0; j<ifmatrix.getactualsize(); j++)
      {
        fmatrix.getvalue(i,j,rv,xyz);
        cout << rv << " ";
      }
    cout << endl;
    }
  // output
  bool x=true;
  for (int i=0; i<nparam; i++) {
    for (int j=0; j<i; j++) {
      double tmp;
      ifmatrix.getvalue(i,j,tmp,x);
      ifisher[i+1][j+1] = tmp;
      if ( i!=j ) ifisher[j+1][i+1] = tmp;
    }
  }

}
