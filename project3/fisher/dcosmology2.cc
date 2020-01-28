/* dcosmology.cc
 * Includes utility functions for all cosmologies.  Many equations are from
 * Barkana & Loeb (2001), henceforth BL01.
 */

#include <math.h>
#include <iostream>
#include <fstream>
#include "dcosmology2.h"
#include "dnumrecipes.h"
using namespace std;
/**********************************************************************
 ********************* Constructor/Destructor *************************
 *********************************************************************/

Cosmology::Cosmology(double om0, double lam0, double omb, double hparam,
		     double s8, double n, double omNu, double w00,double waa)
{
  omega0 = om0;
  lambda0 = lam0;
  omegab = omb;
  h = hparam;
  omeganu = omNu;
  nspec = n;
  sig8 = s8;
  w0 = w00;
  wa = waa;
  omegam = om0 - omb;
  ombhh = omegab*h*h;
  om0hh = omega0*h*h;
  shape = omega0*h*exp(-omegab-omegab/omega0);
  resetPowerSpectrum();
}

Cosmology::~Cosmology()
{
  //  cout << "Destructor called for Cosmology" << endl;
}

/**********************************************************************
 ******************* General cosmology functions **********************
 *********************************************************************/

/* Returns value of omega at zCurrent for an arbitrary cosmology.
 * BL01 eq. 23 */
double Cosmology::omegaZ(double zCurrent)
{
  double temp, lamtemp, tol(1.0e-5);

  zCurrent += 1.0;
  temp = omega0*pow(zCurrent,3.0);
  if ( wa !=0 ) {
    lamtemp = lambda0*pow(zCurrent,3.0*(1.0+w0+wa)) * exp(-3.0*wa*(zCurrent-1.0)/zCurrent);
  } else {
    lamtemp = lambda0*pow(zCurrent,3*(1+w0));
  }
  //double lamtemp = lambda0;
  temp /= temp + lamtemp + (1 - omega0 - lambda0)*zCurrent*zCurrent;
  //temp /= temp + lambda0 + (1 - omega0 - lambda0)*zCurrent*zCurrent;
  return temp;
}

/* Returns equiv of Omega(z) for the cosmological constant */
double Cosmology::lambdaZ(double zCurrent)
{
  double temp, tol(1.0e-5);

  zCurrent += 1.0;
  //temp = lambda0;
  //cout << w0 << " " << temp << endl;

  if ( wa !=0 ) {
    temp = lambda0*pow(zCurrent,3.0*(1.0+w0+wa)) * exp(-3.0*wa*(zCurrent-1.0)/zCurrent);
  } else {
    temp = lambda0*pow(zCurrent,3*(1+w0));
  }
  temp /= temp + omega0*pow(zCurrent,3.0) + ((1 - omega0 - lambda0)*
					     zCurrent*zCurrent);
  return temp;
}

/* Returns value of Hubble constant at zCurrent in units of H0.  BL eq. 7 */
double Cosmology::hubbleZ(double zCurrent)
{
  double temp,lamZ,tol(1e-5);

  zCurrent += 1.0;

  if ( wa !=0 ) {
    lamZ = lambda0*pow(zCurrent,3.0*(1.0+w0+wa)) * exp(-3.0*wa*(zCurrent-1.0)/zCurrent);
  } else {
    lamZ = lambda0*pow(zCurrent,3*(1+w0));
  }

  temp = ((1.0 - omega0 - lambda0)*zCurrent*zCurrent + lamZ + omega0*pow(zCurrent,3.0));
  return sqrt(temp);
}

/* Returns rhocrit(z) in g/cm^3.  BL eq. 5 */
double Cosmology::rhoCritZ(double zCurrent)
{
  double temp;

  temp = pow(1.0+zCurrent,3.0)*omega0/omegaZ(zCurrent);
  temp *= CRITDENSITY*pow(h,2.0);
  return temp;
}

/* cosmicTime calculates the cosmic time corresponding to a particular
 * value of z, returning H0*time.  It comes from Padmanabhan 2.78 for
 * omega < 1; omega = 1 is standard; omega+lambda from Peebles. */
double Cosmology::cosmicTime(double zCurrent)
{
  double temp,ht;

  if (fabs(omega0-1.0) <1e-5) {
    ht = 2.0/3.0/pow(1.0+zCurrent,1.5);
    return ht;
  }
  if (fabs(omega0 + lambda0-1.0) <1e-5) {          // Peebles 1994 13.20
    temp = sqrt((1.0-omega0)/omega0)/pow(1.0+zCurrent,1.5);
    ht = log(temp+sqrt(1.0+temp*temp));
    ht *= 2.0/3.0/sqrt(1.0-omega0);
    return ht;
  }
  if ((omega0 < 1.0) && (fabs(lambda0) < 1.0e-5)) {
    ht = 2.0*sqrt((1.0-omega0)*(omega0*zCurrent+1.0))/omega0/(1.0+zCurrent);
    temp = (omega0*zCurrent-omega0+2.0)/(omega0*zCurrent+omega0);
    ht -= log(temp+sqrt(temp*temp-1.0));
    ht *= omega0/2.0/sqrt(1.0-omega0+TINY)/(1.0-omega0+TINY);
    return ht;
  }
  //  cout << "Invalid cosmology in cosmicTime!\n";
  return 0.0;
}

/* Given a value of H0*t, calculates corresponding z.  Done analytically for
 * flat cosmologies, but must be inverted numerically for an open cosmology.
 * Uses zbrac and rtsafe to do so.  Lambda cosmology from Peebles 13.20. */
double Cosmology::zFromT(double ht0, double zGuess)
{
  void setOpenInvertTime(double z, double *ht, double *deriv, int flag,
			 double htime, double om0, double lam0);
  void openInvertTime(double z, double *ht, double *deriv, int flag);

  double zReal,temp;

  if (fabs(omega0-1.0)<1.0e-5) {
    zReal = pow(2.0/3.0/ht0,2.0/3.0) - 1.0;
    return(zReal);
  }

  if (fabs(omega0 + lambda0-1.0) < 1.0e-5) {
    temp = sinh(3.0/2.0*ht0*sqrt(1-omega0));
    zReal = pow((1.0-omega0)/omega0,1.0/3.0)/pow(temp,2.0/3.0)-1.0;
    return(zReal);
  }

  if ((omega0 < 1.0) && (fabs(lambda0)<1.0e-5)) {
    zReal = zGuess-1.0;
    setOpenInvertTime(0.0,NULL,NULL,2,ht0,omega0,lambda0);
    zbrac(openInvertTime,zGuess,zReal);
    zReal = rtsafe(openInvertTime,zGuess,zReal,1.0e-5);
    return(zReal);
  }

  //  cerr << "Invalid cosmology sent to zFromT\n";
  return(0.0);
}

/* Function used to calculate z(t) in the root finder */
void setOpenInvertTime(double z, double *ht, double *deriv, int flag,
		       double htime, double om0, double lam0)
{
  static double ht0,omega0,lambda0;
  double x,tiny;

  if (flag==2) {
    omega0 = om0;
    lambda0 = lam0;
    ht0 = htime;
    return;
  }

  tiny = 1.0e-30;
  if ((omega0 < 1.0) && (lambda0 ==0.0)) {
    *ht = 2.0*sqrt((1.0-omega0)*(omega0*z+1.0))/omega0/(1.0+z);
    x = (omega0*z-omega0+2.0)/(omega0*z+omega0);
    *ht -= log(x+sqrt(x*x-1.0));
    *ht *= omega0/2.0/pow(1.0-omega0+tiny,1.5);
    *ht -= ht0;
    if (flag==1) {
      *deriv = sqrt((1.0-omega0)/(1.0+z*omega0))/(1.0+z);
      *deriv -= 2.0*sqrt((1.0-omega0)*(1.0+z*omega0))/omega0/(1.0+z)/(1.0+z);
      *deriv += omega0/(1.0+z*omega0)*sqrt((1.0-omega0)*(1.0+z*omega0)/
				   (omega0+z*omega0)/(omega0+z*omega0));
      *deriv *= omega0/2.0/pow(1.0-omega0+tiny,1.5);
    }
    return;
  }

  //  cerr << "Invalid cosmology sent to openInvertTime\n";
  return;
}

/* Dummy function sent to root finder for z(t) */
void openInvertTime(double z, double *ht, double *deriv,
		    int flag)
{
  setOpenInvertTime(z,ht,deriv,flag,0.0,0.0,0.0);
}

/* Returns derivative dzd(H0t), for a flat or open cosmology.  Latter case
 * actually calculates d(H0t)/dz and returns inverse. */
double Cosmology::dzdh0t(double h0time)
{
  double temp,deriv;

  if (omega0 == 1.0) {
    deriv = -1.0/pow(3.0*h0time/2.0,5.0/3.0);
    return deriv;
  }
  if ((omega0 + lambda0) == 1.0) {
    temp = 3.0/2.0*h0time*sqrt(1.0-omega0);
    deriv = -1.0*sqrt(1.0-omega0)*pow((1-omega0)/omega0,1.0/3.0);
    deriv *= cosh(temp)/pow(sinh(temp),5.0/3.0);
    return deriv;
  }
  if ((omega0 < 1.0) && (lambda0 == 0.0)) { // Use 1/(dt/dz)
    temp = zFromT(h0time);
    deriv = -1.0*sqrt(1.0+omega0*temp)*(1.0+temp)*(1.0+temp);
    return deriv;
  }
  //  cout << "Invalid cosmology in dzdH0t!\n";
  return 0.0;
}

/* Returns dr/dz, in units of c/H0, for given cosmology, at redshift z */
double Cosmology::drdz(double z)
{
  return 1.0/hubbleZ(z);
}

/* Returns coordinate distance for given cosmology, to redshift z,
 * in units of c/H0.  Uses conformal time for the calculation. */
double Cosmology::coordDistance(double zCall)
{
  double cDist;

  if (fabs(omega0-1.0) < 1.0e-5 && fabs(lambda0) < 1.0e-5) {
    // Flat CDM universe
    cDist = 2.0*(1.0 - 1.0/sqrt(1.0 + zCall));
  } else {
    if (fabs(omega0 + lambda0 - 1.0) < 1.0e-5) {
      // Flat lambda-CDM model
      cDist = confTime(0.0) - confTime(zCall);
    } else {
      if (omega0 < 1.0 && fabs(lambda0) < 1.0e-5) {
	// Open matter-dominated model
	cDist = 2.0*(omega0*zCall + (omega0 - 2.0)*
		     (sqrt(1.0+omega0*zCall) - 1.0));
	cDist /= omega0*omega0*(1.0 + zCall);
      } else {
	//	cout << "Invalid cosmology in coordDistance!" << endl;
	cDist=coordDistanceNum(zCall);
	return cDist;
      }
    }
  }
  return cDist;
}

//numerical calculation of  coordDistance
double Cosmology::coordDistanceNum(double zCall)
{
	double cDist;
	double tol(1.0e-5);

	setCoordDistanceKernel(0.0,this,1);
  	cDist=qromb(getCoordDistanceKernel,0.0,zCall,tol);

	return cDist;
}

double setCoordDistanceKernel(double z, Cosmology *c1, int iflag)
{
	static Cosmology *c;
	double kernel;

	if(iflag==1){
		c=c1;
		return 0.0;
	}

	kernel=1.0/c->hubbleZ(z);

	return kernel;
}

double getCoordDistanceKernel(double z)
{
	return setCoordDistanceKernel(z,NULL,0);
}


/* Returns luminosity distance for given cosmology, to redshift z,
 * in units of c/H0.  */
double Cosmology::lumDistance(double zCall)
{
  return (1.0+zCall)*coordDistance(zCall);
}

/* Returns angular diameter distance for given cosmology, to redshift
 * z, in units of c/H0.  */
double Cosmology::angDiamDistance(double zCall)
{
  double omegak,cDist;;
  omegak=1.0-omega0+lambda0;

  cDist=coordDistance(zCall);

  if(omegak>1.0e-5){
	 //open
	 return sinh(sqrt(omegak)*cDist)/(1.0+zCall)/sqrt(omegak);
  }else if(omegak<1.0e-5){
	 //closed
	 return sin(sqrt(-omegak)*cDist)/(1.0+zCall)/sqrt(-omegak);
  }else{
	//flat
	 return cDist/(1.0+zCall);
  }

  cout<<"Error in angDiamDistance"<<endl;
  return coordDistance(zCall)/(1.0+zCall);
}

/* Returns conformal time, in units of 1/H0, to zCall for given
 * cosmology.  */
double Cosmology::confTime(double zCall)
{
  double cTime,betaZ,acosArg;

  if (fabs(omega0-1.0) < 1.0e-5 && fabs(lambda0) < 1.0e-5) {
    // Flat CDM universe
    cTime = 2.0/sqrt(1.0+zCall);
  } else {
    if (fabs(omega0 + lambda0 - 1.0) < 1.0e-5) {
      // Flat lambda-CDM model
      acosArg = ((1.0 + zCall + (1.0 - sqrt(3.0))*
		 pow(1.0/omega0 - 1.0,1.0/3.0))/
		 (1.0 + zCall + (1.0 + sqrt(3.0))*
		  pow(1.0/omega0 - 1.0, 1.0/3.0)));
      betaZ = acos(acosArg);
      cTime = (pow(3.0,-0.25)/pow(omega0*omega0*(1.0-omega0),1.0/6.0)*
	       ellf(betaZ,0.965925826));
    } else {
      if (omega0 < 1.0 && fabs(lambda0) < 1.0e-5) {
	// Open matter-dominated model
	acosArg = 1.0 + 2.0*(1.0-omega0)/omega0/(1.0+zCall);
	betaZ = log(acosArg + sqrt(acosArg*acosArg-1.0));
	cTime = betaZ/sqrt(1.0 - omega0);
      } else {
	//	cout << "Invalid cosmology in confTime!" << endl;
	return 0.0;
      }
    }
  }
  return cTime;
}

/* Hydrogen density as a function of z, assuming helium fraction is 0.24
 * by mass */
double Cosmology::nh(double zCall)
{
  return (8.56e-6*ombhh*pow(1.0+zCall,3.0));
}

// Comoving volume element  dV/dZ d\Omega
// Units:  volume   Mpc^3
double Cosmology::volumeComoving(double zCall)
{
  double cH,r,volume;

  cH=SPEEDOFLIGHT_CGS/getH()/H0_CMSMPC;
  r=cH*coordDistance(zCall);
  volume=cH/hubbleZ(zCall)*r*r;    //comoving volume element

  return volume;
}



//comoving volume out to redshift z
double Cosmology::getVolume(double zCall)
{
	double vol;
	double tol(1.0e-5);

	setDVolumeKernel(0.0,this,1);

	vol=qromb(dummyDVolumeKernel,0.0,zCall,tol);

	return vol;
}

double setDVolumeKernel(double z, Cosmology *c1, int iflag)
{
	static Cosmology *c;

	if(iflag==1){
		c=c1;
		return 0.0;
	}

	return c->volumeComoving(z);  //in Mpc for Lambda CDM
}

double dummyDVolumeKernel(double z)
{
	return setDVolumeKernel(z,NULL,0);
}

/**********************************************************************
 ****************** Perturbation theory functions *********************
 *********************************************************************/

/* Returns fitted overdensity of collapsed halos *relative to
 * critical*, according to chosen cosmology.  Taken from Bryan and
 * Norman (1998); BL01 eq. 22. */
double Cosmology::DeltaC(double zCurrent)
{
  double d;

  d = omegaZ(zCurrent)-1;
  if ((lambda0 == 0.0) && (omega0 < 1.0))
    return (18.0*PI*PI + 60.0*d - 32.0*d*d);
  if (fabs(lambda0 + omega0-1.0) < 1.0e-3)
    return (18.0*PI*PI + 82.0*d - 39.0*d*d);
  //  cout << "Improper cosmology passed through function DeltaC\n";
  return 0.0;
}

/* Fit to the *linear* overdensity at which virialization takes place,
 * it is 1.69 with slight cosmology dependence. */
double Cosmology::delCrit0(double z)
{
  double t1,omZ, omL;

  t1 = 0.15*pow(12.0*PI,2.0/3.0);
  omZ = omegaZ(z);
  omL = lambdaZ(z);
  if (omZ > (1.0 - 1.0e-5)) {
    return t1;
  } else {
    if (omL < 1.0e-5) {
    //if (lambda0 < 1.0e-5) {
      t1 *= pow(omZ,0.0185);
      return t1;
    } else {
      t1 *= pow(omZ,0.0055);
      return t1;
    }
  }
}

/* Returns inverse of linear growth from z to 0:
 * D = D1(z)/D1(z=0); from Eisenstein and Hu.  */
/*double Cosmology::growthFac(double z)
{
  double omZ,lamZ,D;

  omZ = omegaZ(z);
  lamZ = lambdaZ(z);
  D = ((1.0 + zEquality)/(1.0 + z)*5.0*omZ/2.0*
       pow(pow(omZ,4.0/7.0) - lamZ + (1.0+omZ/2.0)*(1.0+lamZ/70.0),
	   -1.0));
  D /= ((1.0 + zEquality)*5.0/2.0*omega0*
	pow(pow(omega0,4.0/7.0) - lambda0 + (1.0+omega0/2.0)*
	    (1.0+lambda0/70.0),-1.0));
  return D;
}*/



double Cosmology::growthFac(double z)
{
  double result=growth(z)/growth(0);
  return result;
  //return growth(z);
}

double Cosmology::growth(double z)
{
  double a=1.0/(1.0+z), acc(1e-6);
  setgrowth(0.0,this, 1);
  double result=qromb(getgrowth,1e-4,a, acc);
  //double G0=1e-3, lnG0=log(G0);
  return hubbleZ(z)*5.0*omega0/2.0*result;
  //return exp(result+lnG0)/G0*a;
    //a*exp(result);
}

double setgrowth(double a, Cosmology *c1, int iflag)
{
  static Cosmology *c;
  if (iflag==1) {
    c= c1;
    return 0.0;
  }

  double zCurrent=1.0/a-1.0;
  double result= 1.0/pow(a*c->hubbleZ(zCurrent),3.0);
  return result;
  //return c->GG(a)/a;
}

double getgrowth(double a)
{
  return setgrowth(a, NULL,0);
}

double Cosmology::GG(double a)
{
  double acc(1e-5);
  double zCurrent=1.0/a-1.0;
  setIntG(0.0,this, 1);
  double result=qromb(getIntG,1e-3,a, acc);
  double lamZ=lambdaZ(zCurrent);
  return -0.5*lamZ - 0.25*pow(a,-2.5)*result;
}

double setIntG(double a, Cosmology *c1, int iflag)
{
  static Cosmology *c;
  if (iflag==1) {
    c= c1;
    return 0.0;
  }
  double zCurrent=1.0/a-1.0;
  double lamZ=c->lambdaZ(zCurrent);
  return pow(a,3.0/2.0)*lamZ;
}

double getIntG(double a)
{
  return setIntG(a,NULL,0);
}

/**********************************************************************
 ********************** n-Sigma Perturbations *************************
 *********************************************************************/

/* Given a redshift z0, finds mass scale corresponding to an n-sigma
 * perturbation in the Press-Schechter mass function.  Uses
 * setMassSigma to perform the inversion.  deltax is the threshold for
 * collapse at z=0; if it is less than zero, use critical overdensity
 * (note that it is an optional parameter).  Uses a root finder on
 * setMassSigma to find the answer. */
double Cosmology::findMassAtSigma(double n, double z0, double deltax)
{
  double mlow,mhigh;
  double xlow,xhigh;
  double ans;

  mlow = 1.0e2;
  mhigh = 1.0e3;
  if (deltax <= 0.0)
    deltax = delCrit0(0.0);
  setMassSigma(0.0,n,z0,deltax,this,1);
  while (1) {
    xlow = massSigma(mlow);
    xhigh = massSigma(mhigh);
    if (xlow > 0.0 && xhigh < 0.0)
      break;
    if (xhigh > 0.0) {
      mlow = mhigh;
      mhigh *= 10.0;
    }
    if (xlow < 0.0) {
      mhigh = mlow;
      mlow /= 10.0;
    }
  }
  //  cout << "mlow=" << mlow << " mhigh=" << mhigh << endl;
  ans = zriddrSimp(massSigma,mlow,mhigh,1.0e-4);
  return ans;
}

/* Function passed to root finder to compute n-sigma perturbations.
 *   n = number of sigmas
 *   deltax = z=0 collapse threshold */
double setMassSigma(double mass, double n, double z, double deltax,
		   Cosmology *c1, int flag)
{
  static double N,z0,dx;
  static Cosmology *c;
  double empty,ans;

  if (flag == 1) {
    N = n;
    z0 = z;
    c = c1;
    dx = deltax;
    return 0.0;
  }

  ans = c->sigma0fM(mass,empty,0);
  //  cout << "mass = " << mass << " ans=" << ans << endl;
  ans += -1.0/N*dx*(c->delCrit0(z0)/c->delCrit0(0.0))/c->growthFac(z0);
  //  ans += -1.0/N*dx;
  //  cout << "mass = " << mass << " ans=" << ans << endl;
  return ans;
}

/* Dummy function for finding n-sigma perturbations */
double massSigma(double mass)
{
  return setMassSigma(mass,0.0,0.0,0.0,NULL,0);
}

/**********************************************************************
 ************************ Collapse Fraction ***************************
 *********************************************************************/

/* Collapse fraction above mMin for Press-Schechter mass function.  Note
 * mMin is an optional parameter; if not input it is taken to be the mass
 * above halos with 10^4 K (the minimum cooling mass) */
double Cosmology::fCollPSExact(double z, double mMin)
{
  double sigma,dCritZ,dummy;
  double ans;

  if (mMin < 0.0)
    mMin = coolMass(this,z);
  sigma = sigma0fM(mMin,dummy,0);
  dCritZ = delCrit0(z)/growthFac(z);
  ans = erfc(dCritZ/sqrt(2.0)/sigma);
  return ans;
}

/* Calculates fraction of mass in galaxies by directly integrating the
 * mass function (using setFCollInt).
 *   z = redshift
 *   mMin = minimum mass; if <0 use Tvir=10^4 K (optional)
 *   massFcn = which mass function to use; 0=PS,1=ST,2=Jenkins (optional)
 */
double Cosmology::fColl(double z, double mMin, int massFcn)
{
  double mMax;
  double ans,oldAns;

  setFCollInt(0.0,this,z,massFcn,1);
  ans = 0.0;
  if (mMin < 0.0)
    mMin = coolMass(this,z);
  while (1) {
    mMax = mMin*10.0;
    oldAns = ans;
    ans += qromb(fCollInt,mMin,mMax,1.0e-4);
    if (fabs(ans-oldAns)/ans < 1.0e-4)
      break;
    //    cout << "mMin=" << mMin << " mMax=" << mMax << " fColl=" << ans
    //	 << endl;
    mMin = mMax;
  }
  ans /= CRITDENMSOLMPC*h*h*omega0;
  return ans;
}

/* Integrand for collapse fraction.  fcn is massFcn, says which mass function
 * to use. */
double setFCollInt(double mass, Cosmology *c1, double z1, int fcn,
		  int flag)
{
  static Cosmology *c;
  static double z;
  static int massFcn;
  double result;

  if (flag == 1) {
    c = c1;
    z = z1;
    massFcn = fcn;
    return 0.0;
  }

  if (massFcn == 1)
    result = c->dndlMSheth(z,mass);
  else if (massFcn == 2)
    result = c->dndlMJenkins(z,mass);
  else
    result = c->dndlM(z,mass);
  return result;
}

/* Dummy integrand for collapse fraction */
double fCollInt(double mass)
{
  return setFCollInt(mass,NULL,0.0,0,0);
}

/* Calculates number of objects above mMin.  If mMin < 0.0, assumes it is
 * minimum cooling mass (this is an optional parameter). Only set up
 * for PS mass function. */
double Cosmology::nCollObject(double z, double mMin)
{
  double mMax;
  double ans,oldAns;

  setNCollInt(0.0,this,z,1);
  ans = 0.0;
  if (mMin < 0.0)
    mMin = coolMass(this,z);
  while (1) {
    mMax = mMin*10.0;
    oldAns = ans;
    ans += qromb(nCollInt,mMin,mMax,1.0e-4);
    if (fabs(ans-oldAns)/ans < 1.0e-4)
      break;
    //    cout << "mMin=" << mMin << " mMax=" << mMax << " fColl=" << ans
    //	 << endl;
    mMin = mMax;
  }
  return ans;
}

/* Integrand for finding number of collapsed objects above threshold */
double setNCollInt(double mass, Cosmology *c1, double z1, int flag)
{
  static Cosmology *c;
  static double z;

  if (flag == 1) {
    c = c1;
    z = z1;
    return 0.0;
  }

  return c->dndlM(z,mass)/mass;
}

/* Dummy integrand for finding number of collapsed objects above threshold */
double nCollInt(double mass)
{
  return setNCollInt(mass,NULL,0.0,0);
}

/**********************************************************************
 ************************* Bias Calculations **************************
 *********************************************************************/

/* Linear bias as a function of halo mass for Press-Shechter halos.
 * See, e.g., Cooray & Sheth (2002), s. 3.4. */
double Cosmology::biasPS(double z, double mass)
{
  double deltasc,deltasc0,sig;
  double bias;
  double dummy;

  deltasc = delCrit0(z)/growthFac(z);
  deltasc0 = delCrit0(0.0);
  sig = sigma0fM(mass,dummy,0);
  bias = 1.0 + (deltasc*deltasc/sig/sig - 1.0)/deltasc0;
  return bias;
}

double Cosmology::biasTinker(double z, double mass)
{
  double deltasc,deltasc0,sig, nu;
  double bias;
  double dummy;

  double delta(200), y=log10(delta);
  double a=0.44*y-0.88, A=1.0+0.24*y*exp(-pow(4.0/y,4.0)), b(1.5), B(0.183), cc(2.4), C=0.19+0.107*y+0.19*exp(-pow(4.0/y,4.0));
  deltasc = delCrit0(z)/growthFac(z);
  deltasc0 = delCrit0(0.0);
  sig = sigma0fM(mass,dummy,0);
  nu=deltasc/sig;
  bias = 1.0 - A*pow(nu,a)/(pow(nu,a) + pow(deltasc,a)) + B*pow(nu,b) + C*pow(nu,cc);
  return bias;
}
/**********************************************************************
 ******************** Press-Schechter functions ***********************
 *********************************************************************/
/* Routines to calculate Press-Schecter mass function.
 * Transfer function from Eisenstein and Hu (1998); does not include
 * possibility of hot dark matter or neutrinos (despite the parameters!).
 * Translated from Fortran code provided by Rennan Barkana, transfer
 * function comes ultimately from Wayne Hu's website.
 */

/* Calculates Press-Schechter mass function
 * tM = halo mass (Msun)
 */
double Cosmology::dndlM(double z, double tM)
{
  double dCritZ,sigM,dlsdlM,tdn,dsdM;

  dCritZ = delCrit0(z)/growthFac(z);
  sigM = sigma0fM(tM,dsdM,1);
  dlsdlM = tM*dsdM/sigM;
  tdn = (sqrt(2.0/PI)*dCritZ*fabs(dlsdlM)*
	 exp(-dCritZ*dCritZ/(2.0*sigM*sigM))/(tM*sigM));
  tdn *= CRITDENMSOLMPC*om0hh;
  return tdn;
}

/* Calculates Sheth-Tormen (1999) mass function
 *   tM = halo mass (Msun)
 */
double Cosmology::dndlMSheth(double z, double tM)
{
  const double a(0.707),A(0.322),p(0.3);
  double dCritZ,sigM,dlsdlM,tdn,dsdM,nuc;

  dCritZ = delCrit0(z)/growthFac(z);
  sigM = sigma0fM(tM,dsdM,1);
  dlsdlM = tM*dsdM/sigM;
  nuc = dCritZ/sigM;
  tdn = A*sqrt(2.0*a/PI)*fabs(dlsdlM)/tM*nuc*(1.0+pow(nuc*nuc*a,-p));
  tdn *= exp(-a*nuc*nuc/2.0);
  tdn *= CRITDENMSOLMPC*om0hh;
  return tdn;
}

/* Calculates Jenkins et al. (2001)
 *   tM = halo mass (Msun)
 */
double Cosmology::dndlMJenkins(double z, double tM)
{
  const double a(0.73),A(0.353),p(0.175);
  double dCritZ,sigM,dlsdlM,tdn,dsdM,nuc;

  dCritZ = delCrit0(z)/growthFac(z);
  sigM = sigma0fM(tM,dsdM,1)*growthFac(z);
  dlsdlM = growthFac(z)*tM*dsdM/sigM;
  nuc = dCritZ/sigM;
  //cout<<"tM ="<<tM<<"\t"<<"dCritZ ="<<dCritZ<<"\t"<<"sigM ="<<sigM<<"\t"<<"dlsdlM ="<<dlsdlM<<endl;
  tdn=-(0.316*exp(-1.*pow(fabs(log(sigM)-0.67), 3.82)))*CRITDENMSOLMPC*om0hh*(dlsdlM)/tM;
  //cout << -1.*pow(fabs(log(sigM)-0.67), 3.82)  << endl;

  //tdn = A*sqrt(2.0*a/PI)*fabs(dlsdlM)/tM*nuc*(1.0+pow(nuc*nuc*a,-p));
  //tdn *= exp(-a*nuc*nuc/2.0);
  //tdn *= CRITDENMSOLMPC*om0hh;
  return tdn;
}

/* Calculates Tinker et al. (2010)
 *   tM = halo mass (Msun)
 */
double Cosmology::dndlMTinker(double z, double tM)
{
  const double delta(200);
  double a0(1.47),A0(0.186),b0(2.57), c(1.19);
  //double A0=0.1*log10(delta)-0.05, a0=1.43+pow(log10(delta)-2.3,1.5), b0=1.0+pow(log10(delta)-1.6,-1.5), c=1.2+pow(log10(delta)-2.35, 1.6);
  double tt=pow(0.75/log10(delta/75), 1.2), alpha=pow(10,-tt);
  double A=A0*pow(1+z,-0.14), a=a0*pow(1+z,-0.06), b=b0*pow(1+z,-alpha);
  double dCritZ,sigM,dlsdlM,tdn,dsdM,nuc;

  dCritZ = delCrit0(z)/growthFac(z);
  sigM = sigma0fM(tM,dsdM,1)*growthFac(z);
  dlsdlM = growthFac(z)*tM*dsdM/sigM;
  //tdn=-A*(1+ pow(sigM/b,-a))*exp(-c/(sigM*sigM))  *(dlsdlM);
  tdn=-A*(1+ pow(sigM/b,-a))*exp(-c/(sigM*sigM))*CRITDENMSOLMPC*om0hh*(dlsdlM)/tM;

  return tdn;
}



/* Calculates and returns sigma(tM), and if iDeriv=1, calculates and stores
 * dsdM = d sigma(tM)/d tM.  Procedure comes from Rennan Barkana and is
 * rather complex; could probably be improved.
 *   tM = halo mass (Msun)
 *   dsdM = place to store derivative (if iDeriv=0, ignored)
 */
double Cosmology::sigma0fM(double tM, double &dsdM, int iDeriv)
{
  double klimits[MAXBINS];
  double acc,klim,khi,klhi,klow,kfac;
  double sigma0fm,dsigma0fm;
  int nk;

  // scale in comoving Mpc
  scale = pow(tM*3.0/(4.0*PI*2.77545e11*om0hh),1.0/3.0);

  setSigmatop(0.0,this,1);
  setdsigmatop(0.0,this,1);
  acc = 2.0e-6;
  klimits[0] = fmin(1.0e-3/scale,0.2);
  klim = 2.0e2/scale;
  khi = klim;
  klhi = log(khi);
  if (klhi > klimits[0])
    klim = fmin(klim,khi);
  nk = int(log(klim/klimits[0])/log(10.0)+0.5);
  if (nk < 5) {
    nk = 5;
  } else {
    if (nk > 7) {
      nk = 7;
    }
  }
  //cout << "klimits[0]= " << klimits[0] << endl;
  kfac = log(klim/klimits[0])/(nk - 1.0);
  klimits[0] = log(klimits[0]);
  for (int i=1;i<nk;i++) {
    klimits[i] = klimits[i-1]+kfac;
  }
  klow = 1.0e-9;
  if (khi < klow) {
    sigma0fm = 0.0;
  } else {
    if (klhi < klimits[0]) {
      sigma0fm = qromb(sigmatop2,klow,khi,acc);
      //cout << "klimits1= " << klow << "\t" << khi << endl;
        //      cout << "In low sigma0fm =" << sigma0fm;
    } else {
      sigma0fm = qromb(sigmatop2,klow,exp(klimits[0]),acc);
      //cout << "klimits2= " << klow << "\t" << exp(klimits[0]) << endl;
      //      cout << "In high, sigma0fm = " << sigma0fm << endl;
      for (int j=1;j<nk;j++) {
	sigma0fm += qromb(sigmatop,klimits[j-1],klimits[j],acc);
      //cout <<  " klimits3= " << exp(klimits[j-1]) << "\t" << exp(klimits[j]) << endl;
	//	cout << "In high, sigma0fm = " << sigma0fm << endl;
      }
    }
  }

  sigma0fm = sNorm*sqrt(sigma0fm);
  //  cout << "Final sigm0fm = " << sigma0fm << endl;

  if (iDeriv == 1) {
    if (khi < klow) {
      dsigma0fm = 0.0;
    } else {
      if (klhi < klimits[0]) {
	dsigma0fm = qromb(dsigmatop2,klow,khi,acc);
      } else {
	dsigma0fm = qromb(dsigmatop2,klow,exp(klimits[0]),acc);
	for (int jj=1;jj<nk;jj++) {
	  dsigma0fm += qromb(dsigmatop,klimits[jj-1],klimits[jj],acc);
	}
      }
    }
    dsigma0fm += sigmatop(klhi)/scale;
    dsdM = 1.0/(PI*8.0*2.77545e11*om0hh*scale*scale);
    dsdM = sNorm*sNorm*dsdM*dsigma0fm/sigma0fm;
  }

  return sigma0fm;
}

/* Calculates sigmatophat, assuming k is entered as linear. */
double sigmatop2(double k)
{
  return sigmatop(log(k))/k;
}

/* Calculates derivative of sigmatophat, assuming k is entered as linear. */
double dsigmatop2(double k)
{
  return dsigmatop(log(k))/k;
}

/* Calculates sigmatophat, assuming k is entered as ln(k) */
double setSigmatop(double kl, Cosmology *c1, int flag)
{
  static Cosmology *c;
  double k,x,sigtop;
  double dummy;

  if (flag == 1) {
    c = c1;
    return 0.0;
  }

  k = exp(kl);
  x = c->getScale()*k;

  sigtop = pow(k,3.0+c->getnSpec())*pow(c->TFMaster(k,dummy,0),2.0)*
           pow(3.0*(x*cos(x) - sin(x))/pow(x,3.0),2.0);
  //  sigtop = pow(k,3.0+c->getnSpec())*pow(c->TF_BBKS(k),2.0)*
  //           pow(3.0*(x*cos(x) - sin(x))/pow(x,3.0),2.0);
  //  cout << "kl=" << kl << " sigtop = " << sigtop << endl;
  //  cout << "master= " << c->TFMaster(k,dummy,0) << " k=" << k
  //       << " x=" << x << endl;
  return sigtop;
}

double sigmatop(double kl)
{
  return setSigmatop(kl,NULL,0);
}

/* Calculates derivative of sigmatophat, assuming k is entered as log(k).
 * Derivative calculation; used integration by parts to get a positive
 * integrand */
double setdsigmatop(double kl, Cosmology *c1, int flag)
{
  static Cosmology *c;
  double k,x,t1,dTFdk,dsigmatop;

  if (flag == 1) {
    c = c1;
    return 0.0;
  }

  k = exp(kl);
  x = c->getScale()*k;

  t1 = c->TFMaster(k,dTFdk,1);
  dsigmatop = -(pow(k,3.0+c->getnSpec())*t1*
		(2.0*dTFdk*k+t1*(3.0+c->getnSpec()))*
		pow(3.0*(x*cos(x) - sin(x))/pow(x,3.0),2.0)
		/(c->getScale()));
  return dsigmatop;
}

double dsigmatop(double kl)
{
  return setdsigmatop(kl,NULL,0);
}

/************************************************************************
 *********************** Power Spectrum Functions ***********************
 ***********************************************************************/

/* Return power spectrum, assuming the EH transfer function
 *   k = wavenumber (comoving Mpc^-1) */
double Cosmology::powerSpectrum(double k)
{
  double temp;
  double dummy;

  temp = 2.0*pow(PI*sNorm*TFMaster(k,dummy,0),2.0);
  //  temp = 2.0*pow(PI*sNorm*TF_BBKS(k),2.0);
  temp *= pow(k,nspec);
  return temp;
}

/* Calculates constants necessary for the calculation of the power spectrum
 * and sigma */
void Cosmology::resetPowerSpectrum()
{
  double acc,sigma8;

  thetaCMB = TCMB/2.7;
  // Set up transfer function parameters
  TFSetParameters();
  // Normalize sigma8 at present time
  scale = 8.0/h;
  acc = 1.0e-6;
  setSigmatop(0.0,this,1);
  sigma8 = qromb(sigmatop2,1.0e-6,0.001/scale,acc);
  sigma8 += qromb(sigmatop,log(0.001/scale),log(0.1/scale),acc);
  sigma8 += qromb(sigmatop,log(0.1/scale),log(1.0/scale),acc);
  sigma8 += qromb(sigmatop,log(1.0/scale),log(10.0/scale),acc);
  sigma8 += qromb(sigmatop,log(10.0/scale),log(100.0/scale),acc);
  sigma8 = sqrt(sigma8);
  sNorm = sig8/sigma8;
}

/* Sets up transfer function parameters.  Ultimately from Hu and
 * Eisenstein (1998) online code
 * (http://www.sns.ias.edu/~whu/transfer/transfer.html). */
void Cosmology::TFSetParameters()
{
  double pC,pCB,fC,fCB,fNUB,yD,RDrag,REquality,zDrag,kEquality;
  double fNu,fBaryon;

  fNu = omeganu/omega0;
  fBaryon = omegab/omega0;

  zEquality = 2.50e4*om0hh*pow(thetaCMB,-4.0) - 1.0;
  kEquality = 0.0746*om0hh*pow(thetaCMB,-2.0);
  zDrag = 0.313*pow(om0hh,-0.419)*(1.0+0.607*pow(om0hh,0.674));
  zDrag = 1.0+zDrag*pow(ombhh,0.238*pow(om0hh,0.223));
  zDrag *= 1291.0*pow(om0hh,0.251)/(1.0 + 0.659*pow(om0hh,0.828));

  yD = (1.0 + zEquality)/(1.0+zDrag);
  RDrag = 31.5*ombhh*pow(thetaCMB,-4.0)*1000.0/(1.0+zDrag);
  REquality = 31.5*ombhh*pow(thetaCMB,-4.0)*1000.0/(1.0+zEquality);
  soundHorizon = 2.0/3.0/kEquality*sqrt(6.0/REquality)*
                 log((sqrt(1.0+RDrag) + sqrt(RDrag+REquality))/
		     (1.0+sqrt(REquality)));
  /* Fit from Hu and Eisenstein */
  soundHorizon = 44.5*log(9.83/om0hh)/sqrt(1.0+10.0*pow(ombhh,0.75));

  fC = 1.0 - fNu - fBaryon;
  fCB = 1.0 - fNu;
  pC = -(5.0-sqrt(1.0+24.0*fC))/4.0;
  pCB = -(5.0-sqrt(1.0+24.0*fCB))/4.0;
  fNUB = fNu + fBaryon;

  alphaNu = (fC/fCB)*(2.0*(pC+pCB)+5.0)/(4.0*pCB+5.0);
  alphaNu *= (1.0-0.553*fNUB+0.126*pow(fNUB,3.0));
  alphaNu /= 1.0 - 0.193*sqrt(fNu) + 0.169*fNu;
  alphaNu *= pow(1.0+yD,pC-pCB);
  alphaNu *= (1.0+(pCB-pC)/2.0*(1.0+1.0/(4.0*pC+3.0)/(4.0*pCB+7.0))/
             (1.0+yD));
  betaC = 1.0/(1.0 - 0.949*fNUB);
}

/* Calculates and returns TF(k), and if iDeriv=1, calculates and stores
 * dTF/dk.  Transfer function is from Hu and Eisenstein online code:
 * http://www.sns.ias.edu/~whu/transfer/transfer.html
 * Derivative calculation by Rennan Barkana */
double Cosmology::TFMaster(double k, double &dTFdk, int iDeriv)
{
  double q,gammaEff,qEff,tfMaster,t2,v,dv,dgamedk;
  double t1(0.0),dt1dqe(0.0);

  q = k*thetaCMB*thetaCMB/om0hh;
  gammaEff = (sqrt(alphaNu) +
	      (1.0 - sqrt(alphaNu))/
	      (1.0 + pow(0.43*k*soundHorizon,4.0)));
  qEff = q/gammaEff;
  /*
  cout << "\nk=" << k << "gammaEff=" << gammaEff << " qEff=" << qEff
       << endl;
  */

  tfMaster = log(exp(1.0)+1.84*betaC*sqrt(alphaNu)*qEff);
  if (iDeriv == 1) {
    t1 = tfMaster;
    dt1dqe = 1.84*betaC*sqrt(alphaNu)/exp(tfMaster);
  }
  tfMaster = tfMaster/(tfMaster + qEff*qEff*
		       (14.4+325.0/(1.0+60.5*pow(qEff,1.11))));
  if (iDeriv ==1) {
    t2 = tfMaster;
    v = t1/t2;
    dv = dt1dqe + 2.0*qEff*(14.4+325.0/(1.0+60.5*pow(qEff,1.11))) -
         325.0*60.5*1.11*pow(qEff,2.11)/pow(1.0+60.5*pow(qEff,1.11),2.0);
    dgamedk = -(1.0-sqrt(alphaNu))*4.0*pow(k,3.0)*
              pow(0.43*soundHorizon,4.0)/
              pow(1.0+pow(0.43*soundHorizon*k,4.0),2.0);
    dTFdk = thetaCMB*thetaCMB/om0hh*((v*dt1dqe-t1*dv)/v/v)*
            (1.0/gammaEff-k*dgamedk/gammaEff/gammaEff);
  }
  return tfMaster;
}

/* An alternative, older transfer function from BBKS.  Described in
 * Peacock's book */
double Cosmology::TF_BBKS(double k)
{
  double qBBKS,gammaEff;
  double tf;

  gammaEff = om0hh*exp(-omegab*(1.0+sqrt(2.0*h)/omega0));
  qBBKS = k*thetaCMB*thetaCMB/gammaEff;

  /*  For Adiabatic CDM perturbations, Peacock 15.82 */
  tf = pow(1.0 + 3.89*qBBKS + pow(16.1*qBBKS,2.0) + pow(5.46*qBBKS,3.0) +
	   pow(6.71*qBBKS,4.0),-0.25);
  tf *= log(1.0 + 2.34*qBBKS)/2.34/qBBKS;
  return tf;
}

/************************************************************************
 ****************** Extended Press-Schechter Functions ******************
 ***********************************************************************/

/* Returns rate at which a halo of mass mOrig merges with halos of
 * mass mAccrete at redshift z, per logarithmic interval in mAccrete and
 * per redshift.  Uses extended Press-Schechter formalism.  Taken from
 * Lacey & Cole (1993), eq. 2.18. */
double Cosmology::mergerRate(double mTot, double mOrig, double z)
{
  double mAccrete;
  double dCritZ,sigMTot,sigMOrig,tdn,dsdMTot(0.0);
  double temp(0.0);
  double ddeltacdz;

  dCritZ = delCrit0(z)/growthFac(z);
  ddeltacdz = fabs((delCrit0(z-0.01)/growthFac(z-0.01)-
		    delCrit0(z+0.01)/growthFac(z+0.01))/0.02);
  //  cout << "dCritZ=" << dCritZ << " ddeltacdz=" << ddeltacdz << endl;
  mAccrete = mTot - mOrig;
  sigMTot = sigma0fM(mTot,dsdMTot,1);
  sigMOrig = sigma0fM(mOrig,temp,0);
  //  cout << "sigMTot=" << sigMTot << " sigMOrig=" << sigMOrig << endl;
  //  cout << "dsdMTot=" << dsdMTot << endl;
  temp = 1.0 - pow(sigMTot/sigMOrig,2.0);
  tdn = (sqrt(2.0/PI)*fabs(ddeltacdz)*mAccrete*fabs(dsdMTot)/
	 sigMTot/sigMTot/pow(temp,1.5)*
	 exp(-dCritZ*dCritZ/2.0*(1.0/sigMTot/sigMTot-
				 1.0/sigMOrig/sigMOrig)));
  //  cout << "mTot=" << mTot << " mAcc=" << mAccrete
  //       << " mOrig=" << mOrig << " mergerRate=" << tdn << endl;
  return tdn;
}

/* Returns symmetrized merger kernel.  See Benson et al. (2004), eq. (5). */
double Cosmology::mergerKernelSym(double m1, double m2, double z)
{
  return ((mergerKernel(m1,m2,z) + mergerKernel(m2,m1,z))/2.0);
}

/* Returns merger kernel.  See Benson et al. (2004), eq. (5). */
double Cosmology::mergerKernel(double m1, double m2, double z)
{
  double mTot;
  double dCritZ,sigMTot,sigM1,sigM2;
  double dsdM2(0.0),dsdMTot(0.0);
  double temp(0.0);
  double ddeltacdz;
  double ans;

  dCritZ = delCrit0(z)/growthFac(z);
  ddeltacdz = fabs((delCrit0(z-0.01)/growthFac(z-0.01)-
		    delCrit0(z+0.01)/growthFac(z+0.01))/0.02);
  //  cout << "dCritZ=" << dCritZ << " ddeltacdz=" << ddeltacdz << endl;

  mTot = m1 + m2;
  sigMTot = sigma0fM(mTot,dsdMTot,1);
  sigM1 = sigma0fM(m1,temp,0);
  sigM2 = sigma0fM(m2,dsdM2,1);
  /* Convert to logarithmic derivatives */
  dsdMTot = mTot/sigMTot*dsdMTot;
  dsdM2 = m2/sigM2*dsdM2;
  temp = 1.0 - pow(sigMTot/sigM1,2.0);

  ans = fabs(ddeltacdz/dCritZ)*fabs(dsdMTot)/fabs(dsdM2)/pow(temp,1.5);
  ans *= sigM2/sigMTot*m2*m2/mTot;
  ans *= exp(-dCritZ*dCritZ/2.0*
	     (1.0/sigMTot/sigMTot-1.0/sigM1/sigM1-1.0/sigM2/sigM2));
  ans *= 1.0/CRITDENMSOLMPC/om0hh;
  return ans;
}

/* Evaluates d deltac(z)/d t at redshift z; units are s^-1 */
double Cosmology::deltacderiv(double z)
{
  double step,err(0.0),ans, lamZ,tol(1e-5);

  setDeltacritDeriv(0.0,this,1);
  step = 0.1;
  ans = dfridr(deltaCritDeriv,z,step,&err);
  if (fabs(err/ans) > 1.0e-3*ans) {
    //    cout << "Warning! Possible error in deltacderiv!" << endl;
    //    cout << "z= " << z << " deriv=" << ans << endl;
  }
  //ans *= -h*UNH*(1.0+z)*sqrt(omega0*pow(1.0+z,3.0)+lambda0);

  if ( wa !=0 ) {
    lamZ = lambda0*pow(1+z,3*(1+w0+wa)) * exp(-3.0*wa*z/(1+z));
  } else {
    lamZ = lambda0*pow(1+z,3*(1+w0));
  }

  ans *= -h*UNH*(1.0+z)*sqrt(omega0*pow(1.0+z,3.0)+lamZ);
  return ans;
}

double setDeltacritDeriv(double z, Cosmology *cos, int flag)
{
  static Cosmology *c;
  double ans;

  if (flag == 1) {
    c = cos;
    return 0.0;
  }

  ans = (c->delCrit0(z)/c->growthFac(z));
  return ans;
}

double deltaCritDeriv(double z)
{
  return setDeltacritDeriv(z,NULL,0);
}

/* Evaluates P(zformation < zform), where zformation is the redshift
 * at which largest parent has half the mass of the given halo.
 * Uses (2.28) in Lacey & Cole (1993). */
double Cosmology::formTimeCDF(double zform, double mass, double zfinal)
{
  double omegaform,omegafin;
  double omtilde,temp(0.0);
  double sh,s2;
  double *ans;
  double result,step;
  int goodSteps,badSteps;

  ans = dvector(1,1);
  omegaform = delCrit0(zform)/growthFac(zform);
  omegafin = delCrit0(zfinal)/growthFac(zfinal);
  sh = pow(sigma0fM(mass/2.0,temp,0),2.0);
  s2 = pow(sigma0fM(mass,temp,0),2.0);
  omtilde = (omegaform-omegafin)/sqrt(sh-s2);
  //  cout << "omform=" << omegaform << " omfin=" << omegafin
  //       << "\nsh=" << sh << " s2=" << s2
  //       << "\nomtilde=" << omtilde << endl;
  setFormTimeCDFInt(0.0,ans,ans,this,omtilde,mass,sh,s2,1);
  step = omtilde*omtilde/3.0;
  ans[1] = 0.0;
  odeint(ans,1,1.0,1.0e-10,1.0e-6,step,0.0,&goodSteps,
	 &badSteps,formTimeCDFInt,bsstep);
  // Absolute value is because we integrate backwards in previous step
  result = fabs(ans[1]);
  free_dvector(ans,1,1);
  return result;
}

void setFormTimeCDFInt(double stilde, double y[], double result[],
		       Cosmology *cos, double omt, double m,
		       double shalf, double stwo, int flag)
{
  static Cosmology *c;
  static double omtilde,mfin;
  static double s2,sh;
  double mass;
  double s;

  if (flag == 1) {
    c = cos;
    omtilde = omt;
    mfin = m;
    sh = shalf;
    s2 = stwo;
    return;
  }

  s = stilde*(sh-s2)+s2;
  setInvertSigmaM(0.0,c,sqrt(s),1);
  mass = zriddrSimp(invertSigmaM,mfin/2.0,mfin,1.0e-6);
  //  cout << "m2=" << mfin << " m1=" << mass << endl;
  result[1] = mfin/mass/sqrt(2.0*PI)*omtilde/pow(stilde,1.5);
  result[1] *= exp(-omtilde*omtilde/2.0/stilde);
  //  cout << stilde << "\t" << result[1] << endl;
}

void formTimeCDFInt(double s, double y[], double result[])
{
  setFormTimeCDFInt(s,y,result,NULL,0.0,0.0,0.0,0.0,0);
}

double setInvertSigmaM(double mass, Cosmology *cos, double sig, int flag)
{
  static Cosmology *c;
  static double s;
  double temp(0.0);

  if (flag == 1) {
    s = sig;
    c = cos;
    return 0.0;
  }

  return (c->sigma0fM(mass,temp,0)-s);
}

double invertSigmaM(double mass)
{
  return setInvertSigmaM(mass,NULL,0.0,0);
}

/* Given a mass and a redshift, returns z for which half of the halos
 * had a parent with at least half the final mass, using the Lacey &
 * Cole (1993) analytic method. */
double Cosmology::medianFormTime(double mass, double zfinal)
{
  double ans,temp,zmin,zmax;

  setFindMedian(0.0,mass,zfinal,this,1);
  zmin = zfinal+0.1;
  temp = formTimeCDF(zmin,mass,zfinal);
  if (temp < 0.5)
   zmin = zfinal+0.01;
  zmax = zmin+1.0;
  while (1) {
    temp = formTimeCDF(zmax,mass,zfinal);
    if (temp < 0.5)
      break;
    zmax++;
  }
  ans = zriddrSimp(findMedian,zmin,zmax,1.0e-4);
  return ans;
}

double setFindMedian(double zform, double m, double zfin,
		     Cosmology *cos, int flag)
{
  static Cosmology *c;
  static double mass,zfinal;
  double ans;

  if (flag == 1) {
    c = cos;
    mass = m;
    zfinal = zfin;
    return 0.0;
  }

  ans = (c->formTimeCDF(zform,mass,zfinal)-0.5);
  //  cout << "zform=" << zform << " ans=" << ans << endl;
  return ans;
}

double findMedian(double zform)
{
  return setFindMedian(zform,0.0,0.0,NULL,0);
}

/* Evaluates dp/dtform at zform:  probability density that a halo parent
 * first has m > mass/2 at time corresponding to zform. */
double Cosmology::formTimePDF(double zform, double mass, double zfinal)
{
  double omegaform,omegafin;
  double omtilde,temp(0.0);
  double sh,s2;
  double *ans;
  double result,step;
  int goodSteps,badSteps;

  ans = dvector(1,1);
  omegaform = delCrit0(zform)/growthFac(zform);
  omegafin = delCrit0(zfinal)/growthFac(zfinal);
  sh = pow(sigma0fM(mass/2.0,temp,0),2.0);
  s2 = pow(sigma0fM(mass,temp,0),2.0);
  omtilde = (omegaform-omegafin)/sqrt(sh-s2);
  //  cout << "omform=" << omegaform << " omfin=" << omegafin
  //       << "\nsh=" << sh << " s2=" << s2
  //       << "\nomtilde=" << omtilde << endl;
  setFormTimePDFInt(0.0,ans,ans,this,omtilde,mass,sh,s2,1);
  step = omtilde*omtilde/3.0;
  ans[1] = 0.0;
  odeint(ans,1,1.0,1.0e-10,1.0e-6,step,0.0,&goodSteps,
	 &badSteps,formTimePDFInt,bsstep);
  result = ans[1];
  free_dvector(ans,1,1);
  // Convert from dp/domegatilde to dp/dt
  result *= 1/sqrt(sh-s2)*deltacderiv(zform);
  return result;
}

void setFormTimePDFInt(double stilde, double y[], double result[],
		       Cosmology *cos, double omt, double m,
		       double shalf, double stwo, int flag)
{
  static Cosmology *c;
  static double omtilde,mfin;
  static double s2,sh;
  double mass;
  double s;

  if (flag == 1) {
    c = cos;
    omtilde = omt;
    mfin = m;
    sh = shalf;
    s2 = stwo;
    return;
  }

  s = stilde*(sh-s2)+s2;
  setInvertSigmaM(0.0,c,sqrt(s),1);
  mass = zriddrSimp(invertSigmaM,mfin/2.0,mfin,1.0e-6);
  result[1] = mfin/mass/sqrt(2.0*PI)/pow(stilde,1.5);
  result[1] *= exp(-omtilde*omtilde/2.0/stilde);
  result[1] *= -(1.0-omtilde*omtilde/stilde);
  //  cout << stilde << "\t" << result[1] << endl;
}

void formTimePDFInt(double s, double y[], double result[])
{
  setFormTimePDFInt(s,y,result,NULL,0.0,0.0,0.0,0.0,0);
}

/* Finds redshift at which a fraction FCRIT of the mass of a halo with
 * mass at zfin has assembled into pieces with a size a fraction
 * fInput of mass. */
double Cosmology::collapseRedshift(double mass, double zfin,
				   double fInput)
{
  const double FCRIT(0.5);
  double dcz0;
  double fMass;
  double sig0,sigf0;
  double temp;
  double zLess(0.0),zstep,zh,z;
  double dcz,frac,ratio,df;

  dcz0 = delCrit0(zfin)/growthFac(zfin);

  fMass = fInput*mass;
  sig0 = sigma0fM(mass,temp,0);
  sigf0 = sigma0fM(fMass,temp,0);

  zh = 200.0;
  z = zh;
  zstep = 10.0;

  while (1) {
    z = z - zstep;
    dcz = delCrit0(z)/growthFac(z);
    ratio = (dcz - dcz0)/(sqrt(2.0*(sigf0*sigf0 - sig0*sig0)));
    frac = erfc(ratio);
    df = fabs(frac - FCRIT)/FCRIT;
    if (df < 1.0e-5) {
      return z;
    }
    if (frac < FCRIT) {
      zLess = z;
    } else {
      zstep = zstep/2.0;
      z = zLess + zstep;
    }
  }
}

/************************************************************************
 ************************* Cooling Function *****************************
 ***********************************************************************/

// Calculates cooling rate for the halo, using data from coolingFile,
// to be specified in calling file.  Uses Sutherland and Dopita's table,
// with appropriate interpolation.  Using normalized cooling rate, not net.
// ORDER is default interpolation order, set in numrecipes.h.
/*
double calcCoolingRate(double temperature)
{
  char *coolingFile="/home/sfurlane/include/cooling/cool_met-nil.cie";
  double ltemp,trash;
  int lookup,index;
  double coolingRate,error;
  const int NPTS(91);    //  Make indices begin at 1 in array
  static double lt[NPTS+1],ne[NPTS+1],lambda[NPTS+1];
  static int iflag(0);

  if (iflag==0) {
    ifstream coolData(coolingFile);
    if (!coolData) {
      //      cerr << "coolData could not be opened\n";
      return 0.0;
    }
    for (int i=1;i<=NPTS;i++) {
      coolData >> lt[i] >> ne[i] >> trash >> trash >> trash >> lambda[i]
	       >> trash >> trash >> trash >> trash >> trash >> trash;
    }
    coolData.close();
    iflag = 1;
  }

  ltemp = log(temperature)/2.30259;
  locate(lt,NPTS,ltemp,lookup);
  index = fmin(fmax(lookup-(ORDER-1)/2,1),NPTS+1-ORDER);
  polint(&lt[index-1],&lambda[index-1],ORDER,ltemp,&coolingRate,&error);
  return pow(10.0,coolingRate);
}
*/

/************************************************************************
 ******************* Data Member Retrieval Functions ********************
 ***********************************************************************/

double Cosmology::getOmega0()
{
  return omega0;
}

double Cosmology::getOmegab()
{
  return omegab;
}

double Cosmology::getOmegam()
{
  return omegam;
}

double Cosmology::getOmbhh()
{
  return ombhh;
}

double Cosmology::getOm0hh()
{
  return om0hh;
}

double Cosmology::getLambda0()
{
  return lambda0;
}

double Cosmology::getw0()
{
  return w0;
}

double Cosmology::getH()
{
  return h;
}

double Cosmology::getnSpec()
{
  return nspec;
}

double Cosmology::getScale()
{
  return scale;
}

double Cosmology::getShape()
{
  return shape;
}

////////////////////////////////////////////////////////////////////
//  Misc distance related functions
////////////////////////////////////////////////////////////////////

// Conformal distance in Mpc, between two redshifts z and zp
// Units: r Mpc
double Cosmology::relativeConformalR(double z, double zp)
{
  return SPEEDOFLIGHT_CGS*(confTime(z)-confTime(zp))/H0_CMSMPC/getH();
}


// get the z corresponding to a conformal distance r from a point at redshift
// zbase.
double Cosmology::zOfR(double r, double zbase)
{
  double zmax;
  double tol(1.0e-5);
  double zofr;
  setZOfR(1.0,r,zbase,this,1);

  zmax=1.0001*zbase;
  while(relativeConformalR(zbase,zmax)<=r){
    zmax *=1.01;
    if(zmax>1000.0){
      cout <<"zmax too large in zofR"<<endl;
      return 1000.0;
    }
  }
  zofr= zriddrSimp(funcZOfR,zbase,zmax,tol);
  //  cout <<r <<"\t"<<zmax<<"\t"<<zofr<<endl;
  return zofr;
}

double setZOfR(double zp,double r1, double zbase, Cosmology *c1, int iflag)
{
  double temp;
  static Cosmology *c;
  static double z;
  static double r;
  if(iflag){
    c=c1;
    z=zbase;
    r=r1;
    return 0.0;
  }

  temp=r-c->relativeConformalR(z,zp);

  return temp;
}

//dummy function for rootsolving in zOfR
double funcZOfR(double z)
{
   return setZOfR(z,0.0,0.0,NULL,0);
}


/************************************************************************
 ************************ Halo Characteristics **************************
 ***********************************************************************/

/* Jeans mass; BL01 eq. 41 */
double jeansMass(Cosmology *c, double z)
{
  return (6.2*pow(c->getOm0hh(),-0.5)*pow(c->getOmbhh(),-0.6)*
	  pow(1.0+z,1.5));
}

/* Filtered mass, assuming a continuous Jeans mass and adiabatic
 * expansion past decoupling.  See Gnedin (2000) */
double filterMass(Cosmology *c, double z)
{
  double zt;
  double mass;

  zt = 137.0*pow(c->getOmbhh()/0.022,0.4) - 1.0;
  mass = 0.6*(1.0-2.0/3.0*sqrt((1.0+z)/(1.0+zt)));
  mass += log((1.0+zt)/(1.0+z))-2.0+2.0*sqrt((1.0+z)/(1.0+zt));
  mass = pow(3.0*mass,1.5);
  mass *= jeansMass(c,z);
  return mass;
}

/* Minimum mass for collapse assuming an ionized medium; condition is
 * Tvir > 2x10^5 K.  Uses BL01 relations. */
double minIonMass(Cosmology *c, double z)
{
  double mass;

  mass = 1.3e12/c->getH()/pow(1.0+z,1.5);
  mass *= sqrt(c->omegaZ(z)/c->getOmega0()/c->DeltaC(z));
  return mass;
}

/* Analytic solution for when a halo hits 10^4 K, given z.
 * Uses BL01 relations. */
double coolMass(Cosmology *c, double z)
{
  double m;

  m = 1.5e10/c->getH();
  m /= sqrt(c->getOmega0()/c->omegaZ(z)*c->DeltaC(z));
  m *= pow(1.0+z,-1.5);
  return m;
}

/* Analytic solution for when a halo hits 500 K, given z.
 * Uses BL01 relations. */
double coolMassH2(Cosmology *c, double z)
{
  double m;

  m = 1.69e8/c->getH();  //T=500K
  m /= sqrt(c->getOmega0()/c->omegaZ(z)*c->DeltaC(z));
  m *= pow(1.0+z,-1.5);
  return m;
}

/* Returns virial radius in kpc; mass is in solar masses.  BL01 eq. 24 */
double rvir(Cosmology *c, double mass, double z)
{
  double ans;

  ans = pow(mass*c->omegaZ(z)/c->getOmega0()/pow(c->getH(),2.0)/
	    c->DeltaC(z),1.0/3.0);
  ans *= 9.5e-2/(1.0+z);
  return ans;
}

/* Returns circular velocity in km/s; mass is in solar masses. BL01 eq. 25 */
double vcirc(Cosmology *c, double mass, double z)
{
  double ans;

  ans = pow(mass*mass/c->omegaZ(z)*c->getOmega0()*pow(c->getH(),2.0)*
	    c->DeltaC(z),1.0/6.0);
  ans *= 6.72e-3*sqrt(1.0+z);
  return ans;
}

/* Returns virial temperature in K; mass is in solar masses and mu
 * is mean molecular weight.  BL01 eq. 26 */
double tvir(Cosmology *c, double mass, double z, double mu)
{
  double ans;

  ans = pow(mass*mass/c->omegaZ(z)*c->getOmega0()*pow(c->getH(),2.0)*
	    c->DeltaC(z),1.0/3.0);
  ans *= 2.72e-3*mu*(1.0+z);
  return ans;
}

/* Convert a mass to comoving size */
double RComfromM(double m, Cosmology *c)
{
  double R;

  R = pow(3.0/4.0/PI*m/CRITDENMSOLMPC/c->getOm0hh(),1.0/3.0);
  return R;
}

/* Convert a comoving size to a mass */
double MfromRCom(double R, Cosmology *c)
{
  double m;

  m = 4.0*PI/3.0*R*R*R*CRITDENMSOLMPC*c->getOm0hh();
  return m;
}

/*************************************************************************
 *************************** Utility functions ***************************
 *************************************************************************/

/* The following are interpolated versions of the mass function etc, useful
 * for rapid evaluations. */

/* Halo mass function, in logarithmic units (per comoving Mpc^-3) */
double nm(double tM, double z, Cosmology *c)
{
  double dCritZ,sigM,dlsdlM,tdn,dsdM;

  dCritZ = c->delCrit0(z)/c->growthFac(z);
  sigM = sigm(tM,dsdM,c,1);
  dlsdlM = tM*dsdM/sigM;
  tdn = (sqrt(2.0/PI)*dCritZ*fabs(dlsdlM)*
	 exp(-dCritZ*dCritZ/(2.0*sigM*sigM))/(tM*sigM));
  tdn *= CRITDENMSOLMPC*c->getOm0hh();
  return tdn;
}

/* Calculates Sheth-Tormen (1999) mass function
 * in logarithmic units (per comoving Mpc^-3)
 *   tM = halo mass (Msun)
 */
double nmST(double tM, double z, Cosmology *c)
{
  const double a(0.707),A(0.322),p(0.3);
  double dCritZ,sigM,dlsdlM,tdn,dsdM,nuc;

  dCritZ = c->delCrit0(z)/c->growthFac(z);
  sigM = sigm(tM,dsdM,c,1);
  dlsdlM = tM*dsdM/sigM;
  nuc = dCritZ/sigM;
  tdn = A*sqrt(2.0*a/PI)*fabs(dlsdlM)/tM*nuc*(1.0+pow(nuc*nuc*a,-p));
  tdn *= exp(-a*nuc*nuc/2.0);
  tdn *= CRITDENMSOLMPC*c->getOm0hh();
  return tdn;
}

/* Conditional halo mass function, in logarithmic units
 * (per comoving Mpc^-3) */
double nmcond(double tM, double z, double mBubble, double deltaBubble,
	      Cosmology *c)
{
  double dCritZ,sigM,dlsdlM,tdn,dsdM;
  double sigB,sFac,dFac;

  dCritZ = c->delCrit0(z)/c->growthFac(z);
  dFac = dCritZ - deltaBubble;

  sigB = sigm(mBubble,dsdM,c,0);
  sigM = sigm(tM,dsdM,c,1);
  dlsdlM = tM*dsdM/sigM;
  sFac = sigM*sigM - sigB*sigB;

  tdn = sqrt(2.0/PI)*dFac*fabs(dlsdlM)*exp(-dFac*dFac/(2.0*sFac))/tM;
  tdn *= sigM*sigM/sFac/sqrt(sFac);
  tdn *= CRITDENMSOLMPC*c->getOm0hh();
  return tdn;
}

/* Halo bias (linear) */
double biasm(double mass, double z, Cosmology *c)
{
  double deltasc,deltasc0,sig;
  double bias;
  double dummy;

  deltasc = c->delCrit0(z)/c->growthFac(z);
  deltasc0 = c->delCrit0(0.0);
  sig = sigm(mass,dummy,c,0);
  bias = 1.0 + (deltasc*deltasc/sig/sig - 1.0)/deltasc0;

  return bias;
}

/* Halo bias (linear) */
double biasmST(double mass, double z, Cosmology *c)
{
  double q(0.75), p(0.3);
  double deltasc,deltasc0,sig;
  double bias;
  double dummy,nu;

  deltasc = c->delCrit0(z)/c->growthFac(z);
  deltasc0 = c->delCrit0(0.0);
  sig = sigm(mass,dummy,c,0);
  nu=deltasc*deltasc/sig/sig;
  bias = 1.0 + (q*nu - 1.0)/deltasc0;
  bias += 2.0*p/deltasc0/(1.0+pow(q*nu,p));

  return bias;
}

/* Interpolation of sigma(m).  Returns sigma.  Slightly nonstandard, if:
 *     flag=2, initializes table
 *     flag=3, does memory cleanup
 *     flag=0, just calculates and returns sigma
 *     flag=1, calculates and returns sigma, also dsdm=derivative
 * Parameters are
 *     m = mass (Msun)
 *     dsdm = variable in which to store derivative; ignored if flag
 *            does not equal 2
 * Uses spline interpolation from Numerical Recipes.
 */
double sigm(double m, double &dsdm, Cosmology *c, int flag)
{
  static double *lgmass,*sigm,*sigderiv;
  static double *sigsp,*sigderivsp;
  const double natural(1.0e30);
  int gridpts(300);
  static double mMin,mMax;
  double lgmMin,lgmMax;
  double deltaM;
  int i;
  double mass;
  double lgm;
  double sig,deriv;

  if (flag == 2) {  /* Initialize the spline fit */
    lgmass = dvector(1,gridpts);
    sigm = dvector(1,gridpts);
    sigsp = dvector(1,gridpts);
    sigderiv = dvector(1,gridpts);
    sigderivsp = dvector(1,gridpts);

    mMin = 1.0e3;
    mMax = 1.0e19;
    lgmMin = log(mMin);
    lgmMax = log(mMax);
    deltaM = (lgmMax-lgmMin)/(double)(gridpts);

    for (i=1;i<=gridpts;i++) {
      lgmass[i] = lgmMin + deltaM*(double)(i-1);
      mass = exp(lgmass[i]);
      sigm[i] = c->sigma0fM(mass,deriv,1);
      sigderiv[i] = deriv;
      //      cout << mass << "\t" << sigm[i] << "\t" << sigderiv[i] << endl;
    }

    spline(lgmass,sigm,gridpts,natural,natural,sigsp);
    spline(lgmass,sigderiv,gridpts,natural,natural,sigderivsp);

    return 0.0;
  }

  if (flag == 3) {
    free_dvector(lgmass,1,gridpts);
    free_dvector(sigm,1,gridpts);
    free_dvector(sigsp,1,gridpts);
    free_dvector(sigderiv,1,gridpts);
    free_dvector(sigderivsp,1,gridpts);
    return 0.0;
  }

  if (m > mMin && m < mMax) {
    lgm = log(m);
    splint(lgmass,sigm,sigsp,gridpts,lgm,&sig);
    if (flag == 1)
      splint(lgmass,sigderiv,sigderivsp,gridpts,lgm,&dsdm);
  } else {
    sig = c->sigma0fM(m,dsdm,flag);
  }
  return sig;
}

/**********************************************************************
 *************  Collapse rate *****************************************
 *********************************************************************/

/* Find mean rate at which matter is collapsing into bound objects */
double dfcdzPS(double z, Cosmology *c)
{
  double dfcdz;

  dfcdz = fabs((c->fCollPSExact(z+0.01,coolMass(c,z+0.01))-
                c->fCollPSExact(z-0.01,coolMass(c,z-0.01)))/0.02);
  return dfcdz;
}

/* Find mean rate at which matter is collapsing into bound objects */
double dfcdz(double z, Cosmology *c, int massfunc)
{
  double dfcdz;

  dfcdz = fabs((c->fColl(z+0.01,coolMass(c,z+0.01),massfunc)-
                c->fColl(z-0.01,coolMass(c,z-0.01),massfunc))/0.02);
  return dfcdz;
}

/* Find rate at which matter is collapsing in a region of size sig and
 * overdensity delta */
double dfcdz_region(double z, double sig, double delta, Cosmology *c)
{
  double dfcdz;
  double ml,mh,dcl,dch;
  double sigml,sigmh,dummy;
  double zl,zh;
  double fcl,fch;

  zl = z - 0.01;
  zh = z + 0.01;
  dcl = c->delCrit0(zl)/c->growthFac(zl);
  dch = c->delCrit0(zh)/c->growthFac(zh);

  ml = coolMass(c,zl);
  mh = coolMass(c,zh);
  sigml = sigm(ml,dummy,c,0);
  sigmh = sigm(mh,dummy,c,0);

  fcl = erfc((dcl-delta)/sqrt(2.0*(sigml*sigml-sig*sig)));
  fch = erfc((dch-delta)/sqrt(2.0*(sigmh*sigmh-sig*sig)));

  dfcdz = fabs((fcl-fch)/(zl-zh));
  return dfcdz;
}

// Redshift derivative of halo function
double nmDZ(double tM, double z, Cosmology *c)
{
  double nmdz;
  double dz(0.02);

  nmdz=fabs(nm(tM,z+dz/2.0,c)-nm(tM,z-dz/2.0,c))/dz;
  return nmdz;
}

/////////////////////////////////////////////////////////////////////
// Non-linear scale
/////////////////////////////////////////////////////////////////////
// defined by \sigma(R)=0.5  with R=PI/2/k
double Cosmology::nonlinearScale(double z)
{
	double tol(1.0e-4);
	double nonlin;
	double kmin(1.0e-3), kmax(1.0e5);
	setNonLinearScaleKernel(0.0,z,this,1);
	nonlin=zriddrSimp(getNonLinearScaleKernel,kmin,kmax,tol);

	return nonlin;
}

double setNonLinearScaleKernel(double k, double z1, Cosmology *c1, int iflag)
{
	static double z;
	static Cosmology *c;
	double sigma,dsdm,R,m;

	if(iflag==1){
		c=c1;
		z=z1;
		return 0.0;
	}
 	 R=PI/2.0/k;
	 m=MfromRCom(R,c);
	 sigma=sigm(m,dsdm,c,0)*c->growthFac(z);
	return sigma-0.5;
}

double getNonLinearScaleKernel(double k)
{
	return setNonLinearScaleKernel(k,0.0,NULL,0);
}

//JEong and Komatsu 2006 version
// defined by \Delta^2=0.4  with \Delta=P(k)*k^3/2 pi^2
double Cosmology::nonlinearScaleJK(double z)
{
	double tol(1.0e-4);
	double nonlin;
	double kmin(1.0e-3), kmax(1.0e5);

        return nonlinearScale(z);

      	setNonLinearScaleKernelJK(0.0,z,this,1);
      	nonlin=zriddrSimp(getNonLinearScaleKernelJK,kmin,kmax,tol);

	return nonlin;
}

double setNonLinearScaleKernelJK(double k, double z1, Cosmology *c1, int iflag)
{
	static double z;
	static Cosmology *c;
	double sigma,dsdm,R,m;
        double delta2;
        double nonlevel(0.4);
        double pk;

	if(iflag==1){
		c=c1;
		z=z1;
		return 0.0;
	}
        pk=c->powerSpectrum(k)*pow(c->growthFac(z),2.0);
        delta2=pk*k*k*k/2.0/PI/PI;
   	return delta2-nonlevel;
}

double getNonLinearScaleKernelJK(double k)
{
	return setNonLinearScaleKernelJK(k,0.0,NULL,0);
}
