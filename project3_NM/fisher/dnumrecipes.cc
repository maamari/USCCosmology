// dnumrecipes.cc 
// Contains utility, integration, interpolation, inversion, and ODE functions

#include <math.h>
#include <iostream>
#include <stdlib.h>
#include <string.h>
#include "dnumrecipes.h"
#include "dcomplex.h"
using namespace std;
/****************************************************************************
 ********** Functions for the Numerical Recipes Exception Handler ***********
 ***************************************************************************/

NumRecException::NumRecException(char *errorMessage)
{
  strcpy(message,errorMessage);
}

/****************************************************************************
 ******* Some utility functions used by all numerical recipes functions *****
 ***************************************************************************/

float *vector(long nl, long nh)
{
  float *v;

  v=(float *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(float)));
  return v-nl+NR_END;
}

void free_vector(float *v, long nl, long nh)
{
  free((FREE_ARG) (v+nl-NR_END));
}

float **matrix(long nrl, long nrh, long ncl, long nch)
{
  long i;
  long nrow(nrh-nrl+1);
  long ncol=(nch-ncl+1);
  float **m;

  /* allocate pointers to rows */
  m=(float **) malloc((size_t)((nrow+NR_END)*sizeof(float*)));
  if (!m)
    throw NumRecException("allocation failure 1 in matrix()");
  m += NR_END;
  m -= nrl;

  /* allocate rows and set pointers to them */
  m[nrl]=(float *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float)));
  if (!m[nrl])
    throw NumRecException("allocation failure 2 in matrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for (i=nrl+1;i<=nrh;i++) 
    m[i]=m[i-1]+ncol;

  /* return pointer to array of pointers to rows */
  return m;
}

void free_matrix(float **m, long nrl, long nrh, long ncl, long nch) 
{
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}

double *dvector(long nl, long nh)
{
  double *v;

  v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
  return v-nl+NR_END;
}

void free_dvector(double *v, long nl, long nh)
{
  free((FREE_ARG) (v+nl-NR_END));
}

double **dmatrix(long nrl, long nrh, long ncl, long nch)
{
  long i;
  long nrow(nrh-nrl+1);
  long ncol=(nch-ncl+1);
  double **m;

  /* allocate pointers to rows */
  m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
  if (!m)
    throw NumRecException("allocation failure 1 in dmatrix()");
  m += NR_END;
  m -= nrl;

  /* allocate rows and set pointers to them */
  m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
  if (!m[nrl])
    throw NumRecException("allocation failure 2 in dmatrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for (i=nrl+1;i<=nrh;i++) 
    m[i]=m[i-1]+ncol;

  /* return pointer to array of pointers to rows */
  return m;
}

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch) 
{
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}

//my modifications to steve's file
void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}

int *ivector(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
	int *v;

	v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
	if (!v) nrerror("allocation failure in ivector()");
	return v-nl+NR_END;
}

void free_ivector(int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}


/*********************************************************************
 ********************** Integration Routines *************************
 ********************************************************************/
// Integral from min to max using Romberg integration.  Converges provided
// f is continuous on the interval.  func is double precision.  
// accuracy is desired relative accuracy in the integral.
// Taken ultimately from CMBFAST (Seljak & Zaldarriaga 1996).
double qromb(double (*func)(double), double min, double max, double acc)
{
  double h,gmax,error,g0(0.0),g1,fourj;
  double g[MAXJ+1];
  int jmaxx,nint,i;
  
  h = 0.5*(max - min);
  gmax = h*((*func)(min) + (*func)(max));
  g[1] = gmax;
  nint = 1;
  error = 1.0e20;
  for (i=1;i <= JMAX;i++) {
    if (i > 5 && fabs(error) < acc)
      break;
    
    // Calculate next trapezoidal rule approx for integral
    g0 = 0.0;
    for (int k=1;k<=nint;k++) {
      g0 = g0 + (*func)(min+static_cast<double>(k+k-1)*h);
    }

    g0 = 0.5*g[1]+h*g0;
    h = 0.5*h;
    nint *= 2;
    jmaxx = fmin(i,MAXJ);
    fourj = 1.0;

    // Richardson extrapolation
    for (int j=1; j<=jmaxx; j++) {
      fourj *= 4.0;
      g1 = g0 + (g0-g[j])/(fourj-1.0);
      g[j]=g0;
      g0=g1;
    }
    
    if (fabs(g0) > acc) {
      error = 1.0 - gmax/g0;
    } else {
      error = gmax;
    }
    gmax = g0;
    g[jmaxx+1]=g0;
  }
  if (i > JMAX && fabs(error) > acc)
    throw NumRecException("rombint failed to converge");
  return g0;
}

// Numerical Recipes version of Romberg integration on a closed
// interval.  MAXJ is order of interpolation to be used.
double qrombnr(double (*func)(double), double a, double b, 
	       double accuracy)
{
  int j;
  double ss,dss,h[JMAX+2],s[JMAX+1];
  double sCurrent(0.0);

  h[1]=1.0;
  for (j=1;j<=JMAX;j++) {
    s[j] = trapzd(func,a,b,j,&sCurrent);
    if (j >= MAXJ) {
      polint(&h[j-MAXJ],&s[j-MAXJ],MAXJ,0.0,&ss,&dss);
      if (fabs(dss) <= accuracy*fabs(ss))
	return ss;
    }
    h[j+1] = 0.25*h[j];
  }
  throw NumRecException("Too many steps in routine qrombnr");
  return 0.0;
}

double qsimp(double (*func)(double), double a, double b,
	     double accuracy)
{
  int j;
  double s,st,ost,os;
  double sCurrent(0.0);

  ost = os = -1.0e30;
  for (j=1;j<=JMAX;j++) {
    st = trapzd(func,a,b,j,&sCurrent);
    s = (4.0*st - ost)/3.0;
    if (j > 5) {
      if (fabs(s-os) < accuracy*fabs(os) || (s==0.0 && os == 0.0))
	return s;
    }
    os = s;
    ost = st;
  }
  throw NumRecException("Too many steps in routine qsimp");
  return 0.0;
}

// Performs nth iteration of trapezoid rule; must be called in sequential 
// order of n's.  (*func) is the function to be integrated.
double trapzd(double (*func)(double), double lowerEndpoint, 
	     double upperEndpoint, int n, double *s)
{
  double x,tnm, sum, spacing;
  int numberOfEvals, j;

  if (n==1) {
    return (*s = 0.5*(upperEndpoint - lowerEndpoint)*
	    (((*func)(lowerEndpoint)) + ((*func)(upperEndpoint))));
  } else {
    for (numberOfEvals=1,j=1;j<(n-1);j++) 
      numberOfEvals <<= 1;
    tnm = numberOfEvals;
    spacing = (upperEndpoint - lowerEndpoint)/tnm;
    x = lowerEndpoint + 0.5*spacing;
    for (sum=0.0,j=1;j<=numberOfEvals;j++,x+=spacing)
      sum += (*func)(x);
    *s = 0.5*((*s) + (upperEndpoint - lowerEndpoint)*
			sum/tnm);
    return *s;
  }
}

// Numerical Recipes version of Romberg integration on an open interval
// (used for when one endpoint has a singularity).  MAXJ is order of 
// interpolation to be used; needs a midpoint rule with step-tripling
// passed as choose
double qromo(double (*func)(double), double a, double b, 
	     double (*choose)(double (*)(double), double, double, int, 
			      double (*)),
	     double accuracy)
{
  int j;
  double ss,dss,h[JMAXM+2],s[JMAXM+1];
  double sCurrent(0.0);

  h[1]=1.0;
  for (j=1;j<=JMAXM;j++) {
    s[j] = (*choose)(func,a,b,j,&sCurrent);
    if (j >= MAXJ) {
      polint(&h[j-MAXJ],&s[j-MAXJ],MAXJ,0.0,&ss,&dss);
      if (fabs(dss) <= accuracy*fabs(ss))
	return ss;
    }
    h[j+1] = h[j]/9.0;
  }
  throw NumRecException("Too many steps in routine qromo");
  return 0.0;
}

double qsimpo(double (*func)(double),double a, double b, 
	      double accuracy)
{
  int j;
  double s,st,ost,os;
  double sCurrent(0.0);

  ost = os = -1.0e30;
  for (j=1;j<=JMAXM;j++) {
    st = midpnt(func,a,b,j,&sCurrent);
    s = (9.0*st - ost)/8.0;
    if (j > 5) {
      if (fabs(s-os) < accuracy*fabs(os) || (s==0.0 && os == 0.0))
	return s;
    }
    os = s;
    ost = st;
  }
  throw NumRecException("Too many steps in routine qsimpo");
  return 0.0;
}

// Midpoint evaluation function for qromo, designed for an arbitrary 
// open integrand
double midpnt(double (*func)(double), double a, double b, int n,
	      double *s)
{
  double x,tnm,sum,del,ddel;
  int it,j;

  if (n==1) {
    x = 0.5*(a+b);
    *s = (b-a)*((*func)(x));
    return *s;
  } else {
    for (it=1,j=1;j<n-1;j++)
      it *= 3;
    tnm=it;
    del = (b-a)/(3.0*tnm);
    ddel = del+del;
    x = a + 0.5*del;
    sum = 0.0;
    for (j=1;j<=it;j++) {
      sum += (*func)(x);
      x += ddel;
      sum += (*func)(x);
      x += del;
    }
    *s = ((*s)+(b-a)*sum/tnm)/3.0;
    return *s;
  }
}

// Midpoint evaluation function for qromo, designed for use with an 
// integrand with a 1/sqrt(x-a) singularity in the lower limit
double midsql(double (*func)(double), double aa, double bb, int n,
	      double *s)
{
  double x,tnm,sum,del,ddel,a,b;
  int it,j;

  b = sqrt(bb-aa);
  a = 0.0;
  if (n==1) {
    x = 0.5*(a+b);
    *s = (b-a)*2.0*x*(*func)(aa+x*x);
    return *s;
  } else {
    for (it=1,j=1;j<n-1;j++)
      it *= 3;
    tnm = it;
    del = (b - a)/(3.0*tnm);
    ddel = del + del;
    x = a + 0.5*del;
    sum = 0.0;
    for (j=1;j<=it;j++) {
      sum += 2.0*x*(*func)(aa+x*x);
      x += ddel;
      sum += 2.0*x*(*func)(aa+x*x);
      x += del;      
    }
    *s = ((*s) + (b-a)*sum/tnm)/3.0;
    return *s;
  }
}

// Midpoint evaluation function for qromo, designed for use with an 
// integrand with a 1/sqrt(x-a) singularity in the lower limit
double midsqu(double (*func)(double), double aa, double bb, int n,
	      double *s)
{
  double x,tnm,sum,del,ddel,a,b;
  int it,j;

  b = sqrt(bb-aa);
  a = 0.0;
  if (n==1) {
    x = 0.5*(a+b);
    *s = (b-a)*2.0*x*(*func)(bb-x*x);
    return *s;
  } else {
    for (it=1,j=1;j<n-1;j++)
      it *= 3;
    tnm = it;
    del = (b - a)/(3.0*tnm);
    ddel = del + del;
    x = a + 0.5*del;
    sum = 0.0;
    for (j=1;j<=it;j++) {
      sum += 2.0*x*(*func)(bb-x*x);
      x += ddel;
      sum += 2.0*x*(*func)(bb-x*x);
      x += del;      
    }
    *s = ((*s) + (b-a)*sum/tnm)/3.0;
    return *s;
  }
}

/*********************************************************************
 ****************** Numerical Derivative Routines ********************
 ********************************************************************/

double dfridr(double (*func)(double), double x, double h, double *err)
{
  int i,j;
  double errt,fac,hh,**a,ans(-1.0e30);

  if (h == 0.0) 
    throw NumRecException("h must be nonzero in dfridr");
  a = dmatrix(1,NTAB,1,NTAB);
  hh = h;
  a[1][1] = ((*func)(x+hh) - (*func)(x-hh))/(2.0*hh);
  *err = BIG;
  for (i=2;i<=NTAB;i++) {
    hh /= CON;
    a[1][i] = ((*func)(x+hh)-(*func)(x-hh))/(2.0*hh);
    fac = CON*CON;
    for (j=2;j<=i;j++) {
      a[j][i] = (a[j-1][i]*fac-a[j-1][i-1])/(fac-1.0);
      fac = CON*CON*fac;
      errt=fmax(fabs(a[j][i]-a[j-1][i]),fabs(a[j][i]-a[j-1][i-1]));
      if (errt <= *err) {
	*err = errt;
	ans = a[j][i];
      }
    }
    if (fabs(a[i][i]-a[i-1][i-1]) >= 2.0*(*err))
      break;
  }
  free_dmatrix(a,1,NTAB,1,NTAB);
  return ans;
}

/*********************************************************************
 ********************** Interpolation Routines ***********************
 ********************************************************************/

// Performs polynomial extrapolation using Neville's algorithm.  xa and ya
// are the dependent/independent matrices off of which we extrapolate at
// point x; result is estimated extrapolation and errorEstimate is estimated
// error.  From Numerical Recipes section 3.1.  n is number of points in
// xa and ya (using unit offset - so if you input part of an array, send 
// it in with one less than the real beginning)
void polint(double xa[], double ya[], int n, double x, double *y,
	    double *dy)
{
  int i,m,ns=1;
  double den,dif,dift,ho,hp,w;
  double *c,*d;                   // corrections to apply at each order

  dif = fabs(x-xa[1]);
  c = dvector(1,n);
  d = dvector(1,n);
  for (i=1;i<=n;i++) {           // find ns to x
    if ( (dift = fabs(x-xa[i])) < dif) {
      ns = i;
      dif = dift;
    }
    c[i] = ya[i];
    d[i] = ya[i];
  }
  *y = ya[ns--];       // initial guess
  for (m=1;m<n;m++) {            // loop over columns of tableau
    for (i=1;i<=n-m;i++) {         // update c,d values
      ho = xa[i] - x;
      hp = xa[i+m] - x;
      w = c[i+1] - d[i];
      if ( (den = ho - hp) == 0.0) {
	throw NumRecException("Error in routine polint");
      }
      den = w/den;
      d[i] = hp*den;
      c[i] = ho*den;
    }
    // Next line decides straightest route to take through tableau
    *y += (*dy = (2*ns < (n-m) ? c[ns+1] : d[ns--]));
  }
  free_dvector(d,1,n);
  free_dvector(c,1,n);
}

// Given an array, sorted in either increasing or decreasing order, this 
// finds the lower of the two indices bracketing the value x as j.  Expects 
// array to be in form xx[1...n].
void locate(double xx[], int n, double x, int &j)
{
  int ju,jm,jl;
  int ascnd;

  jl = 0;
  ju = n+1;
  ascnd = (xx[n] >= xx[1]);
  while ((ju-jl) > 1) {
    jm = (ju+jl) >> 1;
    if (x >= xx[jm] == ascnd) {
      jl = jm;
    } else {
      ju = jm;
    }
    //    cout << "jl = " << jl << " ju = " << ju << endl;
  }
  if (x == xx[1]) {
    j = 1;
  }  else {
    if (x == xx[n]) {
      j = n-1;
    } else {
      j = jl;
    }
  }
}

/*********************************************************************
 ********************** Root-finding Routines ************************
 ********************************************************************/

// Given a function func and a guessed interval in which the function 
// changes sign, zbrac expands the interval until it really does change
// sign.  func should have the form
// double func(const double x-value, double *y-value, float *deriv, int flag)
// If flag=0, deriv is not calculated; if flag=1, it is calculated.
// Input x1 and x2 are changed in the function to the bracketing values.
int zbrac(void (*func)(const double, double *, double *, int), double &x1, 
	  double &x2)
{
  int j;
  double f1,f2,dumb;

  if (x1==x2) 
    throw NumRecException("Bad initial range in zbrac");
  (*func)(x1,&f1,&dumb,0);
  (*func)(x2,&f2,&dumb,0);
  for (j=1;j<=MAXIT;j++) {
    if (f1*f2 < 0.0) 
      return(1);
    if (fabs(f1) < fabs(f2)) {
      x1 += FACTOR*(x1-x2);
      (*func)(x1,&f1,&dumb,0);
    } else {
      x2 += FACTOR*(x2-x1);
      (*func)(x2,&f2,&dumb,0);
    }
  }
  return(0);
}

// Same as zbrac, but assumes a simpler input function
int zbracSimp(double (*func)(const double), double &x1, double &x2)
{
  int j;
  double f1,f2;

  if (x1==x2) 
    throw NumRecException("Bad initial range in zbrac");
  f1 = (*func)(x1);
  f2 = (*func)(x2);
  for (j=1;j<=MAXIT;j++) {
    if (f1*f2 < 0.0) 
      return(1);
    if (fabs(f1) < fabs(f2)) {
      x1 += FACTOR*(x1-x2);
      f1 = (*func)(x1);
    } else {
      x2 += FACTOR*(x2-x1);
      f2 = (*func)(x2);
    }
  }
  return(0);
}

double zriddr(void (*func)(const double, double *, double *, int),
	      double x1, double x2, double xacc)
{ 
  int j;
  double ans,fh,fl,fm,fnew,s,xh,xl,xm,xnew,dumb;

  (*func)(x1,&fl,&dumb,0);
  (*func)(x2,&fh,&dumb,0);
  if ((fl > 0.0 && fh < 0.0) || (fl < 0.0 && fh > 0.0)) {
    xl = x1;
    xh = x2;
    ans = -1.11e30;
    for (j=1;j<=MAXIT;j++) {
      xm = 0.5*(xl + xh);
      (*func)(xm,&fm,&dumb,0);
      s = sqrt(fm*fm - fl*fh);
      if (s == 0.0)
	return ans;
      xnew = xm + (xm-xl)*((fl >= fh ? 1.0 : -1.0)*fm/s);
      if (fabs(xnew-ans) <= xacc) 
	return ans;
      ans = xnew;
      (*func)(ans,&fnew,&dumb,0);
      if (fnew == 0.0) 
	return ans;
      if (sign(fm,fnew) != fm) {
	xl = xm;
	fl = fm;
	xh = ans;
	fh = fnew;
      } else {
	if (sign(fl,fnew) != fl) {
	  xh = ans;
	  fh = fnew;
	} else {
	  if (sign(fh,fnew) != fh) {
	    xl = ans;
	    fl = fnew;
	  } else {
	    throw NumRecException("never get here");
	  }
	}
      }
      if (fabs(xh-xl) <= xacc)
	return ans;
    }
    throw NumRecException("Maximum number of iterations exceeded in zriddr");
  } else {
    if (fl == 0.0) 
      return x1;
    if (fh == 0.0)
      return x2;
    throw NumRecException("root must be bracketed in zriddr");
  }
  return 0.0;
}

// Same as zriddr, but function isn't expected to calculate its derivative
double zriddrSimp(double (*func)(const double), 
		  double x1, double x2, double xacc)
{ 
  int j;
  double ans,fh,fl,fm,fnew,s,xh,xl,xm,xnew;

  fl = (*func)(x1);
  fh = (*func)(x2);
  if ((fl > 0.0 && fh < 0.0) || (fl < 0.0 && fh > 0.0)) {
    xl = x1;
    xh = x2;
    ans = -1.11e30;
    for (j=1;j<=MAXIT;j++) {
      xm = 0.5*(xl + xh);
      fm = (*func)(xm);
      s = sqrt(fm*fm - fl*fh);
      if (s == 0.0)
	return ans;
      xnew = xm + (xm-xl)*((fl >= fh ? 1.0 : -1.0)*fm/s);
      if (fabs(xnew-ans) <= xacc) 
	return ans;
      ans = xnew;
      fnew = (*func)(ans);
      if (fnew == 0.0) 
	return ans;
      if (sign(fm,fnew) != fm) {
	xl = xm;
	fl = fm;
	xh = ans;
	fh = fnew;
      } else {
	if (sign(fl,fnew) != fl) {
	  xh = ans;
	  fh = fnew;
	} else {
	  if (sign(fh,fnew) != fh) {
	    xl = ans;
	    fl = fnew;
	  } else {
	    throw NumRecException("never get here");
	  }
	}
      }
      if (fabs(xh-xl) <= xacc)
	return ans;
    }
    throw NumRecException("Maximum number of iterations exceeded in zriddr");
  } else {
    if (fl == 0.0) 
      return x1;
    if (fh == 0.0)
      return x2;
    throw NumRecException("root must be bracketed in zriddr");
  }
  return 0.0;
}

// Same as zriddr, but function isn't expected to calculate its derivative
double zriddrConst(double (*func)(const double), double yval, 
		   double x1, double x2, double xacc)
{ 
  int j;
  double ans,fh,fl,fm,fnew,s,xh,xl,xm,xnew;

  fl = (*func)(x1) - yval;
  fh = (*func)(x2) - yval;
  if ((fl > 0.0 && fh < 0.0) || (fl < 0.0 && fh > 0.0)) {
    xl = x1;
    xh = x2;
    ans = -1.11e30;
    for (j=1;j<=MAXIT;j++) {
      xm = 0.5*(xl + xh);
      fm = (*func)(xm) - yval;
      s = sqrt(fm*fm - fl*fh);
      if (s == 0.0)
	return ans;
      xnew = xm + (xm-xl)*((fl >= fh ? 1.0 : -1.0)*fm/s);
      if (fabs(xnew-ans) <= xacc) 
	return ans;
      ans = xnew;
      fnew = (*func)(ans) - yval;
      if (fnew == 0.0) 
	return ans;
      if (sign(fm,fnew) != fm) {
	xl = xm;
	fl = fm;
	xh = ans;
	fh = fnew;
      } else {
	if (sign(fl,fnew) != fl) {
	  xh = ans;
	  fh = fnew;
	} else {
	  if (sign(fh,fnew) != fh) {
	    xl = ans;
	    fl = fnew;
	  } else {
	    throw NumRecException("never get here");
	  }
	}
      }
      if (fabs(xh-xl) <= xacc)
	return ans;
    }
    throw NumRecException("Maximum number of iterations exceeded in zriddr");
  } else {
    if (fl == 0.0) 
      return x1;
    if (fh == 0.0)
      return x2;
    throw NumRecException("root must be bracketed in zriddr");
  }
  return 0.0;
}

// rtsafe finds a root in the interval x1 -> x2 for funcd, which should have
// same form as func in zbrac, to accuracy xacc.
double rtsafe(void (*funcd)(const double, double *, double *, int), double x1, 
	     double x2, double xacc)
{
  int j;
  double df,dx,dxold,f,fh,fl;
  double temp,xh,xl,rts;

  (*funcd)(x1,&fl,&df,0);
  (*funcd)(x2,&fh,&df,0);
  if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0))
    throw NumRecException("Root must be bracketed in rtsafe");
  if (fl == 0.0) return(x1);
  if (fh == 0.0) return(x2);
  if (fl < 0.0) {
    xl=x1;
    xh=x2;
  } else {
    xh=x1;
    xl=x2;
  }
  rts=0.5*(x1+x2);
  dxold=fabs(x2-x1);
  dx=dxold;
  (*funcd)(rts,&f,&df,1);
  for (j=1;j<=MAXIT;j++) {
    if ((((rts-xh)*df-f)*((rts-xl)*df-f) > 0.0)
	|| (fabs(2.0*f) > fabs(dxold*df))) {
      dxold=dx;
      dx=0.5*(xh-xl);
      rts=xl+dx;
      if (xl == rts) 
	return(rts);
    } else {
      dxold=dx;
      dx=f/df;
      temp=rts;
      rts -= dx;
      if (temp == rts) 
	return(rts);
    }
    if (fabs(dx) < xacc) 
      return(rts);
    (*funcd)(rts,&f,&df,1);
    if (f < 0.0)
      xl=rts;
    else 
      xh=rts;
  }
  throw NumRecException("Max number of iterations exceeded in rtsafe");
  return(0.0);
}

/*********************************************************************
 ********************** Minimization Routines ************************
 ********************************************************************/

double brent(double ax, double bx, double cx, double (*f)(double),
	     double tol, double *xmin)
{
  const int ITMAX(100);
  const double CGOLD(0.3819660),ZEPS(1.0e-10);

  int iter;
  double a,b,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
  double d(1.0),e(0.0);

  a = (ax < cx ? ax : cx);
  b = (ax > cx ? ax : cx);
  x = w = v = bx;
  fw = fv = fx = (*f)(x);
  for (iter=1;iter<=ITMAX;iter++) {
    xm = 0.5*(a+b);
    tol2 = 2.0*(tol1=tol*fabs(x)+ZEPS);
    if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
      *xmin = x;
      return fx;
    }
    if (fabs(e) > tol1) {
      r = (x-w)*(fx-fv);
      q = (x-v)*(fx-fw);
      p = (x-v)*q-(x-w)*r;
      q = 2.0*(q-r);
      if (q > 0.0)
	p = -p;
      q = fabs(q);
      etemp = e;
      e = d;
      if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) ||
	  p >= q*(b-x)) {
	d = CGOLD*(e = (x >= xm ? a-x : b-x));
      } else {
	d = p/q;
	u = x+d;
	if (u-a < tol2 || b-u < tol2)
	  d = sign(tol1,xm-x);
      }
    } else {
      d = CGOLD*(e = (x >= xm ? a-x : b-x));
    }
    u = (fabs(d) >= tol1 ? x+d : x+sign(tol1,d));
    fu = (*f)(u);
    if (fu <= fx) {
      if (u >= x)
	a = x;
      else
	b = x;
      v = w;
      w = x;
      x = u;
      fv = fw;
      fw = fx;
      fx = fu;
    } else {
      if (u < x)
	a = u;
      else
	b = u;
      if (fu <= fw || w == x) {
	v = w;
	w = u;
	fv = fw;
	fw = fu;
      } else if (fu <= fv || v == x || v == w) {
	v = u;
	fv = fu;
      }
    }
  }
  throw NumRecException("Max number of iterations exceeded in brent");
  *xmin = x;
  return fx;
}

/*********************************************************************
 ******************* Elliptic Integral Routines **********************
 ********************************************************************/

// Evaluates Legendre elliptic integral of the first kind by calling rf.
// Argument ranges are 0 <= phi <= pi/2, 0 <= k*sin(phi) <= 1.  Note this
// definition of the elliptical function differs from the 
// Mathematica version!  In the integrand, it is k^2, not k.
double ellf(double phi, double ak)
{
  double s;
  s = sin(phi);
  return (s*rf(cos(phi)*cos(phi),(1.0-s*ak)*(1.0+s*ak),1.0));
}

double ellfsq(double phi, double ak2)
{
  double s;

  s = sin(phi);
  return (s*rf(cos(phi)*cos(phi),(1.0-s*s*ak2),1.0));
}

// Evaluates Legendre elliptic integral of the second kind by calling
// rf and rd.  Argument ranges are 0 <= phi <= pi/2, 0 <= k*sin(phi) <= 1.
double elle(double phi, double ak)
{
  double cc,q,s;

  s = sin(phi);
  cc = pow(cos(phi),2.0);
  q = (1.0 - s*ak)*(1.0 + s*ak);
  return s*(rf(cc,q,1.0)-pow(s*ak,2.0)*rd(cc,q,1.0)/3.0);
}

// Evaluates Carlson's elliptic integral of the first kind.  x,y,z must
// be non-negative, and at most one can be zero.
double rf(const double x, const double y, const double z)
{
  const double c1(1.0/24.0),c2(0.1),c3(3.0/44.0),c4(1.0/14.0);
  double alamb,ave,delx,dely,delz,e2,e3,sqrtx,sqrty,sqrtz;
  double xt,yt,zt;

  if (fmin(fmin(x,y),z) < 0.0 || fmin(fmin(x+y,x+z),y+z) < TINY ||
      fmax(fmax(x,y),z) > BIG)
    throw NumRecException("Invalid Arguments in rf");
  xt = x;
  yt = y;
  zt = z;
  do {
    sqrtx = sqrt(xt);
    sqrty = sqrt(yt);
    sqrtz = sqrt(zt);
    alamb = sqrtx*(sqrty+sqrtz)+sqrty*sqrtz;
    xt = 0.25*(xt+alamb);
    yt = 0.25*(yt+alamb);
    zt = 0.25*(zt+alamb);
    ave = 1.0/3.0*(xt+yt+zt);
    delx = (ave-xt)/ave;
    dely = (ave-yt)/ave;
    delz = (ave-zt)/ave;
  } while (fmax(fmax(fabs(delx),fabs(dely)),fabs(delz)) > ERRTOL);
  e2 = delx*dely - delz*delz;
  e3 = delx*dely*delz;
  return ((1.0+(c1*e2-c2-c3*e3)*e2+c4*e3)/sqrt(ave));
}

double rd(const double x, const double y, const double z)
{
  double alamb,ave,delx,dely,delz,ea,eb,ec,ed,ee,fac;
  double sqrtx,sqrty,sqrtz,sum,xt,yt,zt,ans;
  const double c1(3.0/14.0),c2(1.0/6.0),c3(9.0/22.0);
  const double c4(3.0/26.0),c5(9.0/88.0),c6(4.5/26.0);

  if (fmin(x,y) < 0.0 || fmin(x+y,z) < TINY2 || 
      fmax(fmax(x,y),z) > BIG) 
    throw NumRecException("invalid arguments in rd");
  xt = x;
  yt = y;
  zt = z;
  sum = 0.0;
  fac = 1.0;
  do {
    sqrtx = sqrt(xt);
    sqrty = sqrt(yt);
    sqrtz = sqrt(zt);
    alamb = sqrtx*(sqrty+sqrtz)+sqrty*sqrtz;
    sum += fac/(sqrtz*(zt+alamb));
    fac = 0.25*fac;
    xt = 0.25*(xt+alamb);
    yt = 0.25*(yt+alamb);
    zt = 0.25*(zt+alamb);
    ave = 0.2*(xt+yt+3.0*zt);
    delx = (ave-xt)/ave;
    dely = (ave-yt)/ave;
    delz = (ave-zt)/ave;
  } while (fmax(fmax(fabs(delx),fabs(dely)),fabs(delz)) > ERRTOL);
  ea = delx*dely;
  eb = delz*delz;
  ec = ea - eb;
  ed = ea - 6.0*eb;
  ee = ed + ec + ec;
  ans = 1.0+ed*(-c1+c5*ed-c6*delz*ee);
  ans += delz*(c2*ee+delz*(-c3*ec+delz*c4*ea));
  ans *= fac/(ave*sqrt(ave));
  ans += 3.0*sum;
  return ans;
}

/*********************************************************************
 ********************* Other Special Functions ***********************
 ********************************************************************/

double bessk(int n, double x)
{
  int j;
  double bk,bkm,bkp,tox;

  if (n == 0) {
    return bessk0(x);
  } else if (n == 1) {
    return bessk1(x);
  }
  tox = 2.0/x;
  bkm = bessk0(x);
  bk = bessk1(x);
  for (j=1;j<n;j++) {
    bkp = bkm + j*tox*bk;
    bkm = bk;
    bk = bkp;
  }
  return bk;
}

double bessk0(double x)
{
  double y,ans;

  if (x <= 2.0) {
    y = x*x/4.0;
    ans = ((-log(x/2.0)*bessi0(x)) +
	   -0.57721566 + y*(0.42278420 + 
			    y*(0.23069756 +
			       y*(0.3488590e-1 + 
				  y*(0.262698e-2 + 
				     y*(0.10750e-3 + 
					y*0.74e-5))))));
  } else {
    y = 2.0/x;
    ans = ((exp(-x)/sqrt(x))*
	    (1.25331414 + 
	     y*(-0.7832358e-1 + 
		y*(0.2189568e-1 + 
		   y*(-0.1062446e-1 +
		      y*(0.587872e-2 + 
			 y*(-0.251540e-2 + 
			    y*0.53208e-3)))))));
  }
  return ans;
}

double bessk1(double x)
{
  double y,ans;

  if (x <= 2.0) {
    y = x*x/4.0;
    ans = ((log(x/2.0)*bessi1(x)) + 
	   (1.0/x)*(1.0 + 
		    y*(0.15443144 + 
		       y*(-0.67278579 + 
			  y*(-0.18156897 + 
			     y*(-0.1919402e-1 + 
				y*(-0.110404e-2 + 
				   y*(-0.4686e-4))))))));
  } else {
    y = 2.0/x;
    ans = ((exp(-x)/sqrt(x))*
	   (1.25331414 + 
	    y*(0.23498619 + 
	       y*(-0.3655620e-1 + 
		  y*(0.1504268e-1 + 
		     y*(-0.780353e-2 + 
			y*(0.325614e-2 + 
			   y*(-0.68245e-3))))))));
  } 
  return ans;
}

double bessi1(double x)
{
  double ax,ans;
  double y;

  if ((ax = fabs(x)) < 3.75) {
    y = x/3.75;
    y *= y;
    ans = ax*(0.5 + 
	      y*(0.87890594 + 
		 y*(0.51498869 + 
		    y*(0.15084934 + 
		       y*(0.2658733e-1 + 
			  y*(0.301532e-2 + 
			     y*0.32411e-3))))));
  } else {
    y = 3.75/ax;
    ans = (0.2282967e-1 + y*(-0.2895312e-1 +  
			     y*(0.1787654e-1 - 
				y*0.420059e-2)));
    ans = (0.39894228 + 
	   y*(-0.3988024e-1 + 
	      y*(-0.362018e-2 + 
		 y*(0.163801e-2 + 
		    y*(-0.1031555e-1 + y*ans)))));
    ans *= (exp(ax)/sqrt(ax));
  }
  return (x < 0.0 ? -ans : ans);
}

double bessj0(double x)
{
  double ax,z;
  double xx,y,ans,ans1,ans2;

  if ((ax=fabs(x)) < 8.0) {
    y=x*x;
    ans1=(57568490574.0
	  +y*(-13362590354.0
	      +y*(651619640.7
		  +y*(-11214424.18+y*(77392.33017+y*(-184.9052456))))));
    ans2=(57568490411.0
	  +y*(1029532985.0
	      +y*(9494680.718
		  +y*(59272.64853+y*(267.8532712+y*1.0)))));
    ans=ans1/ans2;
  } else {
    z=8.0/ax;
    y=z*z;
    xx=ax-0.785398164;
    ans1=1.0+y*(-0.1098628627e-2
		+y*(0.2734510407e-4
		    +y*(-0.2073370639e-5+y*0.2093887211e-6)));
    ans2 = (-0.1562499995e-1
	    +y*(0.1430488765e-3
		+y*(-0.6911147651e-5
		    +y*(0.7621095161e-6-y*0.934935152e-7))));
    ans=sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
  }
  return ans;
}

double bessi0(double x)
{
  double ax,ans;
  double y;

  if ((ax=fabs(x)) < 3.75) {
    y = x/3.75;
    y *= y;
    ans=1.0+y*(3.5156229+
	       y*(3.0899424+
		  y*(1.2067492+
		     y*(0.2659732+
			y*(0.360768e-1+y*0.45813e-2)))));
  } else {
    y = 3.75/ax;
    ans=((exp(ax)/sqrt(ax))*
	 (0.39894228+
	  y*(0.1328592e-1+
	     y*(0.225319e-2+
		y*(-0.157565e-2+
		   y*(0.916281e-2
		      +y*(-0.2057706e-1+
			  y*(0.2635537e-1+
			     y*(-0.1647633e-1+y*0.392377e-2)))))))));
  }
  return ans;
}


// Returns complementary error function, 1-erf, using a Chebyshev fit
// from Numerical Recipes
double erfcc(double x)
{
  double t,z,ans;

  z = fabs(x);
  t = 1.0/(1.0+0.5*z);
  ans = t*exp(-z*z-1.26551223+
	      t*(1.00002368+
		 t*(0.37409196+
		    t*(0.09678418+
		       t*(-0.18628806+
			  t*(0.27886807+
			     t*(-1.13520398+
				t*(1.48851587+
				   t*(-0.82215223+
				      t*0.17087277)))))))));
  return (x >= 0.0 ? ans : 2.0-ans);		     
}

// Returns error function, using erfcc above
double erff(double x)
{
  return (1.0-erfcc(x));
}

/* Computes Si(x) and Ci(x); note if x<0 Ci(x) does NOT include the 
 * -i*PI factor. From Numerical Recipes. */
void cisi(double x, double *ci, double *si)
{
  int i,k,odd;
  double a,err,fact(1.0),sign(1.0),sum,sumc,sums(0.0),t,term;
  complex h,b,c,d,del;

  t = fabs(x);
  if (t == 0.0) {
    *si = 0.0;
    *ci = -1.0/TINY;
  }
  if (t > TMIN) {
    b = complex(1.0,t);
    c = complex(1.0/TINY,0.0);
    d = h = 1.0/b;
    for (i=2;i<=MAXIT;i++) {
      a = -(i-1)*(i-1);
      b += 2.0;
      d = 1.0/(a*d+b);
      c = b + a/c;
      del = c*d;
      h *= del;
      if ((fabs(del.real() - 1.0) + fabs(del.imaginary())) < EPS)
	break;
    }
    if (i > MAXIT)
      throw NumRecException("cf failed in cisi");
    h *= complex(cos(t),-sin(t));
    *ci = -h.real();
    *si = PI/2.0 + h.imaginary();
  } else {
    if (t < sqrt(TINY)) {
      sumc = 0.0;
      sums = t;
    } else {
      sum = sums = sumc = 0.0;
      sign = fact = 1.0;
      odd = TRU;
      for (k=1;k<=MAXIT;k++) {
	fact *= t/k;
	term = fact/k;
	sum += sign*term;
	err = term/fabs(sum);
	if (odd) {
	  sign = -sign;
	  sums = sum;
	  sum = sumc;
	} else {
	  sumc = sum;
	  sum = sums;
	} 
	if (err < EPS) 
	  break;
	odd = !odd;
      }
      if (k > MAXIT)
	throw NumRecException("maxits exceeded in cisi");
    }
    *si = sums;
    *ci = sumc + log(t) + EULER;
  }
  if (x < 0.0)
    *si = -(*si);
}

/*********************************************************************
 ******************* Random Number Generators ************************
 ********************************************************************/

float ran1(long *idum)
{
  const int IA(16807),IM(2147483647),IQ(127773),IR(2836);
  const float AM(1.0/IM),REPS(1.2e-7),RNMX(1.0-REPS);
  const int NTAB(32),NDIV(1+(IM-1)/NTAB);

  int j;
  long k;
  static long iy(0);
  static long iv[NTAB];
  float temp;

  if (*idum <= 0 || !iy) {
    if (-(*idum) < 1) 
      *idum=1;
    else 
      *idum = -(*idum);
    for (j=NTAB+7;j>=0;j--) {
      k = (*idum)/IQ;
      *idum = IA*(*idum-k*IQ)-IR*k;
      if (*idum < 0)
	*idum += IM;
      if (j < NTAB)
	iv[j] = *idum;
    }
    iy = iv[0];
  }
  k = (*idum)/IQ;
  *idum = IA*(*idum-k*IQ)-IR*k;
  if (*idum < 0)
    *idum += IM;
  j = iy/NDIV;
  iy = iv[j];
  iv[j] = *idum;
  if ((temp=AM*iy) > RNMX)
    return RNMX;
  else
    return temp;
}

double gasdev(long *idum)
{
  static int iset(0);
  static double gset;
  double fac,rsq,v1,v2;

  if (*idum < 0) iset=0;
  if (iset == 0) {
    do {
      v1 = 2.0*ran1(idum)-1.0;
      v2 = 2.0*ran1(idum)-1.0;
      rsq = v1*v1 + v2*v2;
    } while (rsq >= 1.0 || rsq == 0.0);
    fac = sqrt(-2.0*log(rsq)/rsq);
    gset = v1*fac;
    iset = 1;
    return v2*fac;
  } else {
    iset = 0;
    return gset;
  }
}

/*********************************************************************
 ******************* ODE Integration Routines ************************
 ********************************************************************/

/* Main driver; modified from Numerical Recipes form in order to allow
 * cascading (replacing static variables in bsstep with a private
 * data structure). */
void odeint(double ystart[], int nvar, double x1, double x2, double eps, 
	    double h1, double hmin, int *nok, int *nbad,
	    void  (*derivs)(double, double [], double []),
	    void (*rkqs)(driverData, double [], double [], int, double *, 
			 double, double, double [], double *, double *, 
			 void (*)(double, double [], double [])))
{
  int nstp,i;
  double x,hnext,hdid,h;
  double *yscal,*y,*dydx;

  if (x1 == x2) {
    return;
  }
  yscal=dvector(1,nvar);
  y=dvector(1,nvar);
  dydx=dvector(1,nvar);
  x=x1;
  h=sign(h1,x2-x1);
  *nok = 0;
  *nbad = 0;
  struct driverData bsData;
  bsData.first = 1;
  bsData.epsold = -1.0;
  for (i=1;i<=nvar;i++) {
    y[i]=ystart[i];
  }
  for (nstp=1;nstp<=MAXSTP;nstp++) {
    (*derivs)(x,y,dydx);
    for (i=1;i<=nvar;i++) {
      yscal[i]=fabs(y[i])+fabs(dydx[i]*h)+TINY;
      //      cout << "yscal[i] = " << yscal[i] << endl;
    }
    if ((x+h-x2)*(x+h-x1) > 0.0) 
      h=x2-x;
    try {
      (*rkqs)(bsData,y,dydx,nvar,&x,h,eps,yscal,&hdid,&hnext,derivs);
    }
    catch( NumRecException ex ) {
      for (i=1;i<=nvar;i++)
	ystart[i]=y[i];
      free_dvector(dydx,1,nvar);
      free_dvector(y,1,nvar);
      free_dvector(yscal,1,nvar);
      cout<<"Exception processed in odeint"<<endl;
      cout << "Exception processed in odeint: " << ex.what() << endl;
      throw;
    }
    if (hdid==h) { 
      ++(*nok);
    } else {
      ++(*nbad);
    }
    if ((x-x2)*(x2-x1) >= 0.0) {
      for (i=1;i<=nvar;i++)
	ystart[i]=y[i];
      free_dvector(dydx,1,nvar);
      free_dvector(y,1,nvar);
      free_dvector(yscal,1,nvar);
      return;
    }
    if (fabs(hnext) <= hmin) {
      free_dvector(dydx,1,nvar);
      free_dvector(y,1,nvar);
      free_dvector(yscal,1,nvar);
      cout<<"Error in odeint"<<endl;
      throw NumRecException("Error in odeint");
    }
    h=hnext;
  }
  free_dvector(dydx,1,nvar);
  free_dvector(y,1,nvar);
  free_dvector(yscal,1,nvar);
  cout<<"Too many steps in odeint"<<endl;
  throw NumRecException("Too many steps in odeint");
}  

void rkqs(driverData self, double y[], double dydx[], int n, double *x, 
	  double htry, double eps, double yscal[], double *hdid, 
	  double *hnext, void (*derivs)(double, double [], double []))
{
  void rkck(double y[], double dydx[], int n, double x, double h, 
	    double yout[],double yerr[], 
	    void (*derivs)(double, double [], double []));

  int i;
  double errmax;
  double h, htemp;
  double xnew;
  double *yerr;
  double *ytemp;

  yerr=dvector(1,n);
  ytemp=dvector(1,n);
  h=htry;
  for (;;) {
    rkck(y,dydx,n,*x,h,ytemp,yerr,derivs);
    errmax=0.0;
    for (i=1;i<=n;i++)
      errmax=fmax(errmax,fabs(yerr[i]/yscal[i]));
    errmax /= eps;
    if (errmax <= 1.0) 
      break;
    htemp=SAFETY*h*pow(errmax,PSHRNK);
    if (h >=0.0) {
      h=fmax(htemp,0.1*h);
    } else {
      h=fmin(htemp,0.1*h);
    }
    xnew=(*x)+h;
    if (xnew == *x) {
      for (i=1;i<=n;i++)
	y[i]=ytemp[i];
      *x += (*hdid=h);
      free_dvector(ytemp,1,n);
      free_dvector(yerr,1,n);
      cout<<"Underflow occurred in rkqs"<<endl;
      throw NumRecException("Underflow occurred in rkqs");
    }
  }
  if (errmax > ERRCON) 
    *hnext=SAFETY*h*pow(errmax,PGROW);
  else
    *hnext=5.0*h;
  *x += (*hdid=h);
  for (i=1;i<=n;i++)
    y[i]=ytemp[i];
  free_dvector(ytemp,1,n);
  free_dvector(yerr,1,n);
}

void rkck(double y[], double dydx[], int n, double x, double h, double yout[],
	  double yerr[], void (*derivs)(double, double [], double []))
{
  int i;
  static double a2(0.2),a3(0.3),a4(0.6),a5(1.0),a6(0.875);
  static double b21(0.2);
  static double b31(3.0/40.0),b32(9.0/40.0);
  static double b41(0.3),b42(-0.9),b43(1.2);
  static double b51(-11.0/54.0),b52(2.5),b53(-70.0/27.0),b54(35.0/27.0);
  static double b61(1631.0/55296.0),b62(175.0/521.0),b63(575.0/13824.0);
  static double b64(44275.0/110592.0),b65(253.0/4096.0);
  static double c1(37.0/378.0),c3(250.0/621.0),c4(125.0/594.0);
  static double c6(512.0/1771.0);
  static double dc1(c1-2825.0/27648.0),dc3(c3-18575.0/48384.0);
  static double dc4(c4-13525.0/55296.0),dc5(-277.0/14336.0),dc6(c6-0.25);
  double *ak2,*ak3,*ak4,*ak5,*ak6,*ytemp;

  ak2=dvector(1,n);
  ak3=dvector(1,n);
  ak4=dvector(1,n);
  ak5=dvector(1,n);
  ak6=dvector(1,n);
  ytemp=dvector(1,n);
  for (i=1;i<=n;i++)
    ytemp[i]=y[i]+b21*h*dydx[i];
  (*derivs)(x+a2*h,ytemp,ak2);
  for (i=1;i<=n;i++)
    ytemp[i]=y[i]+h*(b31*dydx[i]+b32*ak2[i]);
  (*derivs)(x+a3*h,ytemp,ak3);
  for (i=1;i<=n;i++)
    ytemp[i]=y[i]+h*(b41*dydx[i]+b42*ak2[i]+b43*ak3[i]);
  (*derivs)(x+a4*h,ytemp,ak4);
  for (i=1;i<=n;i++)
    ytemp[i]=y[i]+h*(b51*dydx[i]+b52*ak2[i]+b53*ak3[i]+b54*ak4[i]);
  (*derivs)(x+a5*h,ytemp,ak5);
  for (i=1;i<=n;i++) 
    ytemp[i]=y[i]+h*(b61*dydx[i]+b62*ak2[i]+b63*ak3[i]
		     +b64*ak4[i]+b65*ak5[i]);
  (*derivs)(x+a6*h,ytemp,ak6);
  for (i=1;i<=n;i++)
    yout[i]=y[i]+h*(c1*dydx[i]+c3*ak3[i]+c4*ak4[i]+c6*ak6[i]);
  for (i=1;i<=n;i++)
    yerr[i]=h*(dc1*dydx[i]+dc3*ak3[i]+dc4*ak4[i]+dc5*ak5[i]+dc6*ak6[i]);
  free_dvector(ytemp,1,n);
  free_dvector(ak6,1,n);
  free_dvector(ak5,1,n);
  free_dvector(ak4,1,n);
  free_dvector(ak3,1,n);
  free_dvector(ak2,1,n);
}

void bsstep(driverData self, double y[], double dydx[],int nv, double *xx, 
	    double htry, double eps, double yscal[], double *hdid, 
	    double *hnext, void (*derivs)(double, double [], double []))
{
  void mmid(double y[], double dydx[], int nvar, double xs, double htot, 
	    int nstep, double yout[], 
	    void (*derivs)(double, double [], double[]));
  void pzextr(driverData bsData, int itest, double xest, double yest[], 
	      double yz[], double dy[], int nv);
  int i,iq,k,kk,km(0);
  double eps1,errmax(TINY),fact,h,red(SAFE2),scale(1.0),work,wrkmin,xest;
  double *err,*yerr,*ysav,*yseq;
  static int nseq[IMAXX+1]={0,2,4,6,8,10,12,14,16,18};
  int reduct,exitflag=0;

  self.d = dmatrix(1,nv,1,KMAXX);
  err = dvector(1,KMAXX);
  self.x = dvector(1,KMAXX);
  yerr = dvector(1,nv);
  ysav = dvector(1,nv);
  yseq = dvector(1,nv);

  if (eps != self.epsold) {
    *hnext = self.xnew = -1.0e29;
    eps1 = SAFE1*eps;
    self.a[1] = nseq[1] + 1;
    for (k=1;k<=KMAXX;k++)
      self.a[k+1] = self.a[k] + nseq[k+1];
    for (iq=2;iq<=KMAXX;iq++) {
      for (k=1;k<iq;k++) {
	self.alf[k][iq] = pow(eps1,(self.a[k+1]-self.a[iq+1])/
			      ((self.a[iq+1]-self.a[1]+1.0)*(2*k+1)));
      }
    }
    self.epsold = eps;
    for (self.kopt=2;self.kopt<KMAXX;self.kopt++) {
      if (self.a[self.kopt+1] > 
	  self.a[self.kopt]*self.alf[self.kopt-1][self.kopt])
	break;
    }
    self.kmax = self.kopt;
  }
  h = htry;
  //  cout << "h=" << h << endl;
  for (i=1;i<=nv;i++)
    ysav[i] = y[i];
  if (*xx != self.xnew || h != (*hnext)) {
    self.first = 1;
    self.kopt = self.kmax;
  }
  reduct = 0;
  for (;;) {
    for (k=1;k<=self.kmax;k++) {
      self.xnew = (*xx) + h;
      if (self.xnew == (*xx)) {
	free_dvector(yseq,1,nv);
	free_dvector(ysav,1,nv);
	free_dvector(yerr,1,nv);
	free_dvector(self.x,1,KMAXX);
	free_dvector(err,1,KMAXX);
	free_dmatrix(self.d,1,nv,1,KMAXX);
	cout << "Underflow occurred!" << endl;
	throw NumRecException("Underflow occurred in bsstep");
      }
      mmid(ysav,dydx,nv,*xx,h,nseq[k],yseq,derivs);
      xest = (h/nseq[k])*(h/nseq[k]);
      pzextr(self,k,xest,yseq,y,yerr,nv);
      if (k != 1) {
	errmax = TINY;
	for (i=1;i<=nv;i++)
	  errmax=fmax(errmax,fabs(yerr[i]/yscal[i]));
	errmax /= eps;
	km = k - 1;
	err[km] = pow(errmax/SAFE1,1.0/(2*km+1));
      }
      //      cout << "errmax=" << errmax << endl;
      if (k != 1 && (k >= self.kopt-1 || self.first)) {
	if (errmax < 1.0) {
	  exitflag = 1;
	  break;
	}
	if (k == self.kmax || k == self.kopt+1) {
	  red = SAFE2/err[km];
	  break;
	}
	else if (k == self.kopt && 
		 self.alf[self.kopt-1][self.kopt] < err[km]) {
	  red = 1.0/err[km];
	  break;
	}
	else if (self.kopt == self.kmax && 
		 self.alf[km][self.kmax-1] < err[km]) {
	  red = self.alf[km][self.kmax-1]*SAFE2/err[km];
	  break;
	}
	else if (self.alf[km][self.kopt] < err[km]) {
	  red = self.alf[km][self.kopt-1]/err[km];
	  break;
	}
      }
    }
    if (exitflag)
      break;
    red = fmin(red,REDMIN);
    red = fmax(red,REDMAX);
    h *= red;
    reduct = 1;
  }
  *xx = self.xnew;
  *hdid = h;
  self.first = 0;
  wrkmin = 1.0e35;
  for (kk=1;kk<=km;kk++) {
    fact = fmax(err[kk],SCALMX);
    work = fact*self.a[kk+1];
    if (work < wrkmin) {
      scale = fact;
      wrkmin = work;
      self.kopt = kk + 1;
    }
  }
  *hnext = h/scale;
  if (self.kopt >= k && self.kopt != self.kmax && !reduct) {
    fact = fmax(scale/self.alf[self.kopt-1][self.kopt],SCALMX);
    if (self.a[self.kopt+1]*fact <= wrkmin) {
      *hnext = h/fact;
      self.kopt++;
    }
  }
  free_dvector(yseq,1,nv);
  free_dvector(ysav,1,nv);
  free_dvector(yerr,1,nv);
  free_dvector(self.x,1,KMAXX);
  free_dvector(err,1,KMAXX);
  free_dmatrix(self.d,1,nv,1,KMAXX);
}

void mmid(double y[], double dydx[], int nvar, double xs, double htot, 
	  int nstep, double yout[], 
	  void (*derivs)(double, double [], double[]))
{
  int n,i;
  double x,swap,h2,h,*ym,*yn;

  ym = dvector(1,nvar);
  yn = dvector(1,nvar);
  h = htot/nstep;
  for (i=1;i<=nvar;i++) {
    ym[i] = y[i];
    yn[i] = y[i] + h*dydx[i];
  }
  x = xs + h;
  (*derivs)(x,yn,yout);
  h2 = 2.0*h;
  for (n=2;n<=nstep;n++) {
    for (i=1;i<=nvar;i++) {
      swap = ym[i] + h2*yout[i];
      ym[i] = yn[i];
      yn[i] = swap;
    }
    x += h;
    (*derivs)(x,yn,yout);
  }
  for (i=1;i<=nvar;i++)
    yout[i] = 0.5*(ym[i]+yn[i]+h*yout[i]);
  free_dvector(yn,1,nvar);
  free_dvector(ym,1,nvar);
}

void pzextr(driverData self, int iest, double xest, double yest[], 
	    double yz[], double dy[], int nv)
{
  int k1,j;
  double q,f2,f1,delta,*c;

  c = dvector(1,nv);
  self.x[iest] = xest;
  for (j=1;j<=nv;j++)
    dy[j] = yz[j] = yest[j];
  if (iest == 1) {
    for (j=1;j<=nv;j++)
      self.d[j][1] = yest[j];
  } else {
    for (j=1;j<=nv;j++)
      c[j] = yest[j];
    for (k1=1;k1<iest;k1++) {
      delta = 1.0/(self.x[iest-k1]-xest);
      f1 = xest*delta;
      f2 = self.x[iest-k1]*delta;
      for (j=1;j<=nv;j++) {
	q = self.d[j][k1];
	self.d[j][k1] = dy[j];
	delta = c[j] - q;
	dy[j] = f1*delta;
	c[j] = f2*delta;
	yz[j] += dy[j];
      }
    }
    for (j=1;j<=nv;j++)
      self.d[j][iest]=dy[j];
  }
  free_dvector(c,1,nv);
}

/* Constructs spline table for 1-D interpolation.  x is dependent variable
 * vector, y is independent variable vector, and table is returned
 * as y2.  n is the array size.  yp1, yp2 are set to 1.0e30 for natural
 * splines. */
void spline(double x[], double y[], int n, double yp1, double ypn,
	    double y2[])
{
  int i,k;
  double p,qn,sig,un,*u;

  u = dvector(1,n-1);
  if (yp1 > 0.99e30)
    y2[1] = u[1] = 0.0;
  else {
    y2[1] = -0.5;
    u[1] = (3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
  }
  for (i=2;i<=n-1;i++) {
    sig = (x[i]-x[i-1])/(x[i+1]-x[i-1]);
    p = sig*y2[i-1] + 2.0;
    y2[i] = (sig-1.0)/p;
    u[i] = (y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
    u[i] = (6.0*u[i]/(x[i+1]-x[i-1]) - sig*u[i-1])/p;
  }
  if (ypn > 0.99e30)
    qn = un = 0.0;
  else {
    qn = 0.5;
    un = (3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
  }
  y2[n] = (un-qn*u[n-1])/(qn*y2[n-1]+1.0);
  for (k=n-1;k>=1;k--)
    y2[k] = y2[k]*y2[k+1]+u[k];
  free_dvector(u,1,n-1);
  return;
}

/* Spline evaluation function.  xa is the dependent variable vector
 * and ya is the independent variable vector.  y2a is the spline table
 * from a previous call to spline.  n is the array size.  x is the point
 * at which the interpolation is to be performed.  Result is returned
 * in y. NOTE xa must be in increasing order! */
void splint(double xa[], double ya[], double y2a[], int n, double x,
	    double *y)
{
  int klo,khi,k;
  double h,b,a;

  klo = 1;
  khi = n;
  while (khi-klo > 1) {
    k = (khi + klo) >> 1;
    if (xa[k] > x)
      khi = k;
    else
      klo = k;
  }
  h = xa[khi] - xa[klo];
  if (h == 0.0)
    throw NumRecException("Bad xa input into routine splint");
  a = (xa[khi] - x)/h;
  b = (x - xa[klo])/h;
  *y = a*ya[klo] + b*ya[khi] + ((a*a*a-a)*y2a[klo] + 
				(b*b*b-b)*y2a[khi])*(h*h)/6.0;
  return;
}

////////////////////////////////////////////////////////////////
// Linear algebra
////////////////////////////////////////////////////////////////
//#include <math.h>
//#define NRANSI
//#include "nrutil.h"
//#define TINY 1.0e-20;

void ludcmp(double **a, int n, int *indx, double *d)
{
	int i,imax,j,k;
	double big,dum,sum,temp;
	double *vv;

	vv=dvector(1,n);
	*d=1.0;
	for (i=1;i<=n;i++) {
		big=0.0;
		for (j=1;j<=n;j++)
			if ((temp=fabs(a[i][j])) > big) big=temp;
		if (big == 0.0) nrerror("Singular matrix in routine ludcmp");
		vv[i]=1.0/big;
	}
	for (j=1;j<=n;j++) {
		for (i=1;i<j;i++) {
			sum=a[i][j];
			for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
		}
		big=0.0;
		for (i=j;i<=n;i++) {
			sum=a[i][j];
			for (k=1;k<j;k++)
				sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
			if ( (dum=vv[i]*fabs(sum)) >= big) {
				big=dum;
				imax=i;
			}
		}
		if (j != imax) {
			for (k=1;k<=n;k++) {
				dum=a[imax][k];
				a[imax][k]=a[j][k];
				a[j][k]=dum;
			}
			*d = -(*d);
			vv[imax]=vv[j];
		}
		indx[j]=imax;
		if (a[j][j] == 0.0) a[j][j]=TINY3;
		if (j != n) {
			dum=1.0/(a[j][j]);
			for (i=j+1;i<=n;i++) a[i][j] *= dum;
		}
	}
	free_dvector(vv,1,n);
}

void lubksb(double **a, int n, int *indx, double b[])
{
	int i,ii=0,ip,j;
	double sum;

	for (i=1;i<=n;i++) {
		ip=indx[i];
		sum=b[ip];
		b[ip]=b[i];
		if (ii){
			for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
		}else if (sum) {
			ii=i;
		}

		b[i]=sum;
	}
	for (i=n;i>=1;i--) {
		sum=b[i];
		for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
		b[i]=sum/a[i][i];
	}
}

void svdcmp(double **a, int m, int n, double w[], double **v)
{
	double pythag(double a, double b);
	int flag,i,its,j,jj,k,l,nm;
	double anorm,c,f,g,h,s,scale,x,y,z,*rv1;

	rv1=dvector(1,n);
	g=scale=anorm=0.0;
	for (i=1;i<=n;i++) {
		l=i+1;
		rv1[i]=scale*g;
		g=s=scale=0.0;
		if (i <= m) {
			for (k=i;k<=m;k++) scale += fabs(a[k][i]);
			if (scale) {
				for (k=i;k<=m;k++) {
					a[k][i] /= scale;
					s += a[k][i]*a[k][i];
				}
				f=a[i][i];
				g = -sign(sqrt(s),f);
				h=f*g-s;
				a[i][i]=f-g;
				for (j=l;j<=n;j++) {
					for (s=0.0,k=i;k<=m;k++) s += a[k][i]*a[k][j];
					f=s/h;
					for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
				}
				for (k=i;k<=m;k++) a[k][i] *= scale;
			}
		}
		w[i]=scale *g;
		g=s=scale=0.0;
		if (i <= m && i != n) {
			for (k=l;k<=n;k++) scale += fabs(a[i][k]);
			if (scale) {
				for (k=l;k<=n;k++) {
					a[i][k] /= scale;
					s += a[i][k]*a[i][k];
				}
				f=a[i][l];
				g = -sign(sqrt(s),f);
				h=f*g-s;
				a[i][l]=f-g;
				for (k=l;k<=n;k++) rv1[k]=a[i][k]/h;
				for (j=l;j<=m;j++) {
					for (s=0.0,k=l;k<=n;k++) s += a[j][k]*a[i][k];
					for (k=l;k<=n;k++) a[j][k] += s*rv1[k];
				}
				for (k=l;k<=n;k++) a[i][k] *= scale;
			}
		}
		anorm=fmax(anorm,(fabs(w[i])+fabs(rv1[i])));
	}
	for (i=n;i>=1;i--) {
		if (i < n) {
			if (g) {
				for (j=l;j<=n;j++)
					v[j][i]=(a[i][j]/a[i][l])/g;
				for (j=l;j<=n;j++) {
					for (s=0.0,k=l;k<=n;k++) s += a[i][k]*v[k][j];
					for (k=l;k<=n;k++) v[k][j] += s*v[k][i];
				}
			}
			for (j=l;j<=n;j++) v[i][j]=v[j][i]=0.0;
		}
		v[i][i]=1.0;
		g=rv1[i];
		l=i;
	}
	for (i=imin(m,n);i>=1;i--) {
		l=i+1;
		g=w[i];
		for (j=l;j<=n;j++) a[i][j]=0.0;
		if (g) {
			g=1.0/g;
			for (j=l;j<=n;j++) {
				for (s=0.0,k=l;k<=m;k++) s += a[k][i]*a[k][j];
				f=(s/a[i][i])*g;
				for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
			}
			for (j=i;j<=m;j++) a[j][i] *= g;
		} else for (j=i;j<=m;j++) a[j][i]=0.0;
		++a[i][i];
	}
	for (k=n;k>=1;k--) {
		for (its=1;its<=30;its++) {
			flag=1;
			for (l=k;l>=1;l--) {
				nm=l-1;
				if ((double)(fabs(rv1[l])+anorm) == anorm) {
					flag=0;
					break;
				}
				if ((double)(fabs(w[nm])+anorm) == anorm) break;
			}
			if (flag) {
				c=0.0;
				s=1.0;
				for (i=l;i<=k;i++) {
					f=s*rv1[i];
					rv1[i]=c*rv1[i];
					if ((double)(fabs(f)+anorm) == anorm) break;
					g=w[i];
					h=pythag(f,g);
					w[i]=h;
					h=1.0/h;
					c=g*h;
					s = -f*h;
					for (j=1;j<=m;j++) {
						y=a[j][nm];
						z=a[j][i];
						a[j][nm]=y*c+z*s;
						a[j][i]=z*c-y*s;
					}
				}
			}
			z=w[k];
			if (l == k) {
				if (z < 0.0) {
					w[k] = -z;
					for (j=1;j<=n;j++) v[j][k] = -v[j][k];
				}
				break;
			}
			if (its == 30) nrerror("no convergence in 30 svdcmp iterations");
			x=w[l];
			nm=k-1;
			y=w[nm];
			g=rv1[nm];
			h=rv1[k];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g=pythag(f,1.0);
			f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x;
			c=s=1.0;
			for (j=l;j<=nm;j++) {
				i=j+1;
				g=rv1[i];
				y=w[i];
				h=s*g;
				g=c*g;
				z=pythag(f,h);
				rv1[j]=z;
				c=f/z;
				s=h/z;
				f=x*c+g*s;
				g = g*c-x*s;
				h=y*s;
				y *= c;
				for (jj=1;jj<=n;jj++) {
					x=v[jj][j];
					z=v[jj][i];
					v[jj][j]=x*c+z*s;
					v[jj][i]=z*c-x*s;
				}
				z=pythag(f,h);
				w[j]=z;
				if (z) {
					z=1.0/z;
					c=f*z;
					s=h*z;
				}
				f=c*g+s*y;
				x=c*y-s*g;
				for (jj=1;jj<=m;jj++) {
					y=a[jj][j];
					z=a[jj][i];
					a[jj][j]=y*c+z*s;
					a[jj][i]=z*c-y*s;
				}
			}
			rv1[l]=0.0;
			rv1[k]=f;
			w[k]=x;
		}
	}
	free_dvector(rv1,1,n);
}

void svbksb(double **u, double w[], double **v, int m, int n, double b[], double x[])
{
	int jj,j,i;
	double s,*tmp;

	tmp=dvector(1,n);
	for (j=1;j<=n;j++) {
		s=0.0;
		if (w[j]) {
			for (i=1;i<=m;i++) s += u[i][j]*b[i];
			s /= w[j];
		}
		tmp[j]=s;
	}
	for (j=1;j<=n;j++) {
		s=0.0;
		for (jj=1;jj<=n;jj++) s += v[j][jj]*tmp[jj];
		x[j]=s;
	}
	free_dvector(tmp,1,n);
}

double pythag(double a, double b)
{
	double absa,absb;
	absa=fabs(a);
	absb=fabs(b);
	if (absa > absb) return absa*sqrt(1.0+SQR(absb/absa));
	else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
}
