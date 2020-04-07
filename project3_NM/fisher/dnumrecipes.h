/* numrecipes.h
 * Header file containing prototypes for functions used from 
 * Numerical Recipes
 *  Originally obtained from Steve Furlanetto in 2004
 */

#ifndef NUM_RECIPES_H
#define NUM_RECIPES_H

// Exception handler for Numerical Recipes methods
class NumRecException {
public:
  NumRecException(char *errorMessage);
  char *what() const { return message; }
private:
  char *message;
};

const double PI(3.14159265359);

// Utility functions
#define NR_END 1 // used by vector
#define FREE_ARG char*

//float *vector(long nl, long nh);
//void free_vector(float *v, long nl, long nh);
float **matrix(long nrl, long nrh, long ncl, long nch);
void free_matrix(float **m, long nrl, long nrh, long ncl, long nch);
double *dvector(long nl, long nh);
void free_dvector(double *v, long nl, long nh);
double **dmatrix(long nrl, long nrh, long ncl, long nch);
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);
int *ivector(long nl, long nh);
void free_ivector(int *v, long nl, long nh);
void nrerror(char error_text[]);

template<class T>
T fmax(T maxarg1, T maxarg2)
{
  return ((maxarg1 > maxarg2) ? maxarg1 : maxarg2);
}
template<class T>
T fmin(T minarg1, T minarg2)
{
  return ((minarg1 > minarg2) ? minarg2 : minarg1);
}

template<class T>
T sign(T a, T b)
{
  return ((b > 0) ? fabs(a) : -fabs(a));
}

template<class T>
T imin(T iminarg1, T iminarg2)
{
 return ((iminarg1 < iminarg2) ? iminarg1 : iminarg2);
}

template<class T>
T SQR(T a)
{
 return (a == 0.0 ? 0.0 : a*a);
}

// Integration functions
const int JMAX(30);   // maximum number of steps for standard
const int JMAXM(14);  // maximum number of steps for open
const int MAXJ(5);    // order of interpolation
const double EPS(1.0e-4);

double trapzd(double (*func)(double), double a, double b, int n, 
	      double *s);
double qromb(double (*func)(double), double a, double b, double acc=EPS);
double qrombnr(double (*func)(double), double a, double b, double acc=EPS);
double qsimp(double (*func)(double), double a, double b,
	     double acc=EPS);
double qromo(double (*func)(double), double a, double b, 
	     double (*choose)(double (*)(double), double, double, int,
			      double (*)),
	     double acc=EPS);
double qsimpo(double (*func)(double), double a, double b, 
	      double acc=EPS);
double midpnt(double (*func)(double), double a, double b, int n,
	      double *s);
double midsql(double (*func)(double), double aa, double bb, int n,
	      double *s);
double midsqu(double (*func)(double), double aa, double bb, int n,
	      double *s);

// Numerical Derivative functions
const int NTAB(10);
const double CON(1.4);

double dfridr(double (*func)(double), double x, double h, double *err);

// Interpolation functions
const int ORDER(5);

void locate(double xx[], int n, double x, int &j);
void polint(double xa[], double ya[], int n, double x, double *y, 
	    double *dy);

//Inversion functions
const int MAXIT(10000);
const double FACTOR(1.6);

int zbrac(void (*func)(const double, double *, double *, int), double &x1, 
	  double &x2);
double rtsafe(void (*funcd)(const double, double *, double *, int), 
	      double x1, double x2, double xacc);
double zriddr(void (*funcd)(const double, double *, double *, int), 
	      double x1, double x2, double xacc);
int zbracSimp(double (*func)(const double), double &x1, double &x2);
double zriddrSimp(double (*func)(const double), double x1, double x2, 
		  double xacc);
double zriddrConst(double (*func)(const double), double yval, double x1, 
		   double x2, double xacc);

// Minimization functions
double brent(double ax, double bx, double cx, double (*f)(double),
	     double tol, double *xmin);

// Elliptic Integrals and Special Functions

const double ERRTOL(0.0015);
const double BIG(4.5e21);
const double TINY2(1.0e-25);

double rf(const double x, const double y, const double z);
double ellf(double phi, double ak);
double ellfsq(double phi, double ak2);
double rd(const double x, const double y, const double z);
double elle(double phi, double ak);

double bessk(int n, double x);
double bessj0(double x);
double bessi0(double x);
double bessk0(double x);
double bessk1(double x);
double bessi1(double x);

double erfcc(double x);
double erff(double x);

#define EULER 0.57721566

void cisi(double x, double *ci, double *si);

const double TMIN(2.0);
const int TRU(1);

// Random number generators

float ran1(long *idum);
double gasdev(long *idum);

// ODE section
const double SAFETY(0.9);
const double PGROW(-0.2);
const double PSHRNK(-0.25);
const double ERRCON(1.89e-4);

const int KMAXX(8);
const int IMAXX(KMAXX+1);
const double SAFE1(0.25);
const double SAFE2(0.7);
const double REDMAX(1.0e-5);
const double REDMIN(0.7);
const double SCALMX(0.1);

const int MAXSTP(10000);
const double ACCURACY(1.0e-6);
const double TINY(1.0e-30);

// "Private" data members for bsstep algorithm (used so cascading works)
struct driverData { 
  double **d,*x;
  double epsold,xnew;
  double a[IMAXX+1];
  double alf[KMAXX+1][KMAXX+1];
  int first,kmax,kopt;
};

void odeint(double ystart[], int nvar, double x1, double x2, double eps, 
	    double h1, double hmin, int *nok, int *nbad,
	    void  (*derivs)(double, double [], double []),
	    void (*rkqs)(driverData, double [], double [], int, double *, 
			 double, double, double [], double *, double *, 
			 void (*)(double, double [], double [])));
void rkqs(driverData self, double y[], double dydx[], int n, double *x, 
	  double htry, double eps, double yscal[], double *hdid, 
	  double *hnext, void (*derivs)(double, double [], double []));
void bsstep(driverData self, double y[], double dydx[],int nv, double *xx, 
	    double htry, double eps, double yscal[], double *hdid, 
	    double *hnext, void (*derivs)(double, double [], double []));


void spline(double x[], double y[], int n, double yp1, double ypn,
	    double y2[]);
void splint(double xa[], double ya[], double y2a[], int n, double x,
	    double *y);

// linear algebra
#define TINY3 1.0e-20;
void ludcmp(double **a, int n, int *indx, double *d);
void lubksb(double **a, int n, int *indx, double b[]);
void svdcmp(double **a, int m, int n, double w[], double **v);
void svbksb(double **u, double w[], double **v, int m, int n, double b[], double x[]);
double pythag(double a, double b);
#endif
