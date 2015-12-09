/*
  ECG for ARDUINO

  ECGSYN: A program for generating a realistic synthetic ECG.
  Copyright (c) 2003 by Patrick McSharry & Gari Clifford. All rights reserved.

  Arduino version by Fernanda Andrade.
  ESPE, 2013
*/


#include <LiquidCrystal.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define PI (2.0*asin(1.0))
#define SWAP(a,b) tempr=a;a=b;b=tempr
#define MIN(a,b) (a < b ? a : b)
#define MAX(a,b) (a > b ? a : b)
#define OFFSET 1
#define ARG1 char*

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

LiquidCrystal lcd(12, 11, 5, 4, 3, 2);

/*--------------------------------------------------------------------------*/
/*    DEFINE PARAMETERS AS GLOBAL VARIABLES                                 */
/*--------------------------------------------------------------------------*/

int N = 1;                   /*  Number of heart beats              */
int sfecg = 64;               /*  ECG sampling frequency             */
int sf = 64;                  /*  Internal sampling frequency        */
double Anoise = 0.0;           /*  Amplitude of additive uniform noise*/
double hrmean = 60.0;          /*  Heart rate mean                    */
double hrstd = 1.0;            /*  Heart rate std                     */
double flo = 0.1;              /*  Low frequency                      */
double fhi = 0.25;             /*  High frequency                     */
double flostd = 0.01;          /*  Low frequency std                  */
double fhistd = 0.01;          /*  High frequency std                 */
double lfhfratio = 0.5;        /*  LF/HF ratio                        */

int Necg = 0;                  /*  Number of ECG outputs              */
int mstate = 3;                /*  System state space dimension       */
double xinitial = 1.0;         /*  Initial x co-ordinate value        */
double yinitial = 0.0;         /*  Initial y co-ordinate value        */
double zinitial = 0.04;        /*  Initial z co-ordinate value        */
int seed = 1;                  /*  Seed                               */
long rseed;
double h;
double *rr,*rrpc;
double *ti,*ai,*bi;

void setup() {
  lcd.begin(16, 2);

  // Run program
  program();
}

void loop() {
  // Empty
}

void dfour1(double data[], int nn, int isign){

  int n,mmax,m,j,istep,i;
  double wtemp,wr,wpr,wpi,wi,theta;
  double tempr,tempi;

	n=nn << 1;
	j=1;
	for (i=1;i<n;i+=2) {
		if (j > i) {
                        SWAP(data[j],data[i]);
			SWAP(data[j+1],data[i+1]);
		}
		m=n>>1;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j+=m;
	}

	mmax=2;
	while (n > mmax) {
		istep=mmax << 1;
		theta=isign*(6.28318530717959/mmax);
		wtemp=sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi=sin(theta);
		wr=1.0;
		wi=0.0;
		for (m=1;m<mmax;m+=2) {
			for (i=m;i<=n;i+=istep) {
				j=i+mmax;
				tempr=wr*data[j]-wi*data[j+1];
				tempi=wr*data[j+1]+wi*data[j];
				data[j]=data[i]-tempr;
				data[j+1]=data[i+1]-tempi;
				data[i] += tempr;
				data[i+1] += tempi;
			}
			wr=(wtemp=wr)*wpr-wi*wpi+wr;
			wi=wi*wpr+wtemp*wpi+wi;
		}
		mmax=istep;
	}

}

float ran1(long *idum){
	int j;
	long k;
	static long iy=0;
	static long iv[NTAB];
	float temp;

	if (*idum <= 0 || !iy) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ;
			*idum=IA*(*idum-k*IQ)-IR*k;
			if (*idum < 0) *idum += IM;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ;
	*idum=IA*(*idum-k*IQ)-IR*k;
	if (*idum < 0) *idum += IM;
	j=iy/NDIV;
	iy=iv[j];
	iv[j] = *idum;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}


/*---------------------------------------------------------------------------*/
/*      ALLOCATE MEMORY FOR VECTOR                                           */
/*---------------------------------------------------------------------------*/

double *mallocVect(long n0, long nx){
        double *vect;

        vect=(double *)malloc((size_t) ((nx-n0+1+OFFSET)*sizeof(double)));
        return vect-n0+OFFSET;
}

/*---------------------------------------------------------------------------*/
/*      FREE MEMORY FOR MALLOCVECT                                           */
/*---------------------------------------------------------------------------*/

void freeVect(double *vect, long n0, long nx){
        free((ARG1) (vect+n0-OFFSET));
}

/*---------------------------------------------------------------------------*/
/*      MEAN CALCULATOR                                                      */
/*---------------------------------------------------------------------------*/

double mean(double *x, int n){
/* n-by-1 vector, calculate mean */
        int j;
        double add;

        add = 0.0;
        for(j=1;j<=n;j++)  add += x[j];

        return (add/n);
}


/*---------------------------------------------------------------------------*/
/*      STANDARD DEVIATION CALCULATOR                                        */
/*---------------------------------------------------------------------------*/

double stdev(double *x, int n){
/* n-by-1 vector, calculate standard deviation */
        int j;
        double add,mean,diff,total;

        add = 0.0;
        for(j=1;j<=n;j++)  add += x[j];
        mean = add/n;

        total = 0.0;
        for(j=1;j<=n;j++)
        {
           diff = x[j] - mean;
           total += diff*diff;
        }

        return (sqrt(total/(n-1)));
}

/*--------------------------------------------------------------------------*/
/*    INTERP                                                                */
/*--------------------------------------------------------------------------*/

void interp(double *y, double *x, int n, int r){
   int i,j;
   double a;

   for(i=1;i<=n-1;i++)
   {
      for(j=1;j<=r;j++)
      {
         a = (j-1)*1.0/r;
         y[(i-1)*r+j] = (1.0-a)*x[i] + a*x[i+1];
      }
   }
}


/*--------------------------------------------------------------------------*/
/*    GENERATE RR PROCESS                                                   */
/*--------------------------------------------------------------------------*/

void rrprocess(double *rr, double flo, double fhi, double flostd, double fhistd, double lfhfratio, double hrmean, double hrstd, double sf, int n){
   int i;
   double xstd,ratio,m;
   double *Sw,*ph,*SwC;

   Sw = mallocVect(1,n);
   ph = mallocVect(1,n);
   SwC = mallocVect(1,2*n);

   for(i=1;i<=n;i++) {

     if(i<=(n/2)){
	Sw[i] = (sf/2.0)*sqrt(lfhfratio*exp(-((i-1)*2.0*PI*(sf/n)-(2.0*PI*flo))*((i-1)*2.0*PI*(sf/n)-(2.0*PI*flo))/(2.0*(2.0*PI*flostd)*(2.0*PI*flostd)))/sqrt(2*PI*(2.0*PI*flostd)*(2.0*PI*flostd)) + 1.0*exp(-((i-1)*2.0*PI*(sf/n)-(2.0*PI*fhi))*((i-1)*2.0*PI*(sf/n)-(2.0*PI*fhi))/(2.0*(2.0*PI*fhistd)*(2.0*PI*fhistd)))/sqrt(2*PI*(2.0*PI*fhistd)*(2.0*PI*fhistd)));

	if(i<=(n/2-1)){
		m=2.0*PI*ran1(&rseed);
		ph[i+1] = m;
		ph[n-i+1] = -m;
	}

      }

      if(i>=(n/2+1)){
	Sw[i] = (sf/2.0)*sqrt(lfhfratio*exp(-((n-i)*2.0*PI*(sf/n)-(2.0*PI*flo))*((n-i)*2.0*PI*(sf/n)-(2.0*PI*flo))/(2.0*(2.0*PI*flostd)*(2.0*PI*flostd)))/sqrt(2*PI*(2.0*PI*flostd)*(2.0*PI*flostd))+ 1.0*exp(-((n-i)*2.0*PI*(sf/n)-(2.0*PI*fhi))*((n-i)*2.0*PI*(sf/n)-(2.0*PI*fhi))/(2.0*(2.0*PI*fhistd)*(2.0*PI*fhistd)))/sqrt(2*PI*(2.0*PI*fhistd)*(2.0*PI*fhistd)));
      }

   }

  /* make complex spectrum */
   for(i=1;i<=n;i++) SwC[2*i-1] = Sw[i]*cos(ph[i]);
   for(i=1;i<=n;i++) SwC[2*i] = Sw[i]*sin(ph[i]);

   /* calculate inverse fft */
   dfour1(SwC,n,-1);

   /* extract real part */
   for(i=1;i<=n;i++) rr[i] = (1.0/n)*SwC[2*i-1];

   xstd = stdev(rr,n);
   ratio = (60.0*hrstd/(hrmean*hrmean))/xstd;

   for(i=1;i<=n;i++) rr[i] *= ratio;
   for(i=1;i<=n;i++) rr[i] += (60.0/hrmean);

   freeVect(Sw,1,n);
   freeVect(ph,1,n);
   freeVect(SwC,1,2*n);
}

/*--------------------------------------------------------------------------*/
/*    THE ANGULAR FREQUENCY                                                 */
/*--------------------------------------------------------------------------*/

double angfreq(double t){
   int i;

   i = 1 + (int)floor(t/h);

   return 2.0*PI/rrpc[i];
}

/*--------------------------------------------------------------------------*/
/*    THE EXACT NONLINEAR DERIVATIVES                                       */
/*--------------------------------------------------------------------------*/

void derivspqrst(double t0,double x[],double dxdt[]){
   int i,k;
   double a0,w0;
   double t,dt,dt2,*xi,*yi,zbase;

   k = 5;
   xi = mallocVect(1,k);
   yi = mallocVect(1,k);

   w0 = angfreq(t0);
   a0 = 1.0 - sqrt((x[1]-0.0)*(x[1]-0.0) + (x[2]-0.0)*(x[2]-0.0))/1.0;;

   for(i=1;i<=k;i++) xi[i] = cos(ti[i]);
   for(i=1;i<=k;i++) yi[i] = sin(ti[i]);

   zbase = 0.005*sin(2.0*PI*fhi*t0);

   t = atan2(x[2],x[1]);
   dxdt[1] = a0*(x[1] - 0.0) - w0*(x[2] - 0.0);
   dxdt[2] = a0*(x[2] - 0.0) + w0*(x[1] - 0.0);
   dxdt[3] = 0.0;
   for(i=1;i<=k;i++) {
      dt = fmod(t-ti[i],2.0*PI);
      dt2 = dt*dt;
      dxdt[3] += -ai[i]*dt*exp(-0.5*dt2/(bi[i]*bi[i]));
   }
   dxdt[3] += -1.0*(x[3] - zbase);

   freeVect(xi,1,k);
   freeVect(yi,1,k);
}

/*--------------------------------------------------------------------------*/
/*    RUNGA-KUTTA FOURTH ORDER INTEGRATION                                  */
/*--------------------------------------------------------------------------*/

void drk4(double y[], int n, double x, double h, double yout[], void (*derivs)(double, double [], double [])){

        int i;
        double *dydx,*dym,*dyt,*yt;

        dydx=mallocVect(1,n);
        dym=mallocVect(1,n);
        dyt=mallocVect(1,n);
        yt=mallocVect(1,n);

        (*derivs)(x,y,dydx);
        for (i=1;i<=n;i++) yt[i]=y[i]+(h*0.5)*dydx[i];

        (*derivs)(x+h*0.5,yt,dyt);

        for (i=1;i<=n;i++) yt[i]=y[i]+(h*0.5)*dyt[i];
        (*derivs)(x+h*0.5,yt,dym);
        for (i=1;i<=n;i++) {
                yt[i]=y[i]+h*dym[i];
                dym[i] += dyt[i];
        }

        (*derivs)(x+h,yt,dyt);
        for (i=1;i<=n;i++)
                yout[i]=y[i]+(h/6.0)*(dydx[i]+dyt[i]+2.0*dym[i]);

        freeVect(dydx,1,n);
        freeVect(dym,1,n);
        freeVect(dyt,1,n);
        freeVect(yt,1,n);

}


/*--------------------------------------------------------------------------*/
/*    MAIN PART OF PROGRAM                                                 */
/*--------------------------------------------------------------------------*/

int program(){
   int i,j,k,q,Nrr,Nt,Nts;
   double tstep,tecg;
   double *x,*xt,*yt,*zt,*xts,*yts,*zts;
   double timev,zmin,zmax,zrange;
      void (*derivs)(double, double [], double []);

   /* perform some checks on input values */
   q = (int)round(sf/sfecg);

   /* declare and initialise the state vector */
   x=mallocVect(1,mstate);
   x[1] = xinitial;
   x[2] = yinitial;
   x[3] = zinitial;

   /* declare and define the ECG morphology vectors (PQRST extrema parameters) */
   ti=mallocVect(1,5);
   ai=mallocVect(1,5);
   bi=mallocVect(1,5);

   /* P            Q            R           S           T        */

   ti[1]=-60.0; ti[2]=-15.0; ti[3]=0.0;  ti[4]=15.0; ti[5]=90.0;
   ai[1]=1.2;   ai[2]=-5.0;  ai[3]=30.0; ai[4]=-7.5; ai[5]=0.75;
   bi[1]=0.25;  bi[2]=0.1;   bi[3]=0.1;  bi[4]=0.1;  bi[5]=0.4;

   /* convert angles from degrees to radians */
   for(i=1;i<=5;i++) ti[i] *= PI/180.0;

   /* adjust extrema parameters for mean heart rate */
   for(i=1;i<=5;i++) bi[i] *= sqrt(hrmean/60.0);
   ti[1]*=sqrt(sqrt(hrmean/60.0));  ti[2]*=sqrt(hrmean/60.0); ti[3]*=1.0; ti[4]*=sqrt(hrmean/60.0); ti[5]*=1.0;

   /* calculate time scales */
   h = 1.0/sf;
   tstep = 1.0/sfecg;

   /* initialise seed */
   rseed = -seed;

   /* select the derivs to use */
   derivs = derivspqrst;

   /* calculate length of RR time series */
   Nrr = (int)pow(2.0, ceil(log10(N*(60/hrmean)*sf)/log10(2.0)))+1;

   /* create rrprocess with required spectrum */
   rr = mallocVect(1,Nrr);
   rrprocess(rr, flo, fhi, flostd, fhistd, lfhfratio, hrmean, hrstd, sf, Nrr);

   /* create piecewise constant rr */
   rrpc = mallocVect(1,2*Nrr);
   tecg = 0.0;
   i = 1;
   j = 1;
   while(i <= Nrr){
      tecg += rr[j];
      j = (int)round(tecg/h);
      for(k=i;k<=j;k++) rrpc[k] = rr[i];
      i = j+1;
   }
   Nt = j;

     /* integrate dynamical system using fourth order Runge-Kutta*/
   zt = mallocVect(1,Nt);

   timev = 0.0;
   for(i=1;i<=Nt;i++){
      zt[i] = x[3];
      drk4(x, mstate, timev, h, x, derivs);
      timev += h;
   }

   /* downsample to ECG sampling frequency */
   zts = mallocVect(1,Nt);

   j=0;
   for(i=1;i<=Nt;i+=q){
      j++;
      zts[j] = zt[i];
   }
   Nts = j;

   /* scale signal to lie between -0.4 and 1.2 mV */
   zmin = zts[1];
   zmax = zts[1];
   for(i=2;i<=Nts;i++){
     if(zts[i] < zmin)       zmin = zts[i];
     else if(zts[i] > zmax)  zmax = zts[i];
   }

   zrange = zmax-zmin;

   /* output ECG */
   for(i=1;i<=Nts;i++){
        zts[i] = (zts[i]-zmin)*(1.6)/zrange - 0.4;
        zts[i] += Anoise*(2.0*ran1(&rseed) - 1.0);

  /* Display in LCD */
	lcd.setCursor(0, 1);
	lcd.print(zts[i]);
	delay(1000);
  }

freeVect(x,1,mstate);
freeVect(rr,1,Nrr);
freeVect(rrpc,1,2*Nrr);
freeVect(ti,1,5);
freeVect(ai,1,5);
freeVect(bi,1,5);
freeVect(zt,1,Nt);
freeVect(zts,1,Nt);

/* END OF DORUN */
}
