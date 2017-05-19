/*code from Numerical Recipes in C*/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "nr.h"
#include "spmain.h"
#define PI 3.141592654

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

float ran1(long *idum)
{
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

#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX



/* (C) Copr. 1986-92 Numerical Recipes Software . */
long int timeseed(void)
/*output a seed for random generator based on system time by  multiplying sec*min*hours or give 1 if 
sec=min=hour=0*/
{
		struct tm *local;
		time_t t;
		long int seed;

		seed=1;
		t = time(NULL);
		local = localtime(&t); 
		if(local->tm_sec) seed=local->tm_sec; 
		if(local->tm_min) seed *= local->tm_min; 
		if(local->tm_hour) seed *= local->tm_hour;
		return seed;
}

float poidev(float xm, long *idum)
/*random deviate from a Poisson distribution of mean xm using ran1(idum) as a source 
of uniform random deviates. from Numerical Recipes*/
{
	float gammln(float xx);
	float ran1(long *idum);
	static float sq,alxm,g,oldm=(-1.0); 
	float em,t,y;
	if (xm < 12.0) {
		if (xm != oldm) { 
			oldm=xm; 
			g=exp(-xm);
		}
		em = -1;
		t=1.0;
		do {
		
			++em;
			t *= ran1(idum); 
		} while (t > g);
	} else {
		if (xm != oldm) { 
			oldm=xm; 
			sq=sqrt(2.0*xm); 
			alxm=log(xm); 
			g=xm*alxm-gammln(xm+1.0);
		}
		do {
			do {
				y=tan(PI*ran1(idum));
				em=sq*y+xm;
			} while (em < 0.0);
			em=floor(em); 
			t=0.9*(1.0+y*y)*exp(em*alxm-gammln(em+1.0)-g);
		} while (ran1(idum) > t);
	}
	return em;
}
#undef PI
/* (C) Copr. 1986-92 Numerical Recipes Software . */
/*******************************************************/
float gammln(float xx)
/*numerical recipes, returns the value of ln(gamma(xx)) for xx>0*/
{
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677, 
		24.01409824083091,-1.231739572450155, 
		0.1208650973866179e-2,-0.5395239384953e-5};
	int j;
	
	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y; 
	return -tmp+log(2.5066282746310005*ser/x);
}
/* (C) Copr. 1986-92 Numerical Recipes Software . */

float bico(int n, int k)
/*returns the binomial coefficient Cn(haut)k(bas) = n!/(k!(n-k)!) as a floating-
point number*/
{
	float factln(int n);
	return floor(0.5+exp(factln(n)-factln(k)-factln(n-k))); 
}

float factln(int n)
/*returns ln(n!)*/
{
	float gammln(float xx);
	static float a[101];
	if (n<0) {
		printf("Negative factorial in routine factln"); 
		exit(1);
	}
	if (n <=1) return 0.0;
	if (n <= 100) return a[n] ? a[n] : (a[n]=gammln(n+1.0));
	else return gammln(n+1.0);
}
/************************************************************************/
/* GAMMALN - natural log of Gamma(xx) */
double gammaln(double xx)
{
	double x, tmp, ser;
	static double cof[6] = {76.18009173, -86.50532033, 24.01409822,
				-1.231739516, 0.120858003e-2, -0.536382e-5}; 
	int j;
	
	x = xx - 1.0;
	tmp = x + 5.5;
	tmp -= (x+0.5)*log(tmp);
	ser = 1.0;
	for(j=0; j<=5; j++) {
		x += 1.0;
		ser += cof[j]/x;
		}
	return -tmp+log(2.50662827465*ser);
} 
/*************************************************************************/
#define PI 3.141592654
/* POISSON - returns Poisson(xm) deviates */
double poisson(double xm, long int *idum)
{
	static double sq, alxm, g, oldm=(-1.0);
	double em, t, y;
	if(xm < 12.0) {
		if(xm != oldm) {
			oldm = xm;
			g = exp(-xm);
			}
		em = -1;
		t = 1.0;
		do {
			em += 1.0;
			t *= ran1(idum); 
		} while(t > g);
	} else {
		if(xm != oldm) {
			oldm = xm;
			sq = sqrt(2.0*xm);
			alxm = log(xm);
			g = xm*alxm - gammln(xm+1.0);
		}
		do {
			do {
				y = tan(PI*(ran1(idum)));
				em = sq*y + xm;
			} while(em < 0.0);
			em = floor(em);
			t = 0.9*(1.0 + y*y)*exp(em*alxm - gammln(em+1.0) - g); 
		} while((ran1(idum)) > t);
	}
	return em;
}
#undef PI

/************************************************************************************/
#define FUNC(x) ((*func) (x))

#define EPS 1.0e-6	        /*fractional accuracy desired*/
#define JMAX 30		        /*limits the total number of steps*/
#define JMAXP (JMAX+1)
#define K 5		                /*number of points used in the extrapolation*/

float qromb(double (*func)(float),float a, float b)
/*returns the integral of the function func from a to b. Integration is performed
by Romberg's method of order 2K, where K=2 is Simpson's rule, Numerical Recipes p.140*/
{
		float ss,dss;
		double s[JMAXP+1],h[JMAXP+1];
		int j;
		h[1]=1.0;
		for (j=1;j<=JMAX;j++) {
			s[j]=trapzd(func,a,b,j);
			if (j >= K) {
				polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
				if (fabs(dss) < EPS*fabs(ss)) return ss;
			}
			s[j+1]=s[j];
			h[j+1]=0.25*h[j];
		}
printf ("\ntoo many steps in integration routine qromb");
 return 9999;
		/*exit(1);*/
}
#undef EPS
#undef JMAX
#undef JMAXP
#undef K

/*********************************************************************/
double trapzd(double (*func)(float), float a, float b, int n)
/*extended trapezoidal rule for integration, Numerical recipes pl 137*/ 
{
		float x,tnm,del;
		double sum;
		static double s;
/*		int it,j;  */
		long it,j;
		if (n == 1) {
			return (s=0.5*(b-a)*( ((*func)(a)) + ((*func)(b)) ));
		} else {
			for (it=1,j=1;j<n-1;j++) it <<= 1;
			tnm=it;
			del=(b-a)/tnm;
			x=a+0.5*del;
			for (sum=0.0,j=1;j<=it;j++,x+=del) sum +=  ((*func)(x)); 
			s=0.5*(s+(b-a)*sum/tnm);
			return s;
		}
}
/*******************************************************************
**/
void polint(double xa[],double ya[],int n, float x,float *y, float *dy) /*polynomial 
interpolation from Numerical recipes p.109*/
{
	int i,m,ns=1;
	double den,dif,dift,ho,hp,w;
	float *c, *d;
	dif=fabs(x-xa[1]);
	c=vector(1,n);
	d=vector(1,n);
	for (i=1;i<=n;i++) {
		if ( (dift=fabs(x-xa[i])) < dif) 
		{ 
			ns=i;
			dif=dift;
		}
		c[i]=ya[i];
		d[i]=ya[i];
	}
	*y=ya[ns--];
	for (m=1;m<n;m++) {
		for (i=1;i<=n-m;i++) { 
			ho=xa[i]-x; 
			hp=xa[i+m]-x; 
			w=c[i+1]-d[i];
			if ( (den=ho-hp) == 0.0) {
				printf("\nerror in routine polint"); 
				exit(1);
			}
			den=w/den;
			d[i]=hp*den;
			c[i]=ho*den;
		}
		*y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
	}
	free_vector(d,1,n);
	free_vector(c,1,n);
}
/*******************************************************************
**/
#define FACTOR	1.6
#define NTRY 50
int zbrac(double (*func)(float), float *x1, float *x2)
{
	int j;
	double f1,f2;
	if (*x1 == *x2) {
		printf("Bad initial range in zbrac"); exit(1);
	}
	f1=(*func)(*x1);
	f2=(*func)(*x2);
	for (j=1;j<=NTRY;j++) {
		if (f1*f2 < 0.0) return 1;
		if (fabs(f1) < fabs(f2)) f1=(*func)(*x1 += FACTOR*(*x1-*x2));
		else f2=(*func)(*x2 += FACTOR*(*x2-*x1));
	}
	return 0;
}
#undef FACTOR
#undef NTRY
/* (C) Copr. 1986-92 Numerical Recipes Software . */
/*******************************************************************
**/
#define NRANSI
#define ITMAX 100
#define EPS 3.0e-8
float zbrent(double (*fx) (float),float x1, float x2, float tol)
{
	int iter;
	float a=x1,b=x2,c=x2,d,e,min1,min2;
	double fa=(*fx)(a),fb=(*fx)(b),fc,p,q,r,s,tol1,xm;
	if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
	{ 
		write(ERRORFILE,"\nError in zbrent: Root must be bracketed in zbrent"); 
		return (a+b/2.0);
		/*exit(1);*/
	}
	fc=fb;
	for (iter=1;iter<=ITMAX;iter++) {
		if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) { c=a;
			fc=fa;
			e=d=b-a;
		}
		if (fabs(fc) < fabs(fb)) {
			a=b;
			b=c;
			c=a;
			fa=fb;
			fb=fc;
			fc=fa;
		}
		tol1=2.0*EPS*fabs(b)+0.5*tol;
		xm=0.5*(c-b);
		if (fabs(xm) <= tol1 || fb == 0.0) return b;
		if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) 
		{ 
			s=fb/fa;
			if (a == c) 
			{
				p=2.0*xm*s;
				q=1.0-s;
			} else {
				q=fa/fc;
				r=fb/fc; 
				p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0)); 
				q=(q-1.0)*(r-1.0)*(s-1.0); 
			}
			if (p > 0.0) q = -q;
			p=fabs(p);
			min1=3.0*xm*q-fabs(tol1*q); 
			min2=fabs(e*q);
			if (2.0*p < (min1 < min2 ? min1 : min2)) 
			{ 
				e=d;
				d=p/q;
			} else {
				d=xm;
				e=d;
			}
		} else {
			d=xm;
			e=d;
		}
		a=b;
		fa=fb;
		if (fabs(d) > tol1) b += d;
		else b += (xm > 0.0 ? fabs(tol1) : -fabs(tol1)); 
		fb=(*fx)(b);
	}
printf("Maximum number of iterations exceeded in zbrent"); exit(1);
	/*return 0.0;*/
}
#undef ITMAX
#undef EPS
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software . */
/*******************************************************************
**/
void zbrak(double (*fx) (float), float x1, float x2, int n, float xb1[],float xb2[], int *nb) {
	int nbb,i;
	double x,fp,fc,dx;
	nbb=-1;
	dx=(x2-x1)/n;
	fp=(*fx)(x=x1);
	for (i=0;i<n;i++) {
		fc=(*fx)(x += dx);
		if (fc*fp < 0.0) { 
			xb1[++nbb]=x-dx; 
			xb2[nbb]=x;
			if(*nb-1 == nbb) return; 
		}
		fp=fc;
	}
	*nb = nbb+1;
}
/* (C) Copr. 1986-92 Numerical Recipes Software . */
/*******************************************************************
**/
#define NRANSI
#define ITMAX 100
#define CGOLD 0.3819660
#define ZEPS 1.0e-10
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
float brent(float ax, float bx, float cx,float tol,float *xmin,double (*func)(float)) {
	int iter;
	double a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
	float e=0.0;
	a=(ax < cx ? ax : cx);
	b=(ax > cx ? ax : cx);
	x=w=v=bx;
	fw=fv=fx=(*func)(x);
for (iter=1;iter<=ITMAX;iter++) {
	xm=0.5*(a+b); 
	tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
	if (fabs(x-xm) <= (tol2-0.5*(b-a))) { 
			*xmin=x;
			return fx;
	}
	if (fabs(e) > tol1) 
	{ 
		r=(x-w)*(fx-fv); 
		q=(x-v)*(fx-fw);
			p=(x-v)*q-(x-w)*r;
			q=2.0*(q-r);
			if (q > 0.0) p = -p;
			q=fabs(q);
			etemp=e;
			e=d;
if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x)) d=CGOLD*(e=(x >= xm ? a-x : b-x));
			else {
				d=p/q;
				u=x+d;
				if (u-a < tol2 || b-u < tol2) d=SIGN(tol1,xm-x);
			}
		} else {
			d=CGOLD*(e=(x >= xm ? a-x : b-x));
		}
		u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
		fu=(*func)(u);
		if (fu <= fx) {
			if (u >= x) a=x;
			else b=x;
			SHFT(v,w,x,u)
			SHFT(fv,fw,fx,fu)
		} else {
			if (u < x) a=u; else b=u;
if (fu <= fw || w == x) {
v=w;
w=u;
fv=fw;
fw=fu;
} else if (fu <= fv || v == x || v == w) { v=u;
fv=fu;
			}
		}
	}
printf("Maximum number of iterations exceeded in zbrent"); 
exit(1);
	*xmin=x;
	return fx;
}
#undef ITMAX
#undef CGOLD
#undef ZEPS
#undef SHFT
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software	. */
/*******************************************************************
**/
#define NRANSI
#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY 1.0e-20
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
void mnbrak(float *ax, float *bx, float *cx, float *fa, float *fb, float *fc, double 
(*func)(float))
{
	double ulim,u,r,q,fu,dum;
	*fa=(*func)(*ax);
	*fb=(*func)(*bx);
	if (*fb > *fa) {
		SHFT(dum,*ax,*bx,dum)
		SHFT(dum,*fb,*fa,dum)
	}
*cx=(*bx)+GOLD*(*bx-*ax);
	*fc=(*func)(*cx);
while (*fb > *fc) { 
r=(*bx-*ax)*(*fb-*fc); 
q=(*bx-*cx)*(*fb-*fa); 
u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/(2.0*SIGN(FMAX(fabs(q-r),TINY),q-r)); 
ulim=(*bx)+GLIMIT*(*cx-*bx);
if ((*bx-u)*(u-*cx) > 0.0) { 
fu=(*func)(u);
			if (fu < *fc) {
				*ax=(*bx);
				*bx=u;
				*fa=(*fb);
				*fb=fu;
				return;
			} else if (fu > *fb) {
				*cx=u;
				*fc=fu;
				return;
			}
			u=(*cx)+GOLD*(*cx-*bx);
fu=(*func)(u);
} else if ((*cx-u)*(u-ulim) > 0.0) {
fu=(*func)(u);
if (fu < *fc) { 
SHFT(*bx,*cx,u,*cx+GOLD *(*cx-*bx)) 
SHFT(*fb,*fc,fu,(*func)(u))
			}
		} else if ((u-ulim)*(ulim-*cx) >= 0.0) {
			u=ulim;
fu=(*func)(u);
		} else {
u=(*cx)+GOLD*(*cx-*bx);
			fu=(*func)(u);
		}
		SHFT(*ax,*bx,*cx,u)
		SHFT(*fa,*fb,*fc,fu)
	}
}
#undef GOLD
#undef GLIMIT
#undef TINY
#undef SHFT
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software . */
/*********************************************************************/

float erff(float x)
/*returns the error function erf(x)
uses gammp,gser,gcf and gammln*/
{
	float gammp(float a, float x);
	
	return x < 0.0 ? -gammp(0.5,x*x) : gammp(0.5,x*x);
}
/*********************************************************************/

float gammp(float a, float x)
/* returns the incomplete gamma function P(a,x)*/
{
	float gamser,gammcf,gln;
	
	if (x < 0.0 || a <= 0.0) {
		write("erreur","\n Invalid arguments in routine gammp");
		exit(1);
	}
	if (x < (a+1.0)) {
		gser(&gamser,a,x,&gln);
		return gamser;
	} else {
		gcf(&gammcf,a,x,&gln);
		return 1.0-gammcf;
	}
}
/*********************************************************************/
#define ITMAX 100
#define EPS 3.0e-7

void gser(float *gamser,float a, float x, float *gln)
{
	int n;
	float sum,del,ap;
	
	*gln=gammln(a);
	if (x <= 0.0) {
		if (x<0.0) {
			write("erreur","\n x less than 0 in routine gser");
			exit(1);
		}
		*gamser=0.0;
		return;
	} else {
		ap=a;
		del=sum=1.0/a;
		for(n=1;n<=ITMAX;n++) {
			++ap;
			del *= x/ap;
			sum += del;
			if (fabs(del) < fabs(sum)*EPS) {
				*gamser=sum*exp(-x+a*log(x)-(*gln));
				return;
			}
		}
		write("erreur","\n a too large, ITMAX too small in routine gser");
		exit(1);
	}
}
		
#undef EPS
#undef ITMAX

/*********************************************************************/

#define ITMAX 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30

void gcf(float *gammcf,float a,float x,float *gln)
{
	int i;
	float an,b,c,d,del,h;
	
	*gln=gammln(a);
	b=x+1.0-a;
	c=1.0/FPMIN;
	d=1.0/b;
	h=d;
	for(i=1;i<=ITMAX;i++) {
		an = -i*(i-a);
		b += 2.0;
		d = an*d+b;
		if(fabs(d) < FPMIN) d=FPMIN;
		c=b+an/c;
		if(fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		del=d*c;
		h *= del;
		if(fabs(del-1.0) < EPS) break;
	}
	if(i > ITMAX) {
		write("erreur","\n a too large, ITMAX too small in routine gcf");
		exit(1);
	}
	*gammcf=exp(-x+a*log(x)-(*gln))*h;
}
		
#undef EPS
#undef ITMAX
#undef FPMIN

/*****************************************************************************/		

float betai(float a, float b, float x)
/*incomplete beta function or confluent hypergeometric function???p.227 Num recipe*/
{
	float bt;
	char smess[SMAX];

	if (x < 0.0 || x > 1.0) {
		sprintf(smess,"Error in betai: bad value of x = %7.5f ",x);
		write(ERRORFILE,smess);
	}
 	if (x == 0.0 || x == 1.0) bt = 0.0;
	else	bt=exp(gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.0-x));
	if (x < (a+1.0)/(a+b+2.0))
		return bt*betacf(a,b,x)/a;
	else
		return 1.0-bt*betacf(b,a,1.0-x)/b;

}

#define MAXIT 100
#define EPS 3.0e-7
#define	FPMIN 1.0e-30

float betacf(float a,float b,float x)
/*used by betai*/
{
	int m,m2;
	float aa,c,d,del,h,qab,qam,qap;

	qab=a+b;
	qap=a+1.0;
	qam=a-1.0;
	c=1.0;
	d=1.0-qab*x/qap;
	if (fabs(d) < FPMIN) d=FPMIN;
	d=1.0/d;
	h=d;
	for(m=1;m<=MAXIT;m++) {
		m2=2*m;
		aa=m*(b-m)*x/((qam+m2)*(a+m2));
		d=1.0+aa*d;
		if(fabs(d) < FPMIN) d=FPMIN;
		c=1.0+aa/c;
		if(fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		h *= d*c;
		aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
		d=1.0+aa*d;
		if(fabs(d) < FPMIN) d=FPMIN;
		c=1.0+aa/c;
		if(fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		del=d*c;
		h *= del;
		if(fabs(del-1.0) < EPS) break;
	}
	if(m>MAXIT) {
		write(ERRORFILE,"Error in betacf: a or b too big, or MAXIT too small");
	}
	return h;
}

#undef EPS
#undef MAXIT
#undef FPMIN

/***********************************************************************************/

double factrl(long n)
/*returns the value of n!*/
{
	static int ntop=4;
	static float a[33]={1.0,1.0,2.0,6.0,24.0};
	int j;

	if(n<0) {
		write(ERRORFILE,"error in factrl: negative value of n");
		exit(1);
	}
	if(n > 32) return exp(gammln(n+1.0));
	while (ntop<n) {
		j=ntop++;
		a[ntop]=a[j]*ntop;
	}
	return (double) a[n];
}

/***********************************************************************************/

