float ran1(long  *idum);
float ran3(long *idum);
float poidev(float xm, long  *idum); 
float gammln(float xx);
long int timeseed(void);
void polint(double xa[],double ya[],int n, float x,float *y, float *dy); 
double trapzd(double (*func)(float), float a, float b, int n);
float qromb(double (*func)(float),float a, float b);
double gammaln(double xx);
double poisson(double xm, long int *idum);
int zbrac(double (*func)(float), float *x1, float *x2);
float zbrent(double (*fx) (float),float x1, float x2, float tol);
void zbrak(double (*fx) (float), float x1, float x2, int n, float xb1[],float xb2[], int *nb);
float bico(int n, int k);
float brent(float ax, float bx, float cx,float tol,float *xmin,double (*func)(float)); 
void mnbrak(float *ax, float *bx, float *cx, float *fa, float *fb, float *fc,
	double (*func)(float));
#define SIGN(a,b) ((b) > 0.0 ? fabs(a) : -fabs(a))
static float maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ? (maxarg1) : (maxarg2))
static float sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg) 
float erff(float x);
float gammp(float a, float x);
void gser(float *gamser,float a, float x, float *gln);
void gcf(float *gammcf,float a,float x,float *gln);
float betai(float a, float b, float x);
float betacf(float a,float b,float x);
double factrl(long n);




