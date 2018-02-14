

#include "gamma.h"

static void gser(double *gamser, const double a, const double x, double *gln);
static void gcf(double *gammcf, const double a, const double x, double *gln);


static double betacf(const double a, const double b, const double x);


double gammln(const double xx)
{	int j;
	double x,y,tmp,ser;
	static const double cof[6] = {76.18009172947146,-86.50532032941677,24.01409824083091,-1.231739572450155,0.1208650973866179e-2,-0.5395239384953e-5};

	y = x = xx;
	tmp = x + 5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for(j=0;j<6;j++) ser += cof[j]/++y;
	return -tmp + log(2.5066282746310005*ser/x);
}


static void gser(double *gamser, const double a, const double x, double *gln)
{	const int ITMAX = 1000;
	double EPS = 1e-100;
	int n;
	double sum, del, ap;

	EPS = (EPS>DBL_MIN) ? EPS : DBL_MIN;

	//*gln = gammln(a);
	*gln = lgamma(a);
	if(x<=0.0){
		if(x<0.0) {fprintf(stderr,"x less than 0 in routine gser\n");exit(1);}
		*gamser = 0.0;
		return;
		}else{
		ap = a;
		del = sum = 1.0/a;
		for(n=0;n<ITMAX;n++){
			++ap;
			del *=x/ap;
			sum += del;
			if(fabs(del)<fabs(sum)*EPS){
				*gamser = sum*exp(-x+a*log(x)-(*gln));
				return;
				}
			}
		fprintf(stderr,"a too large, ITMAX too small in routine gser\n");
		}

	return;
}


static void gcf(double *gammcf, const double a, const double x, double *gln)
{	const int ITMAX = 1000;
	const double EPS = 1e-100;
	const double FPMIN = DBL_MIN/EPS;
	int i;
	double an, b, c, d, del, h;

	//*gln = gammln(a);
	*gln = lgamma(a);
	b = x+1.0-a;
	c = 1.0/FPMIN;
	d = 1.0/b;
	h = d;
	for(i=1;i<=ITMAX;i++){
		an = -i*(i-a);
		b += 2.0;
		d=an*d+b;
		if(fabs(d)<FPMIN) d=FPMIN;
		c=b+an/c;
		if(fabs(c)<FPMIN) c=FPMIN;
		d=1.0/d;
		del=d*c;
		h *= del;
		if(fabs(del-1.0)<=EPS) break;
		}
	if(i>ITMAX) {fprintf(stderr,"a too large, ITMAX too small in gcf\n");}

	*gammcf = exp(-x+a*log(x)-(*gln))*h;
}

double gammp(const double a, const double x)
{	double gamser, gammcf, gln;
	if(x<0.0||a<=0){
		fprintf(stderr,"Invalide arguments in routine gammp\n");
		exit(1);
		}
	if(x < a+1.0){ 
		gser(&gamser,a,x,&gln);
		return gamser;
		} else{
		gcf(&gammcf,a,x,&gln);
		return 1-gammcf;
		}
}


double gammq(const double a, const double x) // upper incomplete function Q(a,x) = 1 - P(a,x)
{	double gamser, gammcf, gln;
        if(x<0.0||a<=0){
                fprintf(stderr,"Invalide arguments in routine gammq\n");
                exit(1);
                }
	if(x < a+1.0){
		gser(&gamser,a,x,&gln);
		return 1-gamser;
		}else{
		gcf(&gammcf,a,x,&gln);
		return gammcf;
		}
}


double pchisq(int k, const double x, int lower_tail)
{	double p,a,q;
	a = ((double) k)/2;
 	if(k<=0){
		fprintf(stderr,"degree of freedom (%d) of chi-sq distribution must be postive\n",k);
	exit(1);
	}
	q = (x<0.0 ? 0.0 : x/2);

	if(lower_tail==1){
		p = gammp(a,q);
		}else{
		p = gammq(a,q);
		}
	return p;
}


double qchisq_newton(int k, const double p_in, int lower_tail)
{	double q, q1,p=p_in,fprime=1, sf,a, sign;
	int MAX_ITR = 10000,i;
	const double EPS = 1e-10;

	if(k<=0){
		fprintf(stderr,"degree of freedom (%d) of chi-sq distribution must be postive\n",k);
		exit(1);
		}
	if(p<0||p>1){
		fprintf(stderr,"probability must be between 0 and 1\n");
		exit(1);
		}

	if(lower_tail==1){ sign = 1;} else {sign=-1;}

	a = q = (double) k;
	q1 = q+1;
	i = 0;
	//sf = -(a/2.0)*log(2.0)-gammln(a/2.0);
	sf = -(a/2.0)*log(2.0)-lgamma(a/2.0);
	while(fabs(q-q1)/q1>EPS&&i<MAX_ITR){
		q1 = q;
		fprime = (a/2.0-1.0)*log(q1) - q1/2.0 +sf;
		fprime = sign*exp(fprime);
		q = q1 - (pchisq(k,q1,lower_tail)-p)/fprime;
		i++;
		}
	if(i>=MAX_ITR && fabs(q-q1)/q1>EPS) fprintf(stderr,"Newton-Raphson didn't converge\n");

	return q;
}


double qchisq_bisection(int k, const double p_in, int lower_tail)
{       double q1, q0,q_mid,p = p_in, mu, sd;
        int MAX_ITR = 10000,i;
        const double EPS = 1e-10;

        if(k<=0){
                fprintf(stderr,"degree of freedom (%d) of chi-sq distribution must be postive\n",k);
                exit(1);
                }
        if(p<0.0||p>1.0){
                fprintf(stderr,"probability must be between 0 and 1\n");
                exit(1);
                }

	mu = (double) k;
	sd = sqrt(2*mu);
	
        if(lower_tail==1){ 
		if(p==0.0) return 0.0; else if(p==1.0) return DBL_MAX; //INFINITY;
		q1 = mu + 3*sd; i = 0;
		while(pchisq(k,q1,0) > 1- p) {q1 += 2*sd;}
		
		q0 = 0.0;i=0;
		while(fabs(q0-q1)/q1>EPS&&i<MAX_ITR){
			q_mid = q0 + (q1-q0)/2.0;
			if(pchisq(k,q_mid,1)<p) q0 = q_mid;
			else if(pchisq(k,q_mid,1)>p) q1 = q_mid;
			else {q0 = q_mid; q1 = q_mid;  break;}
			i++;
			}
		}else{
		if(p==1.0) return 0.0; else if(p==1.0) return DBL_MAX;  //INFINITY;
		q1 = mu + 3*sd;
		while(pchisq(k,q1,0) > p) q1 += 2*sd;
		q0 = 0.0;i=0;
		while(fabs(q0-q1)/q1>EPS&&i<MAX_ITR){
			q_mid = q0 + (q1-q0)/2.0;
			if(pchisq(k,q_mid,0)<p) q1 = q_mid;
			else if(pchisq(k,q_mid,0)>p) q0 = q_mid;
			else {q0 = q_mid;q1 = q_mid; break;}
			i++;
			}
		}

	if(fabs(q0+q1)/2.0<DBL_MIN*10) {q0=0.0;q1=0.0;i=0;} /*too small, treat it as zero*/
	if(i>=MAX_ITR) fprintf(stderr,"Error in qchisq: MAX_ITR too small\n");
	fprintf(stderr,"fabs(q0-q1)/q1=%g\n",fabs(q0-q1)/q1);
	fprintf(stderr,"fabs(q0-q1)=%g\n",fabs(q0-q1));
        return (q0+q1)/2.0;
}


double qchisq(int k, const double p, int lower_tail)
{       double q; 
	if(k<=0){
                fprintf(stderr,"degree of freedom (%d) of chi-sq distribution must be postive\n",k);
                exit(1);
                }
        if(p<0.0||p>1.0){
                fprintf(stderr,"probability must be between 0 and 1\n");
                exit(1);
                }
	if(p<1e-20||(1-p)<1e-10)  q = qchisq_bisection(k,p,lower_tail); /*p or 1-p too small, Newton-Raphson may not work*/
	else q = qchisq_newton(k,p,lower_tail);

	return q;
}





static double betacf(const double a, const double b, const double x)
{	const int MAXIT = 1000;
	const double EPS = 1e-50;
	const double FPMIN = DBL_MIN/EPS;
	int m, m2;
	double aa, c, d, del, h, qab, qam, qap;

	qab = a+b;
	qap = a+1.0;
	qam = a-1.0;
	c = 1.0;
	d = 1.0 - qab*x/qap;
	if(fabs(d)<FPMIN) d = FPMIN;
	d = 1.0/d;
	h = d;
	for(m=1;m<=MAXIT;m++){
		m2=2*m;
		aa = m*(b-m)*x/((qam+m2)*(a+m2));
		d = 1.0 + aa*d;
		if(fabs(d)<FPMIN) d = FPMIN;
		c = 1.0 + aa/c;
		if(fabs(c)<FPMIN) c = FPMIN;
		d = 1.0/d;
		h *= d*c;
		aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
		d = 1.0 + aa*d;
		if(fabs(d)<FPMIN) d = FPMIN;
		c = 1.0 + aa/c;
		if(fabs(c)<FPMIN) c = FPMIN;
		d = 1.0/d;
		del = d*c;
		h *=del;
		if(fabs(del-1.0)<=EPS) break;
		}
	if(m > MAXIT) fprintf(stderr,"a or b too big, or MAXIT too small in betacf\n");


	return h;
}


double betai(const double a, const double b, const double x)
{	double bt, rslt;
	if(x < 0.0 || x> 1.0) {fprintf(stderr,"Bad x in routine betai\n");exit(1);}
	if(x==0.0||x==1.0) bt = 0.0;
	else{
		//bt = exp(gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.0-x));
		bt = exp(lgamma(a+b)-lgamma(a)-lgamma(b)+a*log(x)+b*log(1.0-x));
		}

	if(x < (a+1.0)/(a+b+2.0)) rslt = bt*betacf(a,b,x)/a;
	else rslt = 1.0 - bt*betacf(b,a,1.0-x)/b;

	return rslt;
}


double pbinom(const double k, const double n, const double p, const int lower_tail)
{	double prob;
	double k1;
	
	if(p< 0.0 || p> 1.0) {fprintf(stderr,"Bad probability parameter in routine pbinom\n");exit(1);}
	if(n<=0) {fprintf(stderr,"Bad size parameter in routine pbinom\n");exit(1);}
	if(lower_tail==1){
		k1 = n - k;
		if(k1<=0) prob = 1.0;
		else if(k1>n) prob = 0.0;
		else prob = betai(k1,n-k1 +1.0, 1-p);
		}else{
		if(k<=0) prob = 1.0;
		else if(k>n) prob = 0.0;
		else prob =  betai(k, n-k+1.0, p);
		}

	return prob;
}






static int idum=-849345; /*the seed of the random generator rand_lp*/

void seed_set(int seed)
        { idum = seed;}

int genSeed()
{	struct timeval tv;
	int seed;
	gettimeofday(&tv, NULL);
	seed = tv.tv_sec*1000000 + tv.tv_usec;

	return seed;
}

double rand_lp()
        { const int IM1 = 2147483563, IM2 = 2147483399;
          const int IA1 = 40014, IA2 = 40692, IQ1 = 53668, IQ2 = 52774;
          const int IR1 = 12211, IR2 = 3791, NTAB = 32, IMM1 = IM1-1;
          const int NDIV = 1+IMM1/NTAB;
          const double EPS = 3.0e-16, RNMX = 1.0- EPS, AM = 1.0/((double) IM1);
          static int idum2 = 123456789, iy=0;
          static int iv[32];
          int j,k;
          double temp;

          //iv = (int *) malloc(sizeof(int)*(NTAB)); // initialize iv as a vector

          if(idum<=0){ //Initialize
                  idum = (idum==0?1:-idum);
                  idum2 = idum;
                  for(j=NTAB+7;j>=0;j--){ // Load the shuffle talbe (after 8 warm ups)
                          k = idum/IQ1;
                          idum = IA1*(idum-k*IQ1)-k*IR1;
                          if(idum<0) idum += IM1;
                          if(j< NTAB) iv[j] = idum;
                        }
                  iy = iv[0];
                }
          k = idum/IQ1; // Start here when not initializing.
          idum = IA1*(idum-k*IQ1)-k*IR1;
          if(idum<0) idum +=IM1;
          k = idum2/IQ2;
          idum2=IA2*(idum2-k*IQ2)-k*IR2;
          if(idum2<0) idum2 +=IM2;
          j = iy/NDIV;
          iy = iv[j]-idum2;
          iv[j] = idum;
          if(iy<1) iy += IMM1;
          if((temp=AM*iy)>RNMX) return RNMX;
          else return temp;
        }



void db_shuffle(double *array, int n)
        { int i,k;
          double tmp;
          k = n;
          while(k>1)
                {
                  i = (int) floor(k*rand_lp());
                  k--;
                  tmp = array[k];
                  array[k] = array[i];
                  array[i] = tmp;
                }
          return;
        }

void db_shuffle_int(int *array, int n)
        { int i,k;
          int tmp;
          k = n;
          while(k>1)
                {
                  i = (int) floor(k*rand_lp());
                  k--;
                  tmp = array[k];
                  array[k] = array[i];
                  array[i] = tmp;
                }
          return;
        }

/*
double rgamma1(double alpha, double beta)
        { double U0, U1;
          double x, y;
          const double e = 2.718281828459045; //Euler's number

          if(alpha<=0||beta<=0) {printf("Error in rgamm1: invalid parameters.\n"); return -1;}

          if(alpha==1.0) //when alpha=1, gamma distribution is exponential distribution; simple generating method can be used
                { U0 = rand_lp();
                  x = -log(1-U0);
                  return beta*x;
                }


          do{
                  U0 = rand_lp();
                  U1 = rand_lp();
                  if(U0>e/(e+alpha))
                        { x = -log((alpha+e)*(1.0-U0)/(alpha*e));
                          y = pow(x,alpha-1.0);
                        }
                  else
                        { x = pow((alpha+e)*U0/e,1.0/alpha);
                          y = exp(-x);
                        }
                }while(U1>=y);

          return beta*x;
        }*/


double rgamma1(double alpha_in, double beta)
{	double c1, c2 ,c3, c4, c5;
	double U1,U2, W,U;
	double x=1.0;
	double alpha;
	int flag = 0;
	int i=0;

	if(alpha_in<=0||beta<=0) {fprintf(stderr,"Error in rgamm1: invalid parameters.\n"); return -1;}

        if(alpha_in==1.0) //when alpha=1, gamma distribution is exponential distribution; simple generating method can be used
                { U = rand_lp();
                  x = -log(1-U);
                  return beta*x;
                }


	if(alpha_in<1.0){
		alpha = alpha_in +1.0;
		}else{
		alpha = alpha_in;
		}

	c1 = alpha-1.0;
	c2 = (alpha-(1.0/(6.0*alpha)))/c1;
	c3 = 2.0/c1;
	c4 = 1.0+c3;
	c5 = 1.0/sqrt(alpha);

	flag = 0;
	do{
		do{	//fprintf(stderr,"i=%d\n",i); i++;
			U1 = rand_lp();
			U2 = rand_lp();
			if(alpha>2.5) U1 = U2 + c5*(1-1.86*U1);
			//fprintf(stderr,"U1=%g\nU2=%g\n",U1,U2);
			i++;
			} while(U1<=0.0||U1>=1.0);
		W = c2*U2/U1;
		if(c3*U1+W+1.0/W <=c4 || c3*log(U1)-log(W)+W<=1.0) {x=c1*W; flag = 1;}
		}while(flag==0);


	if(alpha_in<1.0){
		U = rand_lp();
		x = x*pow(U,1.0/alpha_in);
		}

	return beta*x;
}


double rbeta(double alpha, double beta)
{	double x,y;
	if(alpha<=0.0||beta<=0.0) {fprintf(stderr,"Error in rbeta: invalid parameters\n");return -1;}

	x = rgamma1(alpha,1);
	y = rgamma1(beta,1);

	return x/(x+y);
}



int rDirichlet(double *alpha, double *x,int n)
        { double sum;
          int i;

          if(alpha==NULL||x==NULL) {printf("Error in rDirichlet.\n"); return -1;}
          sum=0.0;
          for(i=0;i<n;i++)
                { x[i] = rgamma1(alpha[i],1.0);
                  if(x[i]<0) return -1;
                  sum = sum+x[i];
                }
          for(i=0;i<n;i++)
                { x[i] = x[i]/sum;}
          return 1;
        }
                    


double rpois(const double xm)
{	const double PI = 3.141592653589793238;
	static double sq, alxm, g, oldm = (-1.0);
	double em, t, y; 

	if(xm<=0.0) {fprintf(stderr,"Error in rpois: the parameter (=%g) must be postive\n",xm); exit(1);}

	if(xm<12.0){ /*use directy method*/
		if(xm!=oldm){
			oldm = xm;
			g = exp(-xm);
			}

		em = -1;
		t = 1.0;
		do{
			++em;
			t *= rand_lp();
			}while(t > g);
		}else{
		if(xm!=oldm){
			oldm = xm;
			sq = sqrt(2.0*xm);
			alxm = log(xm);
			g = xm*alxm - lgamma(xm+1.0);
			}
		do{
			do{
				y = tan(PI*rand_lp());
				em = sq*y+xm;
				}while(em<0.0);
			em = floor(em);
			t = 0.9*(1.0+y*y)*exp(em*alxm-lgamma(em+1.0)-g);
			}while(rand_lp()>t);
		}
	return em;
}



double rbinom(const double pp, const int n)
{	const double PI = 3.141592653589793238;
	int j;
	static int nold=(-1);
	double am, em, g, angle, p, bnl, sq, t, y;
	static double pold = (-1.0), pc, plog, pclog, en, oldg;

	p = (pp <=0.5? pp : 1.0-pp);
	am = n*p;
	if(n<25){ /*when n is small use direct method to generate random number from binomial distribution*/
		bnl = 0.0;
		for(j=0;j<n;j++){
			if(rand_lp() < p) ++bnl;
			}
		}else if(am<1.0){ /*If fewer than one event is expected out of 25 or more trials, then the distribution is quite accuratly Poisson, Use direct Poisson method*/
		g = exp(-am);
		t = 1.0;
		for(j=0;j<=n;j++){
			t *= rand_lp();
			if(t<g) break;
			}
		bnl = (j <=n ? j:n);
		}else{ /*use rejection method*/
		if(n != nold){
			en = n;
			//oldg = gammln(en+1.0);
			oldg = lgamma(en+1.0);
			nold = n;
			}
		if(p != pold){
			pc = 1.0 -p;
			plog = log(p);
			pclog = log(pc);
			pold = p;
			}

		sq = sqrt(2.0*am*pc);
		do{
			do{
				angle = PI*rand_lp();
				y = tan(angle);
				em = sq*y+am;
				}while(em < 0.0 || em >= (en+1.0)); /*Reject*/
				em = floor(em); /*Trick for integer-valued distribution*/
				t = 1.2*sq*(1.0+y*y)*exp(oldg-lgamma(em+1.0)-lgamma(en-em+1.0)+em*plog+(en-em)*pclog);
			}while(rand_lp()>t);
			bnl = em;
		}
		
		if(p != pp) bnl = n-bnl;
		return bnl;
}


double rnegbinom(const double n, const double p)
{	double alpha, beta;
	double y, em;

	if(p>=1.0||p<0.0) {fprintf(stderr,"Error in rnegbinom, p=%g is out of expected range\n",p); exit(1);}
	if(n<=0.0) {fprintf(stderr,"Error in rnegbinom, p=%g is out of expected range\n",p);exit(1);}
	if(p==0.0) return 0.0; 
	alpha = n;
	//beta = (1.0-p)/p;
	beta = p/(1.0-p);

	//fprintf(stderr,"alpha=%g\nbeta=%g\n",alpha,beta);
	y = rgamma1(alpha,beta);
	//fprintf(stderr,"y=%g\n",y);
	em = rpois(y);
	//if((int) em < 0) fprintf(stderr,"alpha=%g,beta=%g,y=%g,em=%g\n",alpha,beta,y,em);
	return em;
}

double rnegbinom_mv(const double m, const double v)
{	double n, p;
	double em;
	
	if(m<0.0||v<m) {fprintf(stderr,"Error in rnegbinom_mv: invalid parameters in rnegbinom_mv: mean = %g, variance = %g\n",m,v);} 

	n = m*m/(v-m);
	p = (v-m)/v;

	//fprintf(stderr,"n=%g\np=%g\n",n,p);
	em = rnegbinom(n,p);
	//if((int) em< 0) {fprintf(stderr,"mean=%g,variance = %g, n= %g, p = %g , em = %g\n",m,v,n,p,em);}
	return em;
}

double rnorm(const double mu, const double sd) /*see page 43 Robert and Casella, Monte Carlo Statistical Methods*/
{	double U1, U2;
	static double x1=0.0, x2=0.0;
	static int flag=0;
	const double PI = 3.141592653589793238;
	
	if(sd<0.0) {fprintf(stderr,"Error in rnorm: standard deviation of normal distribution must be positive\n");exit(1);}
	if(flag==0){
		U1 = rand_lp();
		U2 = rand_lp();
		x1 = sqrt(-2.0*log(U1))*cos(2.0*PI*U2);
		x2 = sqrt(-2.0*log(U1))*sin(2.0*PI*U2);
		flag = 1;
		return sd*x1+mu; /*Box-Muller returns two independent standard normal variable, return one, save the other for the next call*/
		}
	else {	flag = 0;
		return sd*x2 + mu;
		}
}
