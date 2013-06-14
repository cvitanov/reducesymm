#include <stdio.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <complex.h>


int seed;

void zgeTranspose( double complex *Transposed, double complex *M ,int n)
{

int i,j;
for(i=0;i<n;i++)

  for(j=0;j<n;j++) Transposed[i+n*j] = M[i*n+j];
}

void dgeTranspose( double *Transposed, double *M ,int n)
{

int i,j;
for(i=0;i<n;i++)

  for(j=0;j<n;j++) Transposed[i+n*j] = M[i*n+j];
}

//......................................................................................................
//  MatrixComplexEigensystem: computes the eigenvectors and eigenValues of input matrix A
//  The eigenvectors are stored in columns
//.....................................................................................................
void MatrixComplexEigensystem( double *eigenvectorsVL, double *eigenvaluesWr, double *eigenvaluesWi, double *A, int N)
{

 int i;
  
  double  *AT = (double *) malloc( N*N*sizeof(double) );

  double *VL = (double*) malloc(N*N*sizeof(double));

  dgeTranspose( AT, A , N);
  
  char JOBVL ='V';   // Compute Left eigenvectors

  char JOBVR ='N';   // Do not compute Right eigenvectors

  double complex VT[1];

  int LDVL = N; 
  int LDVR = 1;

 int LWORK = 4*N; 

 double complex *WORK =  (double complex*)malloc( LWORK*sizeof(double complex));

 //double complex *RWORK = (double complex*)malloc( 2*N*sizeof(double complex));

int INFO;

 dgeev_( &JOBVL, &JOBVR, &N, AT ,  &N , eigenvaluesWr , eigenvaluesWi, 

   eigenvectorsVL, &LDVL,
   VT, &LDVR, 
   WORK, 
   &LWORK, &INFO );

//  zgeTranspose( VT, eigenvectorsVR , N);

// for(i=0;i<N*N;i++) eigenvectorsVR[i]=VT[i];

 
  free(WORK);
//  free(RWORK);
  free(AT);
}

/* Auxiliary routine: printing eigenvectors */
void print_eigenvectors( char* desc, int n, double* wi, double* v, double* a, int ldv ) {
        int i, j,l,m,bound;
	double dummy;
	FILE *f1;
	f1 = fopen("lozi_dens.dat","w");
	dummy = sqrt((double)n);
	bound = (int)dummy;
	l = 0;
	m = 0;
        printf( "\n %s\n", desc );
   for( i = 0; i < n; i++ ) {
      	j = 0;
	m++;
	if(m==bound){
		m=0;
		l++;
		fprintf(f1,"\n");} // if 
        fprintf(f1, "%6.2f %6.2f" ,a[l],a[m]);
      while( j < n ) {
         if( wi[j] == (double)0.0 ) {
             fprintf(f1, " %6.2f", fabs(v[i+j*ldv]) );
            j++;
         } else {
             fprintf(f1, " (%6.2f,%6.2f)", v[i+j*ldv], v[i+(j+1)*ldv] );
             fprintf(f1, " (%6.2f,%6.2f)", v[i+j*ldv], -v[i+(j+1)*ldv] );
            j += 2;
         }
      }
       fprintf(f1, "\n" );
   }

fclose(f1);

}




main(){

	const int N = 4096 ;
	int i,j,k,l,m,rounds,n,s,ii,jj,count,nmax,ctr[200][200];
	int t1,t2,t3,t4,t5,t6,dummy;
//	int  *ctr = malloc(N*N * sizeof(*ctr));

	extern int seed;

  	double *dmn = malloc(N*N * sizeof(*dmn));
  	double *eigenVectors = malloc(N*N * sizeof(*eigenVectors));
  	double *eigenValuesr = malloc(N*N * sizeof(*eigenValuesr));
  	double *eigenValuesi = malloc(N*N * sizeof(*eigenValuesi));
//	double ranf();
	double p,q,ppr,qpr,pdp,qdp,qtr,ptr,eps,pi,r,ang,alp,kap,a[400],dx,x,y,lam,b;

	struct tm t;
	time_t ts, tp;
	
	FILE *f2;
	f2 = fopen("lozi_evs.dat","w");

	rounds=90000;
	pi = acos(-1.);
	eps=0.1;
	kap=0.001;
	alp=1.85;
	lam = 6.;
	b = .3;
	ts = time(&tp);
  	t  = *gmtime(&tp);
  	t1 = t.tm_sec;
  	t2 = t.tm_min;
  	t3 = t.tm_hour;
  	t4 = t.tm_mday;
  	t5 = t.tm_mon;
  	t6 = t.tm_year;
  	seed = t6+70*(t5+12*(t4+31*(t3+23*(t2+59*t1))));
  	if (seed%2 == 0) seed = seed-1;
	printf("S = %d\n", seed);
	srand(time(NULL));
	dummy = (double)N;
	k=(int)sqrt(dummy);
	nmax=2;
	for(n=1;n<nmax;n++){
		count=0;
		dx=2.6/k;
		a[0]= -1.3 + dx/2;
		for(i=1;i<k+1;i++) // {
			a[i]=a[i-1]+dx;
//			printf("a=%g\n",a[i]);}
		for(i=0;i<k;i++){
			for(j=0;j<k;j++){
				for(l=0;l<k;l++){
					for(m=0;m<k;m++)
						ctr[l][m]=0; // m
					} // l
				for(m=0;m<rounds;m++){
					x=rand()/((RAND_MAX + 1.0))*dx;
					y=rand()/((RAND_MAX + 1.0))*dx;
//				printf("q=%g\n",q);
					q =a[i]+x;
					p =a[j]+y;
                                	qpr = 1- alp*fabs(q) + b*p;
                                	ppr = q;
					q = qpr;
                                        p = ppr;
					for(ii=0;ii<k;ii++){
						for(jj=0;jj<k;jj++){
							if(p>a[jj] && p<a[jj+1] && q>a[ii] && q<a[ii+1])
									ctr[ii][jj]++;
							} // jj
						} // ii
					} // m
			printf("i=%d j=%d ctr=%d\n",i,j,ctr[i][j]);
			for(ii=0;ii<k;ii++){
				for(jj=0;jj<k;jj++){
					dmn[count] = ctr[ii][jj]/(double)rounds;
//				        printf("dmn=%g count=%d\n",dmn[count],count);
					count++;
					} // jj
				} // ii
			} // j
		} // i 

		printf("i=%d j=%d ctr=%d\n",i,j,ctr[i][j]);


		MatrixComplexEigensystem( eigenVectors, eigenValuesr, eigenValuesi, dmn, N);

	        printf("\nEigenvectors\n");

 /* 		for(i=0;i<k;i++){
    			for(j=0;j<k;j++) printf(" (%g,%g) \t", creal(eigenVectors[i*k + j]), cimag(eigenVectors[i*k + j]));

    			printf("\n");
  			} // i  */

		 print_eigenvectors( "Left eigenvectors", N, eigenValuesi, eigenVectors,a, N );

  		fprintf(f2,"\nEigenvalues \n");
  		for(i=0;i<N;i++) fprintf(f2,"\n (%d %g +%g i) \t",i+1, eigenValuesr[i], eigenValuesi[i]);


		printf("\n------------------------------------------------------------\n"); 


		k*=2;
		} // n 

	
	


	return 0;



	fclose(f2);
	} // main







/* double ranf() {
// Function to generate a uniform random number in [0,1]
// following x(i+1)=a*x(i) mod c with a=pow(7,5) and
// c=pow(2,31)-1.  Here the seed is a global variable.

  const int a = 16807,  c = 2147483647, q = 127773, r = 2836;
  int l, h, t;
  double cd = c;
  extern int seed;

  h = seed/q;
  l = seed%q;
  t = a*l - r*h;
  if (t > 0) seed = t;
  else seed = c + t;
  return seed/cd;
} 
*/
