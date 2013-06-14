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
        int i, j;
        printf( "\n %s\n", desc );
   for( i = 0; i < n; i++ ) {
      j = 0;
      printf( "%6.2f" ,a[i]);
      while( j < n ) {
         if( wi[j] == (double)0.0 ) {
            printf( " %6.2f", v[i+j*ldv] );
            j++;
         } else {
            printf( " (%6.2f,%6.2f)", v[i+j*ldv], v[i+(j+1)*ldv] );
            printf( " (%6.2f,%6.2f)", v[i+j*ldv], -v[i+(j+1)*ldv] );
            j += 2;
         }
      }
      printf( "\n" );
   }
}




main(){

	const int N = 8 ;
	int i,j,k,l,m,rounds,n,s,ii,jj,ctr[400],count,nmax;
	int t1,t2,t3,t4,t5,t6;

	extern int seed;

  	double *dmn = malloc(N*N * sizeof(*dmn));
  	double *eigenVectors = malloc(N*N * sizeof(*eigenVectors));
  	double *eigenValuesr = malloc(N*N * sizeof(*eigenValuesr));
  	double *eigenValuesi = malloc(N*N * sizeof(*eigenValuesi));
//	double ranf();
	double p,q,ppr,qpr,pdp,qdp,qtr,ptr,eps,pi,r,ang,alp,kap,a[400],dx,x,y,lam,b;

	struct tm t;
	time_t ts, tp;

/*	FILE *f1;
	f1 = fopen("PF.dat","w");*/

	rounds=300000;
	pi = acos(-1.);
	eps=0.1;
	kap=0.001;
	alp=0.6;
	lam = 6.;
	b = 0.6;
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
	k=4;
	nmax=3;
	for(n=1;n<nmax;n++){
		count=0;
		a[0]=0.;
		dx=1./k;
		for(i=1;i<k+1;i++){
			a[i]=a[i-1]+dx;
			printf("a=%g\n",a[i]);}
		for(i=0;i<k;i++){
			for(j=0;j<k;j++)
				ctr[j]=0;
			for(m=0;m<rounds;m++){
				q=rand()/((RAND_MAX + 1.0))*dx;
//				printf("q=%g\n",q);
				x =a[i]+q;
				x = lam*x*(1-b*x)*(1-x);
//				l=0;
				for(j=0;j<k;j++){
						if(x>a[j] && x<a[j+1])
							ctr[j]++;}
//						l++;} // j
					} // m 
			for(j=0;j<k;j++){
				dmn[count] = ctr[j]/(double)rounds;
//				printf("%g\n",dmn[count]);
				count++;
				} // j
			} // i 

		MatrixComplexEigensystem( eigenVectors, eigenValuesr, eigenValuesi, dmn, k);

	        printf("\nEigenvectors\n");

 /* 		for(i=0;i<k;i++){
    			for(j=0;j<k;j++) printf(" (%g,%g) \t", creal(eigenVectors[i*k + j]), cimag(eigenVectors[i*k + j]));

    			printf("\n");
  			} // i  */

		 print_eigenvectors( "Left eigenvectors", k, eigenValuesi, eigenVectors,a, k );

  		printf("\nEigenvalues \n");
  		for(i=0;i<k;i++) printf("\n (%g +%g i) \t",  eigenValuesr[i], eigenValuesi[i]);


		printf("\n------------------------------------------------------------\n"); 


		k*=2;
		} // n 

	
	


	return 0;



//	fclose(f1);
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
