#include "stdio.h"
#include "stdlib.h"
#include "math.h"

#define FOURIER 1
/* whether using fft to get the nolinear terms */

#define ORBIT_ONLY 0
/* conditional compiling */

#define periodTime 10.25
/* the period of the orbit */

#define oldmethod 0
/* if using the old method */ 

#define cutoffNum 32
/* define cuttoff number for the Fourier modes */

#define stepSize 0.25
/* define the step size of the simulation */

#define maxTime 30000
/* define the max time step of the simulation */

//#define L 3.501408808
#define L 22.0
/* period of x in KS equation */

#define FourierN 1024
/* the total segmentation of the integral */

class complx
{
	public:
		double real, img;

	complx(double _real, double _img)
	{
		real = _real;
		img = _img;
	}

	complx()
	{
		real = 0.0;
		img = 0.0;
	}

	complx(const complx &c)
	{
		real = c.real;
		img = c.img;
	}

	complx(double phase)
	{
		real = cos(phase);
		img = sin(phase);
	}

	complx operator +(const complx &c)
	{
		return complx(real + c.real, img + c.img);
	}

	complx operator +(double d)
	{
		return complx(real + d, img);
	}

	complx operator -(const complx &c)
	{
		return complx(real - c.real, img - c.img);
	}

	complx operator -(double d)
	{
		return complx(real - d, img);
	}

	complx operator *(const complx &c)
	{
		return complx(real * c.real - img * c.img, real * c.img + img * c.real);
	}

	complx operator *(double d)
	{
		return complx(real * d, img * d);
	}

	complx operator /(double l)
	{
		return complx(real / l, img / l);
	}

	complx couple()
	{
		return complx(real, -img);
	}
};

class cVector
{
	public:
		complx terms[cutoffNum];

		cVector()
		{
			for(int i = 0; i < cutoffNum; i++)
				terms[i] = complx(0, 0);
		}
		
		cVector(complx * series)
		{
			for(int i = 0; i < cutoffNum; i++)
				terms[i] = series[i];
		}

		cVector(cVector* v1, const cVector& v2, int sign)
		{
			switch (sign){
				case 1:
					for(int i = 0; i < cutoffNum; i++)
						terms[i] = v1->terms[i] + v2.terms[i];
					break;
				case 2:
					for(int i = 0; i < cutoffNum; i++)
						terms[i] = v1->terms[i] - v2.terms[i];
					break;
				case 3:
					for(int i = 0; i < cutoffNum; i++)
						terms[i] = v1->terms[i] * v2.terms[i];
					break;
				case 4:
					for(int i = 0; i < cutoffNum; i++)
						terms[i] = complx(exp(v2.terms[i].real), 0);
					break;
			}
		}

		cVector(cVector* v, double d, int sign)
		{
			switch(sign){
				case 1:
					for(int i = 0; i < cutoffNum; i++)
						terms[i] = v->terms[i] + complx(d, 0);
					break;
				case 2:
					for(int i = 0; i < cutoffNum; i++)
						terms[i] = v->terms[i] * complx(d, 0);
					break;
			}
		}

		cVector(const cVector& vec)
		{
			for(int i = 0; i < cutoffNum; i++)
				this->terms[i] = vec.terms[i];
		}
			
		void base(int i)
		{
			terms[i] = complx(1, 0);
		}
		
		#if(oldmethod == 1)
		void plus(cVector &v)
		{
			for(int i = 0; i < cutoffNum; i++)
				terms[i] = terms[i] + v.terms[i];
		}

		void plusp(cVector *v)
		{
			for(int i = 0; i < cutoffNum; i++)
				terms[i] = terms[i] + v->terms[i];
			free(v);
		}

		void minus(cVector *v)
		{
			for(int i = 0; i < cutoffNum; i++)
				terms[i] = terms[i] - v->terms[i];
			free(v);
		}
		#endif
		
		cVector operator +(const cVector &v)
		{
			return cVector(this, v, 1);
		}
		
		cVector operator +(double p)
		{
			return cVector(this, p, 1);
		}
		
		cVector operator -(const cVector &v)
		{
			return cVector(this, v, 2);
		}
		
		cVector operator *(const cVector &v)
		{
			return cVector(this, v, 3);
		}
		
		cVector operator *(double m)
		{
			return cVector(this, m, 2);
		}

		complx operator ^(const cVector &v)
		{
			complx result;
			for(int i = 0; i < cutoffNum; i++)
			{
				complx couple = v.terms[i];
				result = result + (terms[i] * couple.couple());
			}
			double real = result.real;
			double img = result.img;
			return complx(real, img);
		}
	
		#if(oldmethod == 1)
		cVector * multiV(const complx &c)
		{
			cVector * result = new cVector();
			for (int i = 0; i < cutoffNum; i++)
				result->terms[i] = (this->terms[i]) * c;
			return result;
		}
		#endif

		double abs()
		{
			double absolute = 0;
			for(int i = 0; i < cutoffNum; i++)
				absolute += (this->terms[i] * (this->terms[i].couple())
				).real;
			return sqrt(absolute);
		}

		void divide(double l)
		{
			for(int i = 0; i < cutoffNum; i++)
				terms[i] = terms[i] / l;
		}

		void copy(const cVector &v)
		{
			for(int i = 0; i < cutoffNum; i++)
				terms[i] = v.terms[i];
		}

		void copygroup(cVector * origin, cVector * paste)
		{
			for(int i = 0; i < cutoffNum; i++)
				paste[i].copy(origin[i]);
		}

		void bezero()
		{
			for(int i = 0; i < cutoffNum; i++)
				terms[i] = complx(0, 0);
		}
		
		cVector expo(const cVector& v)
		{
			return cVector(NULL, v, 4);
		}
};


#if (oldmethod == 1)
currently not useful in this ETDRK4 scheme.
complx bchrTab[6][6] = {{complx(0.25, 0)},
						{complx(0.09375, 0), complx(0.28125, 0)},
						{complx(0.879381, 0), complx(-3.277196, 0), complx(3.320892, 0)},
						{complx(2.0324074, 0), complx(-8, 0), complx(7.1734893, 0), complx(-0.208090, 0)},
						{complx(-0.296296, 0), complx(2, 0), complx(-1.381676, 0), complx(0.4529727, 0), complx(-0.275, 0)},
						{complx(0.1185185, 0), complx(0, 0), complx(0.518986, 0), complx(0.50613149, 0), complx(-0.18, 0), complx(0.0363636, 0)}};
#endif

float fwdTime;
/* set the forward time of the simulation */

double h = stepSize;
/* stepsize of the simulation for short. */

cVector k[3];
/* k coefficients of the Runge-Kutta algorithms */

cVector vec[9];
/* the vectors used to complish the computation */

cVector a;
/* fourier modes */

complx stab[cutoffNum][cutoffNum];
/* stability matrix */

cVector crntGSVs[cutoffNum], tempGSVs[cutoffNum];  
/* define the current GSVs */

complx jacbterm;
/* currently using jacobian term */

complx nolinear[cutoffNum];
/* the nolinear part of the calculation */

double x[FourierN];
/* the value after reverse fourier transform, the value is u^2 */

FILE * readdata; 
FILE * resultData;
/* to read and save the data */

double sumLypExpo[cutoffNum];
/* cumulated logrithsm of the diagonal term of matrix R in QR decomposition */


void init()
{
	fwdTime = 0;
	readdata = fopen("orbit.txt", "r");
	resultData = fopen("result.txt", "w+");
	if(resultData == NULL || readdata == NULL)
	{
		printf("error: cannot open the file\n");
		exit(0);
	}
	int k = 0;
	float tempr, tempi;
	while(k < cutoffNum)
	{
		fscanf(readdata, "%f\n%f", &tempr, &tempi);
		a.terms[k].real = tempr;
		a.terms[k].img = tempi;
		k++;
	}
	fclose(readdata);

	for(int i = 0; i < cutoffNum; i++)
		crntGSVs[i].base(i);
	/* initiate for the L, Lh, Exp[Lh/2], Exp[Lh] ... etc. */
	
	for(int i = 0; i < cutoffNum; i++)
	{
		vec[0].terms[i] = complx((pow((double)(i / L), 2) - pow((double)(i / L) , 4)), 0);
		vec[1].terms[i] = complx(1.0 / (pow((double)(i / L), 2) - pow((double)(i / L) , 4)), 0);
	}
	/* initiate L and L^(-1)*/
	vec[2] = vec[0].expo(vec[0] * (h / 2.0));
	vec[3] = vec[1] * (vec[2] + (-1.0));
	vec[4] = vec[0].expo(vec[0] * h);
	vec[5] = (vec[1] * vec[1] * vec[1]) * (1.0 / (h * h));
	vec[6] = vec[5] * (vec[0] * (-h) + vec[4] * (vec[0] * vec[0] * h * h + vec[0] * (-h) * 3.0 + 4.0) + (-4.0));
	vec[7] = vec[5] * (vec[4] * (vec[0] * h + (-2.0)) + vec[0] * h + 2.0) * 2.0;
	vec[8] = vec[5] * (vec[4] * (vec[0] * (-h) + 4.0) - vec[0] * vec[0] * h * h - vec[0] * h * 3.0 + (-4.0));
	/* modification */
	vec[3].terms[0] = complx(0.1250, 0);	vec[3].terms[22] = complx(0.1250, 0);	
	vec[6].terms[0] = complx(0.0417, 0);	vec[6].terms[22] = complx(0.0417, 0);
	vec[7].terms[0] = complx(0.0833, 0);	vec[7].terms[22] = complx(0.0833, 0);
	vec[8].terms[0] = complx(0.0417, 0);	vec[8].terms[22] = complx(0.0417, 0);
	/* modification */
}

#if (FOURIER == 1)
complx fft(int index, double kz)
{
	complx xF(0, 0);
	for(int i = 0; i < FourierN; i++)
		xF = xF + complx(-kz * i * index) * x[i];
	return xF;
}
#endif

cVector nLinear(cVector& v)
{
	#if (FOURIER == 0)
	complx tempa, tempb;
	for(int i = 0; i < cutoffNum; i++)
	{	
		nolinear[i] = complx(0, 0);
		for(int j = cutoffNum - 1; (i - j) < cutoffNum; j--)
		{
			int cp = i - j;
			if(j >= 0 && cp >= 0){
				tempa = v.terms[j];
				tempb = v.terms[cp];
			}
			else if(j < 0 && cp >= 0){
				tempa = v.terms[-j].couple();
				tempb = v.terms[cp];
			}
			else if(j >= 0 && cp < 0){
				tempa = v.terms[j];
				tempb = v.terms[-cp].couple();
			}
			else if(j < 0 && cp < 0){
				tempa = v.terms[-j].couple();
				tempb = v.terms[-cp].couple();
			}
			nolinear[i] = nolinear[i] + tempa * tempb;
		}
		nolinear[i] = nolinear[i] * 0.5;
		nolinear[i] = nolinear[i] * complx(0, -(double) i / L);
	}
	return cVector(nolinear);
	#endif
	
	#if (FOURIER == 1)
	double kz = 2.0 * 3.1415926 / FourierN;
	for(int i = 0; i < FourierN; i++)
	{
		x[i] = 0.0;
		complx temp(0, 0);
		for(int j = 0; j < cutoffNum; j++)
			temp = temp + v.terms[j] * complx(kz * i * j) + v.terms[j].couple() * complx(kz * i * (FourierN - j));
		x[i] =  pow((temp.real / FourierN), 2);
	}
	/* First simulate the series of x. */

	for(int i = 0; i < cutoffNum; i++)
	{
		nolinear[i] = fft(i, kz) * complx(0, -(double) i / L);
		nolinear[i] = nolinear[i] * 0.5;
	}
	/* Then, calculate the fourier series. */
	return cVector(nolinear);
	#endif
}

void sim()
{
	while(fwdTime <= periodTime)
	{
		k[0] = vec[2] * a + vec[3] * nLinear(a);
		k[1] = vec[2] * a + vec[3] * nLinear(k[0]);
		k[2] = vec[2] * k[0] + vec[3] *(nLinear(k[1]) * 2.0 - nLinear(a));
		a = vec[4] * a + vec[6] * nLinear(a) + vec[7] * (nLinear(k[0]) + nLinear(k[1])) + vec[8] * nLinear(k[2]);
		fwdTime += h;
		printf("%f\n", fwdTime);
	}
	while(1);
}

#if (oldmethod == 1)
complx fODE(int j, int i)
	/* the ode function to calculate k coefficient */
{
	for(int s = 0; s < cutoffNum; s++)
		tempa[s] = a[s];
	/* use a temp */
	for(int m = 0; m < j; m++)
		for(int n = 0; n < cutoffNum; n++)
			tempa[n] = tempa[n] + bchrTab[j][m] * k[n][m];
	/* addup using former result */
	complx sigma(0, 0);
	for(int x = 0; x < cutoffNum; x++)
	{
		if((abs(i - x)) <= cutoffNum)
			sigma = sigma + tempa[x] * tempa[i - x + cutoffNum];
	}
	i -= cutoffNum;
	complx part2 = ((complx(0, 1) * (i / L * 0.5)) * sigma);
	complx part1 = (tempa[i + cutoffNum] * (pow((double)(i / L), 2) - pow((double)(i / L) , 4)));
	return part1 - part2;
}
#endif

#if (ORBIT_ONLY == 1)
void lyaExpo()
	/* to get the lyapunov exponent of the system */
{
	fprintf(resultData, "Time is %d\n", fwdTime);
	for(int i = 0; i < cutoffNum; i++)
		fprintf(resultData, "The %dth LE: %f\n", i, sumLypExpo[i] / (stepSize * fwdTime));
}
#endif

#if (ORBIT_ONLY == 1)
complx jacobian(int i, int j)
{
	if(i == j)
	{
		jacbterm = complx(exp((pow((double)(i / L), 2) - pow((double)(i / L) , 4))) * stepSize, 0) - complx(0, 1) * (a[cutoffNum] * (double)(i / L) * stepSize);
		return jacbterm;
	}
	else
		if(abs(i - j) > cutoffNum)
			return complx(1, 0);
		else
			return complx(0, 1) * (a[i - j + cutoffNum] * (double)(i / L) * stepSize);
}
#endif

#if (ORBIT_ONLY == 1)
void fwdSim()
{
	while(fwdTime++ < maxTime)
	{
		for(int j = 0; j < 5; j++) /* calculate the k coefficients */
			for(int i = 0; i < cutoffNum; i++)
				k[i][j] = fODE(j, i) * stepSize;

		for(int i = 0; i < cutoffNum; i++) /* calculate using k coefficients */
			for(int c = 0; c < 5; c++)
				a[i] = a[i] + k[i][c] * bchrTab[5][c];
		 
		#if(ORBIT_ONLY == 1) /* if 0, only compute the orbit part. */
		
			crntGSVs[0].copygroup(crntGSVs, tempGSVs); /* copy the current GSVs to temp GSVs */
			for(int i = 0; i < cutoffNum; i++) /* calculate the e_i = L_ij * e_j */
			{
				crntGSVs[i].bezero();
				for(int j = 0; j < cutoffNum; j++)
					crntGSVs[i].plusp(tempGSVs[j].multiV(jacobian(i - cutoffNum, j - cutoffNum))); 
			}	

			for(int n = 0; n < cutoffNum; n++) /* GS process to calculate the lyapunov exponents */
			{
				cVector temp;
				temp.copy(crntGSVs[n]);
				for(int m = 0; m < n; m++)
					crntGSVs[n].minus(crntGSVs[m].multiV(temp * crntGSVs[m]));
				crntGSVs[n].divide(crntGSVs[n].abs());

				double diagonal = pow((temp * crntGSVs[n]).real, 2) + pow((temp * crntGSVs[n]).img, 2);
				sumLypExpo[n] = sumLypExpo[n] + 0.5 * log(diagonal); /* sum over the lyapunove exponents */
			}

			if((fwdTime % 100) == 0) /* generate the result */
			{
				printf("check the time: %d\n", fwdTime);
				lyaExpo();
			}
			
		#else
			
			if(fwdTime == Period)
				for(int i = 0; i < cutoffNum; i++)
					fprintf(resultData, "a[%d] = %f\n", i, a[i]);
			
		#endif
	}
	fclose(resultData);
}
#endif
int main()
{
	init();
	sim();
}
