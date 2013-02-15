/*
Compute Lyapunov exponents & vectors
*/

#define OMPSW 0		//for using openMP (0: no openMP, 1: for many vectors, 2: for many DOFs)
#define OMPRSW 0	//for reproducibility in openMP (0:reproducible (required for loop division), 1:irreproducible);
#define ARMASW 0	//use library Armadillo (required for calculating angle between subspaces)
#define FFTWSW 0	//use FFTW (required for 2d-KS);

#define VERSION -4
#define NUMDOSV 8
#define NUMSUPPLPTCP 4
#define PICONST 3.141592653589793

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cfloat>
#include <vector>
#include <new>
#include "mt19937ar-kaz.cpp"
#include "fft.h"	//NOTE: x should be read as -x for code using fft.h
#if FFTWSW
	#include <complex>
	#include <fftw3.h>
#endif
#if OMPSW
	#include "omp.h"
#endif
#if ARMASW
	#include "armadillo"
	using namespace arma;
#endif
using namespace std;


#define DEFAULTMEMORY 256	//(in MB);
#define MAXFILESIZE 2048	//(in MB);
//#define MAXFILESIZE 1792

#define MAX(a,b) ((a)>(b) ? (a) : (b))
#define MIN(a,b) ((a)<(b) ? (a) : (b))
#define SINGLE(a) (static_cast<double>(static_cast<float>(a)))
#define INTSIGN(a) ((a)>0 ? +1 : ((a)<0 ? -1 : 0))
#define DBLSIGN(a) ((a)>0.0 ? +1 : ((a)<0.0 ? -1 : 0))

//specified file exists or not;
bool isFileExisting(const char *FileName) {
	bool out;
	ifstream File;
	File.open(FileName);
	out = File;
	File.close();
	return out;
}

//namespace for OpenMP;
namespace omp {
	int NumThreads, ThreadID;		//number & ID of thread;
#if OMPSW
	#pragma omp threadprivate(ThreadID)
	void Init() {
		//Init OpenMP;
		omp_set_dynamic(0);
		#pragma omp parallel
		{
			#pragma omp single
			NumThreads = omp_get_num_threads();
			ThreadID = omp_get_thread_num();
		}
	}
#else
	void Init() {
		NumThreads = 1;
		ThreadID = 0;
	}
#endif
}

//End program with displaying error message;
void ErrorExit(string Mes) {
	cout << Mes << endl;
	exit(1);
}

// Code for random shuffling;

#define SWAP(a,b) itemp=(a);(a)=(b);(b)=itemp;
#define M 7
#define NSTACK 100

static void indexx(int n, int arr[], int indx[])
{
  int i,indxt,ir=n,itemp,j,k,l=1;
  int jstack=0,*istack,a;

  istack = new int[NSTACK+1];
  for(j=1;j<=n;j++) indx[j]=j;
  for(;;){
    if(ir-l<M){
      for(j=l+1;j<=ir;j++){
        indxt=indx[j];
        a=arr[indxt];
        for(i=j-1;i>=1;i--){
          if(arr[indx[i]]<=a) break;
          indx[i+1]=indx[i];
        }
        indx[i+1]=indxt;
      }
      if(jstack==0) break;
      ir=istack[jstack--];
      l=istack[jstack--];
    }else{
      k=(l+ir) >> 1;
      SWAP(indx[k],indx[l+1]);
      if(arr[indx[l+1]]>arr[indx[ir]]){
        SWAP(indx[l+1],indx[ir])
      }
      if(arr[indx[l]]>arr[indx[ir]]){
        SWAP(indx[l],indx[ir])
      }
      if(arr[indx[l+1]]>arr[indx[l]]){
        SWAP(indx[l+1],indx[l])
      }
      i=l+1;
      j=ir;
      indxt=indx[l];
      a=arr[indxt];
      for(;;){
        do i++; while(arr[indx[i]]<a);
        do j--; while(arr[indx[j]]>a);
        if(j<i) break;
        SWAP(indx[i],indx[j])
      }
      indx[l]=indx[j];
      indx[j]=indxt;
      jstack+=2;
      if(jstack>NSTACK)
        ErrorExit("NSTACK too small in indexx.");
      if(ir-i+1>=j-l){
        istack[jstack]=ir;
        istack[jstack-1]=i;
        ir=j-1;
      }else{
        istack[jstack]=j-1;
        istack[jstack-1]=l;
        l=i;
      }
    }
  }
  delete[] istack;
}
#undef M
#undef NSTACK
#undef SWAP


//abstract class for Maps;
class map {
	protected:
		int SizeDX;				//Size for each dX;
		double *MatR;			//Upper-triangular matrix R (or alpha) in GS (see note);
		double *TmpMatR;		//Temporarily used matrix R;
		double *StoredX;		//Stored X;
		double *StoredDX;		//Stored dX;
		double *StoredMatR;		//Stored MatR;
		//
	public:
		MersenneTwister RNG;	//Random Number Generator;
		int Ver;				//Version of the program (in which the present data is obtained);
		double *pStoredX, *pStoredDX, *pStoredMatR;	//Pointers for them;
		bool isContinued;		//Simulation from saved data or not;
		bool isFromBeginning;	//Start from the beginning (right after the transient);
		bool isCheckingCVExp;	//Check exponents for covariant vectors or not;
		bool isDoingBackward;	//Do backward evolution or not;
		bool isNowBackward;		//During backward evolution or not;
		bool SkipReseeding;		//Skip making a new seed;
		bool isRecordingOutputQuantity;
		bool ExistOutputQuantity;		//true if OutputQuantity is recorded in forward process;
		bool ExistOutputQuantityBack;	//true if OutputQuantity is recorded in backward process;
		bool isNowRecordingOutputQuantity;		//true at the moment OutputQuantity is recorded;
		bool isRecordingQuantityHist;	//true if QuantityHist is recorded;
		bool isQuantityHistBack;	//true if QuantityHist is recorded in backward process;
		bool isRecordingPowerSpectra, isRecordingPowSpecDyn;
		bool isRecordingTsCumulants;
		int BackwardMode;		//Backward mode;
		bool isDividingLoop;	//whether trajectory is divided into short forward-backward loops;
		bool isDuringCVExp;		//true if df called by CheckCVExp;
		const int MapID;		//ID number of map;
		const string MapName;	//Name of map;
		const int ParamNum;		//Number of parameters;
		const int LocalDOF;		//degree of freedom per each site;
		int Dim;				//Number of spacial dimensions;
		int Size, RealSize;		//Number of local elements;
		int SizeMatR;			//Size of each MatR;
		int NumLyap;			//Number of computed Lyapunov exponents & vectors;
		int t;					//Time;
		int ThermalizationTime;	//Transient time for thermalization;
		int TransientTime;		//Transient time for vectors (both forward and backward);
		int RecordPeriod;		//Record period;
		int SavedRecordPeriod;	//Record period of saved data;
		int TransientEndTime;	//Time until transient;
		int RecordEndTime;		//Time until which data are recorded;
		int StoreEndTime;		//Time until which data are stored;
		int LoopPeriod;			//Period for each forward-backward loop;
		int LoopByte;			//Byte for each loop;
		int NumLoop;			//Number of forward-backward loops;
		int NumRecLoop;			//Number of forward-backward loops excluding final transient;
		int NumTransLoop;		//Number of forward-backward loops for transient;
		int LoopInd;			//Current loop;
		int SavedLoopInd;		//Saved index of loop;
		int *LoopTimes;			//Times of loop junctions;
		int RecordPerFile;		//Number of records per file;
		int CurrTmpFileID;		//ID of current temporary file;
		int OutputQuantityNum;		//Size of array OutputQuantity for forward process;
		int OutputQuantityNumBack;	//Size of array OutputQuantity for backward process;
		int QuantityHistNum, QuantityHistBinNum;	//Number of quantity histograms and bins of each histogram;
		int QuantityHistCt, QuantityHistRecInt;		//Counter and interval for quantity histograms;
		int QuantityHistSupplNum;	//Number of QuantityHistSuppl;
		unsigned int *QuantityHist;
		int *QuantityHistSupplDistID;
		ofstream OutFileQuantityHist;				//File stream for output QuantityHist;
		double MaxMem;			//Maximum occupied memory;
		double TimeStepLength;	//Length of each time step;
		double SignDynamics;	//+1 for forward dynamics, -1 for backward dynamics (used to calculate Lyap. exps.);
		double QuantityHistMin, QuantityHistMax, InvQuantityHistBinWidth;
		fstream OutFileTmp;		//Temporary file for data storing;
		string OutFileName;		//File name;
		ifstream *InFile;		//Stream for reading save data;
		int LyapRecCount;		//Recording times of Lyapunov exponents;
		int CVExpRecCount;		//Recording times of exponents for covariant vectors;
		int GSInterval;			//Time interval between orthonormalizing vectors;
		int NextGSTime;			//Next time of orthonormalization;
		unsigned long Seed;		//Seed of random number generator;
		double *X;				//Local dynamical variables;
		double *dX;				//Lyapunov vectors (Gram-Schmidt);
		double *CV;				//Covariant vectors;
		double *MatC;			//Upper-triangular matrix C for backward vectors;
		double *MeanX;			//Mean field of x;
		double *LocalExpRate;	//Local expansion rate;
		double *LyapExp;		//Lyapunov exponents;
		double *CVExp;			//exponents for covariant vectors;
		double *LyapExpStd;		//Standard deviation of local Lyapunov exponents;
		double *CVExpStd;		//Standard deviation of local exponents for covariant vectors;
		double *LocalPtcpRatio;	//Instanteneous value of paticipation ratio;
		double *GSPtcpRatio;	//Participation ratio for GS vectors;
		double *CVPtcpRatio;	//Participation ratio for covariant vectors;
		double *Param;			//Parameters;
		double *OutputQuantity;	//Some quantity to output;
		double *PowerSpectra, *PSFFTdata;	//Power spectra & temporary values for FFT;
		double *PowSpecDyn;
		double *QuantityHistSuppl;	//Averaged quantities associated with QuantityHist;
		double *TsCumulants;	//Cumulants for mean-field time series;
		string *ParamName;		//Name of parameters;
		//
		map(int TmpMapID, string TmpMapName, int TmpParamNum,
			int TmpLocalDOF) :
			RNG(0), MapID(TmpMapID), MapName(TmpMapName), ParamNum(TmpParamNum),
			LocalDOF(TmpLocalDOF) {
			//Constructor;
			//Make random number generator and set values of constants;
			//Allocate memory for parameters;
			Param = new double[ParamNum];
			ParamName = new string[ParamNum];
		};
		~map() {
			//Destructor;
			//Free memory for arrays;
			delete[] X;
			delete[] MeanX;
			delete[] dX;
			delete[] CV;
			delete[] MatR, TmpMatR;
			delete[] MatC;
			delete[] LocalExpRate, LocalPtcpRatio;
			delete[] LyapExp, LyapExpStd;
			delete[] CVExp, CVExpStd;
			delete[] GSPtcpRatio;
			delete[] CVPtcpRatio;
			switch (BackwardMode) {
				case 1:
				case 4:
				case 5:
					delete[] StoredX;
					delete[] StoredDX;
					delete[] StoredMatR;
					break;
				case 3:
					delete[] StoredX;
					break;
			}
			delete[] Param;
			delete[] ParamName;
			if (isDividingLoop) delete[] LoopTimes;
			if (isContinued) delete InFile;
			if (isRecordingQuantityHist) delete[] QuantityHist;
			if (isRecordingQuantityHist && QuantityHistSupplNum > 0) delete[] QuantityHistSuppl, QuantityHistSupplDistID;
			if (isRecordingPowerSpectra) delete[] PowerSpectra;
			if (isRecordingPowSpecDyn) delete[] PowSpecDyn;
			if (isRecordingPowerSpectra || isRecordingPowSpecDyn) delete[] PSFFTdata;
			//Close temporary file;
			if (isDividingLoop) {
				int i;
				char buf[10];
				OutFileTmp.close();
				if ((BackwardMode == 1 && !isContinued) || (BackwardMode == 3) || (BackwardMode == 4)) {
					for (i=0;i<(NumLoop-1) / RecordPerFile + 1;i++) {
						if (i==0) strcpy(buf,"\0"); else sprintf(buf,"%d",i);
						if (isFileExisting((OutFileName + ".var" + buf).c_str()))
							remove((OutFileName + ".var" + buf).c_str());
					}
				}
			}
			//Close data file;
			if (isRecordingQuantityHist && QuantityHistRecInt > 0) OutFileQuantityHist.close();
		}
		virtual void InitPowerSpectra() {
			//Initialize power spectra;
			int i, MaxInd;
			if (isRecordingPowerSpectra) {
				MaxInd = NumLyap * LocalDOF * (RealSize / 2 + 1);
				PowerSpectra = new double[MaxInd];
				for (i=0;i<MaxInd;i++) PowerSpectra[i] = 0.0;
			}
			if (isRecordingPowSpecDyn) {
				MaxInd = LocalDOF * (RealSize / 2 + 1);
				PowSpecDyn = new double[MaxInd];
				for (i=0;i<MaxInd;i++) PowSpecDyn[i] = 0.0;
			}
			if (isRecordingPowerSpectra || isRecordingPowSpecDyn) {
				PSFFTdata = new double[RealSize];
			}
		}
		virtual void Init() {
			//Initialize map;
			int i;
			//Set variables;
			Size = LocalDOF;
			for (i=0;i<Dim;i++) Size *= RealSize;		//Size = RealSize^Dim * LocalDOF;
			SizeDX = Size * NumLyap;
			SizeMatR = NumLyap * ( NumLyap + 1 ) / 2;
			OutputQuantityNum = 0;
			OutputQuantityNumBack = 0;
			isNowRecordingOutputQuantity = true;
			isDuringCVExp = false;
			//Set transient time;
			if (!isContinued) ThermalizationTime = (ThermalizationTime > 0 ? ThermalizationTime * Size : abs(ThermalizationTime));
			if (!isContinued) TransientTime = (TransientTime > 0 ? TransientTime * Size : abs(TransientTime));
			TransientEndTime = ThermalizationTime + TransientTime;
			RecordEndTime = TransientEndTime + RecordPeriod;
			StoreEndTime = RecordEndTime + TransientTime;
			//Check if GSInterval is a divisor of TransientTime;
			if (TransientTime % GSInterval != 0)
				ErrorExit("Improper GS interval.");
			//Calculate period of forward-backward loop (see note P.3);
			isDoingBackward = (BackwardMode == 1) || (BackwardMode == 4) || (BackwardMode == 5);
			if (BackwardMode == 1 || BackwardMode == 2 || BackwardMode == 5) {
				LoopPeriod = static_cast<int>(MaxMem / static_cast<double>((Size + NumLyap/2) * NumLyap * sizeof(double)));
				LoopPeriod = (LoopPeriod / GSInterval) * GSInterval;
				if (BackwardMode == 2 && RecordPeriod % LoopPeriod != 0) {
					//Set RecordPeriod a multiple of LoopPeriod;
					RecordPeriod = ((RecordPeriod / LoopPeriod) + 1) * LoopPeriod;
					cout << "Recording period has been changed to " << RecordPeriod << endl;
					RecordEndTime = TransientEndTime + RecordPeriod;
					StoreEndTime = RecordEndTime + TransientTime;
				}
				if (LoopPeriod == 0) {
					//Error;
					ErrorExit("Larger memory is required.");
				} else if (BackwardMode == 1 && LoopPeriod > RecordPeriod + TransientTime) {
					//Without division;
					isDividingLoop = false;
					LoopPeriod = RecordPeriod;
					cout << "Trajectory not devided." << endl;
				} else {
					//With division;
					isDividingLoop = true;
					NumRecLoop = (RecordPeriod-1) / LoopPeriod + 1;
					if (BackwardMode == 1 || BackwardMode == 5) NumLoop = NumRecLoop + (TransientTime-1) / LoopPeriod + 1;
						else NumLoop = NumRecLoop;
					LoopTimes = new int[NumLoop+1];
					for (i=0;i<NumRecLoop;i++) LoopTimes[i] = TransientEndTime + i * LoopPeriod;
					for (i=NumRecLoop;i<NumLoop;i++) LoopTimes[i] = RecordEndTime + (i-NumRecLoop) * LoopPeriod;
					if (BackwardMode == 1 || BackwardMode == 5) LoopTimes[NumLoop] = StoreEndTime;
						else LoopTimes[NumLoop] = RecordEndTime;
					LoopByte = sizeof(double) * (Size + SizeDX);
					cout << "Trajectory divided into " << NumLoop << " loops." << endl;
				}
			} else if (BackwardMode == 3) {
				LoopPeriod = static_cast<int>(MaxMem / static_cast<double>(Size * sizeof(double)));
				LoopPeriod = (LoopPeriod / GSInterval) * GSInterval;
				if (LoopPeriod == 0) {
					//Error;
					ErrorExit("Larger memory is required.");
				} else if (LoopPeriod > RecordPeriod + TransientTime) {
					//Without division;
					isDividingLoop = false;
					LoopPeriod = RecordPeriod;
					cout << "Trajectory not devided." << endl;
				} else {
					//With division;
					isDividingLoop = true;
					NumRecLoop = (RecordPeriod-1) / LoopPeriod + 1;
					NumLoop = NumRecLoop + (TransientTime-1) / LoopPeriod + 1;
					LoopTimes = new int[NumLoop+1];
					for (i=0;i<NumRecLoop;i++) LoopTimes[i] = TransientEndTime + i * LoopPeriod;
					for (i=NumRecLoop;i<NumLoop;i++) LoopTimes[i] = RecordEndTime + (i-NumRecLoop) * LoopPeriod;
					LoopTimes[NumLoop] = StoreEndTime;
					LoopByte = sizeof(double) * Size;
					cout << "Trajectory divided into " << NumLoop << " loops." << endl;
				}
			} else if (BackwardMode == 4) {
				LoopPeriod = static_cast<int>(MaxMem / static_cast<double>((Size + NumLyap/2) * NumLyap * sizeof(double)));
				LoopPeriod = (LoopPeriod / GSInterval) * GSInterval;
				if (LoopPeriod == 0) {
					//Error;
					ErrorExit("Larger memory is required.");
				} else if (LoopPeriod > RecordPeriod + TransientTime) {
					//Without division;
					isDividingLoop = false;
					LoopPeriod = RecordPeriod;
					cout << "Trajectory not devided." << endl;
				} else {
					//With division;
					isDividingLoop = true;
					NumTransLoop = (TransientTime-1) / LoopPeriod + 1;
					NumRecLoop = NumTransLoop + (RecordPeriod-1) / LoopPeriod + 1;
					NumLoop = NumRecLoop + (TransientTime-1) / LoopPeriod + 1;
					LoopTimes = new int[NumLoop+1];
					for (i=0;i<NumTransLoop;i++) LoopTimes[i] = ThermalizationTime + i*LoopPeriod;
					for (i=NumTransLoop;i<NumRecLoop;i++) LoopTimes[i] = TransientEndTime + (i-NumTransLoop) * LoopPeriod;
					for (i=NumRecLoop;i<NumLoop;i++) LoopTimes[i] = RecordEndTime + (i-NumRecLoop) * LoopPeriod;
					LoopTimes[NumLoop] = StoreEndTime;
					LoopByte = sizeof(double) * (Size + SizeDX);
					cout << "Trajectory divided into " << NumLoop << " loops." << endl;
				}
			} else {
				isDividingLoop = false;
			}
			//Set configurations about file io;
			if (isDividingLoop) {
				RecordPerFile = MAXFILESIZE / (LoopByte / (1024*1024) + 1);
				cout << "Required disk space = " << static_cast<double>(LoopByte) * static_cast<double>(NumLoop) / (1024.0*1024.0) << " MB" << endl;
				cout << (NumLoop-1) / RecordPerFile + 1 << " temporary file(s) will be created." << endl;
			}
			//Set seed for random number generator;
			RNG.init_genrand(Seed);
			//Allocate memory for dynamical variables & Lyapunov vectors;
			X = new double[Size];
			MeanX = new double[LocalDOF];
			dX = new double[SizeDX];
			CV = new double[SizeDX];
			MatR = new double[SizeMatR];
			TmpMatR = new double[NumLyap];
			MatC = new double[SizeMatR];
			LocalExpRate = new double[NumLyap];
			LyapExp = new double[NumLyap];
			CVExp = new double[NumLyap];
			LyapExpStd = new double[NumLyap];
			CVExpStd = new double[NumLyap];
			LocalPtcpRatio = new double[NumLyap];
			GSPtcpRatio = new double[NumLyap];
			CVPtcpRatio = new double[NumLyap];
			try{
				switch (BackwardMode) {
					case 1:
					case 5:
						if (isDividingLoop) {
							StoredX = new double[LoopPeriod * Size];
							StoredDX = new double[LoopPeriod * SizeDX];
							StoredMatR = new double[LoopPeriod * SizeMatR];
						} else {
							StoredX = new double[RecordPeriod * Size];
							StoredDX = new double[RecordPeriod * SizeDX];
							StoredMatR = new double[(RecordPeriod + TransientTime / GSInterval) * SizeMatR];
						}
						//Set pointers;
						pStoredX = StoredX; pStoredDX = StoredDX; pStoredMatR = StoredMatR;
						break;
					case 3:
						if (isDividingLoop) StoredX = new double[(LoopPeriod + 1) * Size];
							else StoredX = new double[(RecordPeriod + TransientTime + 1) * Size];
						//Set pointers;
						pStoredX = StoredX;
						break;
					case 4:
						if (isDividingLoop) {
							StoredX = new double[(LoopPeriod + 2) * Size];
							StoredDX = new double[LoopPeriod * SizeDX];
							StoredMatR = new double[LoopPeriod * SizeMatR];
						} else {
							StoredX = new double[(RecordPeriod + 2 * TransientTime + 1) * Size];
							StoredDX = new double[RecordPeriod * SizeDX];
							StoredMatR = new double[(RecordPeriod + TransientTime / GSInterval) * SizeMatR];
						}
						//Set pointers;
						pStoredX = StoredX; pStoredDX = StoredDX; pStoredMatR = StoredMatR;
						break;
				}
			} catch(bad_alloc) {
				ErrorExit("Failed to allocate memory.");
			}
			//Allocate memory for power spectra;
			InitPowerSpectra();
		}
		void ResetPStored() {
			//Reset pointers for stored data;
			pStoredX = StoredX; pStoredDX = StoredDX; pStoredMatR = StoredMatR;
		}
		void CheckTmpFileID(int NewID) {
			//Check ID of currently opened temporary file and change it if necessary;
			if (NewID == CurrTmpFileID) return;
			//Change temporary file;
			char buf[10];
			OutFileTmp.close();
			if (NewID==0) strcpy(buf,"\0"); else sprintf(buf,"%d",NewID);
			OutFileTmp.open((OutFileName + ".var" + buf).c_str(),ios::in | ios::out | ios::binary | (isFileExisting((OutFileName + ".var" + buf).c_str()) ? ios::ate : ios::trunc));
			CurrTmpFileID = NewID;
		}
		virtual void WriteVar() {
			//Write variables at junctions of forward-backward loops;
			CheckTmpFileID(LoopInd/RecordPerFile);
			OutFileTmp.seekp((LoopInd - CurrTmpFileID * RecordPerFile) * LoopByte, ios::beg);
			OutFileTmp.write(reinterpret_cast<char*>(X),Size*sizeof(double));
			OutFileTmp.write(reinterpret_cast<char*>(dX),SizeDX*sizeof(double));
			if (!OutFileTmp.good()) ErrorExit("IO Error!");
			LoopInd++;
		}
		virtual void ReadVar() {
			//Read variables at junctions of forward-backward loops;
			int i,j;
			LoopInd--;
			t = LoopTimes[LoopInd];
			CheckTmpFileID(LoopInd/RecordPerFile);
			OutFileTmp.seekg((LoopInd - CurrTmpFileID * RecordPerFile) * LoopByte, ios::beg);
			OutFileTmp.read(reinterpret_cast<char*>(X),Size*sizeof(double));
			OutFileTmp.read(reinterpret_cast<char*>(dX),SizeDX*sizeof(double));
		}
		virtual void WriteVarX() {
			//Write variable X at junctions of forward-backward loops;
			CheckTmpFileID(LoopInd/RecordPerFile);
			OutFileTmp.seekp((LoopInd - CurrTmpFileID * RecordPerFile) * LoopByte, ios::beg);
			OutFileTmp.write(reinterpret_cast<char*>(X),Size*sizeof(double));
			if (!OutFileTmp.good()) ErrorExit("IO Error!");
			LoopInd++;
		}
		virtual void ReadVarX() {
			//Read variables at junctions of forward-backward loops;
			int i,j;
			LoopInd--;
			t = LoopTimes[LoopInd];
			CheckTmpFileID(LoopInd/RecordPerFile);
			OutFileTmp.seekg((LoopInd - CurrTmpFileID * RecordPerFile) * LoopByte, ios::beg);
			OutFileTmp.read(reinterpret_cast<char*>(X),Size*sizeof(double));
		}
		virtual void AppendVarDX() {
			//Append variables dX at junctions of forward-backward loops;
			//Do not change LoopInd;
			CheckTmpFileID(LoopInd/RecordPerFile);
			OutFileTmp.seekp((LoopInd - CurrTmpFileID * RecordPerFile) * LoopByte + sizeof(double) * Size, ios::beg);
			OutFileTmp.write(reinterpret_cast<char*>(dX),SizeDX*sizeof(double));
			if (!OutFileTmp.good()) ErrorExit("IO Error!");
		}
		virtual void ReadVarNoDX() {
			//Read variables except dX at junctions of forward-backward loops;
			int i,j;
			LoopInd--;
			t = LoopTimes[LoopInd];
			CheckTmpFileID(LoopInd/RecordPerFile);
			OutFileTmp.seekg((LoopInd - CurrTmpFileID * RecordPerFile) * LoopByte, ios::beg);
			OutFileTmp.read(reinterpret_cast<char*>(X),Size*sizeof(double));
			OutFileTmp.seekg(SizeDX*sizeof(double), ios::cur);
		}
		virtual void CalcMeanField() {
			//Calculate mean field;
			int i,j;
			for (i=0;i<LocalDOF;i++) {
				MeanX[i] = 0.0;
				for (j=i;j<Size;j+=LocalDOF) MeanX[i] += X[j];
				MeanX[i] /= static_cast<double>(Size / LocalDOF);
			}
		}
		void CalcTsCumulants() {
			//Calculate cumulants of mean-field;
			int i,j;
			double x;
			for (i=0;i<LocalDOF;i++) {
				x = 1.0;
				for (j=0;j<4;j++) {
					x *= MeanX[i];
					TsCumulants[4*i+j] += x;
				}
			}
		}
		virtual void SetInitialX() {
			//Set initial conditions;
			int i,j;
			double *pMdX;
			//Random initial conditions for X are set by child class;
			//for (i=0;i<Size;i++) X[i] = MinX + ( MaxX - MinX ) * RNG.genrand_real1();
			for (i=0;i<SizeDX;i++) dX[i] = 2.0 * RNG.genrand_real1() - 1.0;
		}
		virtual void SetInitialVar() {
			//Other initializations;
			int i;
			double a;
			if (!isContinued) {
				//Start from beginning;
				GramSchmidt();
				//other initializations;
				LoopInd = 0;
				t = 0;
				SkipReseeding = false;
				LyapRecCount = 0; CVExpRecCount = 0;
				for (i=0;i<NumLyap;i++) {
					LyapExp[i] = 0.0;
					CVExp[i] = 0.0;
					LyapExpStd[i] = 0.0;
					CVExpStd[i] = 0.0;
					GSPtcpRatio[i] = 0.0;
					CVPtcpRatio[i] = 0.0;
				}
			} else {
				//Continue from saved data;
				LoopInd = (!isFromBeginning ? SavedLoopInd : 1);
				ReadVar();
				SkipReseeding = true;
				//other initializations;
				LyapRecCount = (!isFromBeginning ? SavedRecordPeriod : 0);
				CVExpRecCount = 0;
				for (i=0;i<NumLyap;i++) {
					CVExp[i] = 0.0;
					CVExpStd[i] = 0.0;
					CVPtcpRatio[i] = 0.0;
				}
				if (!isFromBeginning) {
					InFile->read(reinterpret_cast<char*>(LyapExp), NumLyap * sizeof(double));
					InFile->read(reinterpret_cast<char*>(GSPtcpRatio), NumLyap * sizeof(double));
					if (Ver <= -2) {
						InFile->read(reinterpret_cast<char*>(LyapExpStd), NumLyap * sizeof(double));
					} else {
						//Non-sense, but to avoid errors;
						for (i=0;i<NumLyap;i++) LyapExpStd[i] = 0.0;
					}
					if (Ver == -3) InFile->ignore(3 * NumLyap * sizeof(double));
					InFile->close();
				} else {
					for (i=0;i<NumLyap;i++) {
						LyapExp[i] = 0.0;
						LyapExpStd[i] = 0.0;
						GSPtcpRatio[i] = 0.0;
					}
				}
			}
			//other initializations;
			isNowBackward = false;
		}
		virtual void SingleDF(const double *pX, double *ppdX) {};
			//Time evolution of a single vector;
		virtual void SingleBackDF(const double *pX, double *ppdX) {};
			//Backward time evolution of a signle vector;
		virtual void df(double *pX, double *pdX) {
			//Time evolution of vectors;
			int i;
			//Compute dx;
#if OMPSW==1
			#pragma omp parallel private(i)
			#pragma omp for
#endif
			for (i=0;i<NumLyap;i++) SingleDF(pX, pdX + i * Size);
		}
		virtual void BackDF(double *pdX) {
			//Backward time evolution of vectors;
			//Values of X are restored here;
			//
			//As usual;
			t--;
		}
		virtual void f() {
			//Time evolution;
			//As usual;
			t++;
		}
		virtual double VectorNorm(double *pdX) {
			//return norm of a given vector;
			int i;
			double Norm = 0.0;
#if OMPSW==2 && OMPRSW
			#pragma omp parallel private(i)
			#pragma omp for reduction(+:Norm)
#endif
			for (i=0;i<Size;i++) Norm += pdX[i] * pdX[i];
			return sqrt(Norm);
		}
		virtual double VectorDotProduct(double *pdX1, double *pdX2) {
			//return dot product of a given pair of vectors;
			int i;
			double Norm = 0.0;
			for (i=0;i<Size;i++) Norm += pdX1[i] * pdX2[i];
			return Norm;
		}
#if OMPSW
		virtual void GramSchmidt() {
			//Do Gram-Schmidt orthonormalization (openMP multithread version);
			int i,j,k;
			double Norm,Norm2;
			double *pdX,*pdX2;	//vector y' and x, respectively, in my note;
			double *pMatR=MatR;	//pointer for MatR;
			#pragma omp parallel private(i,j,k,pdX2,Norm2)
			for (i=0;i<NumLyap;i++) {
				#pragma omp single
				{
				//Set pointer;
				pdX = dX + i * Size;
				//Calculate Norm;
				Norm = VectorNorm(pdX);
				*pMatR = Norm; pMatR++;
				}
				//Normalization;
				#pragma omp for
				for (k=0;k<Size;k++) pdX[k] /= Norm;
				//Orthogonalization;
				#pragma omp for
				for (j=i+1;j<NumLyap;j++) {
					//Set pointer;
					pdX2 = dX + j * Size;
					//Calculate dot products;
					Norm2 = VectorDotProduct(pdX, pdX2);
					pMatR[j-i-1] = Norm2;
					//Subtract vector x projection on y' from x;
					for (k=0;k<Size;k++) pdX2[k] -= Norm2 * pdX[k];
				}
				#pragma omp single
				{
				pMatR += NumLyap - 1 - i;
				}
			}
		}
#else
		virtual void GramSchmidt() {
			//Do Gram-Schmidt orthonormalization;
			int i,j,k;
			double Norm;
			double *pdX,*pdX2;	//vector y' and x, respectively, in my note;
			double *pMatR=MatR;	//pointer for MatR;
			for (i=0;i<NumLyap;i++) {
				//Set pointer;
				pdX = dX + i * Size;
				//Calculate Norm;
				Norm = VectorNorm(pdX);
				*pMatR = Norm; pMatR++;
				//Normalization;
				for (k=0;k<Size;k++) pdX[k] /= Norm;
				//Orthogonalization;
				for (j=i+1;j<NumLyap;j++) {
					//Set pointer;
					pdX2 = dX + j * Size;
					//Calculate dot products;
					Norm = VectorDotProduct(pdX, pdX2);
					*pMatR = Norm; pMatR++;
					//Subtract vector x projection on y' from x;
					for (k=0;k<Size;k++) pdX2[k] -= Norm * pdX[k];
				}
			}
		}
#endif
		void CalcLyap() {
			//Calculate Lyapunov exponents;
			int i;
			double *pMatR;
			pMatR = MatR;
			LyapRecCount++;
			for (i=0;i<NumLyap;i++) {
				LocalExpRate[i] = log(*pMatR);
				LyapExp[i] += LocalExpRate[i];
				LyapExpStd[i] += LocalExpRate[i] * LocalExpRate[i];
				pMatR += NumLyap - i;
			}
		}
		virtual void CalcPtcpRatio(double *pdX, double *pPtcpRatio) {
			//Calculate participation ratio;
			int i,j,k;
			double a,b;
			double *ppdX;
			for (i=0;i<NumLyap;i++) {
				ppdX = pdX + i*Size;
				LocalPtcpRatio[i] = 0.0;
				for (j=0;j<Size;j+=LocalDOF) {
					a = 0.0;
					for (k=0;k<LocalDOF;k++) {
						b = ppdX[j+k] * ppdX[j+k];
						a += b;
					}
					LocalPtcpRatio[i] += a*a;
				}
				pPtcpRatio[i] += LocalPtcpRatio[i];
			}
		}
		void StoreX() {
			//Store X;
			int i;
			for (i=0;i<Size;i++) pStoredX[i] = X[i];
			pStoredX += Size;
		}
		void StoreDX() {
			//Store dX;
			int i;
			for (i=0;i<SizeDX;i++) pStoredDX[i] = dX[i];
			pStoredDX += SizeDX;
		}
		virtual void StoreXdX() {
			//Store X and dX;
			int i;
			for (i=0;i<Size;i++) pStoredX[i] = X[i];
			for (i=0;i<SizeDX;i++) pStoredDX[i] = dX[i];
			pStoredX += Size;
			pStoredDX += SizeDX;
		}
		virtual void StoreMatR() {
			//Store MatR;
			int i;
			for (i=0;i<SizeMatR;i++) pStoredMatR[i] = MatR[i];
			pStoredMatR += SizeMatR;
		}
		virtual void SetInitialMatC() {
			//Set initial conditions for backward vectors;
			int i,k;
			double Sum;
			double *pMatC;	//pointer for MatC;
			for (i=0;i<SizeMatR;i++) MatC[i] = 2.0 * RNG.genrand_real1() - 1.0;
			for (k=0;k<NumLyap;k++) {
				//Set pointer;
				pMatC = MatC + k*(k+1)/2;
				//Normalization;
				Sum = 0.0;
				for (i=0;i<=k;i++) Sum += pMatC[i] * pMatC[i];
				Sum = 1.0 / sqrt(Sum);
				for (i=0;i<=k;i++) pMatC[i] *= Sum;
			}
		}
		virtual void BackwardEvolution() {
			//Evolve MatC backward (compatible with OpenMP);
			//See note and "Numerical Recipes in C" P.41 for algorithm;
			int i,j,k;
			double Sum;
			double *pMatC;	//pointer for MatC;
			//Set pointer;
			pStoredMatR -= SizeMatR;
			//Solve C_m = R_m C_{m-1} (see note);
#if OMPSW
			#pragma omp parallel private(i,j,k,pMatC,Sum)
			#pragma omp for schedule(dynamic)
#endif
			for (k=0;k<NumLyap;k++) {
				//Set pointer;
				pMatC = MatC + k*(k+1)/2;
				//Main part;
				for (i=k;i>=0;i--) {
					//Sum = 0.0;
					//for (j=i+1;j<=k;j++) Sum += pStoredMatR[i*NumLyap + j - i*(i+1)/2] * pMatC[j];
					//pMatC[i] = (pMatC[i] - Sum) / pStoredMatR[i*NumLyap + i - i*(i+1)/2];
					for (j=i+1;j<=k;j++) pMatC[i] -= pStoredMatR[i*NumLyap + j - i*(i+1)/2] * pMatC[j];
					pMatC[i] /= pStoredMatR[i*NumLyap + i - i*(i+1)/2];
				}
				//Normalization;
				Sum = 0.0;
				for (i=0;i<=k;i++) Sum += pMatC[i] * pMatC[i];
				Sum = 1.0 / sqrt(Sum);
				for (i=0;i<=k;i++) pMatC[i] *= Sum;
			}
		}
		virtual void RestoreX() {
			//Set pointers of X to restore them;
			pStoredX -= Size;
		}
		virtual void ReRestoreX() {
			//Set pointers of X to restore them (for 2nd forward evolution);
			pStoredX += Size;
		}
		void RestoreDX() {
			//Set pointers of dX to restore them;
			pStoredDX -= SizeDX;
		}
		virtual void RestoreXdX() {
			//Set pointers of X and dX to restore them;
			pStoredX -= Size; pStoredDX -= SizeDX;
		}
		virtual void CovariantVectors() {
			//Calculate covariant vectors (OpenMP multithread version);
			int i,j,k;
			double *pCV, *pMatC, *ppStoredDX;
			//Calculate covariant vectors;
			// v_m^j = \sum_{i=0}^j C_m^{ij} * e_m^i;
#if OMPSW
			#pragma omp parallel private(i,j,k,ppStoredDX,pCV,pMatC)
			{
				#pragma omp for
#endif
				for (i=0;i<SizeDX;i++) CV[i] = 0.0;
#if OMPSW
				#pragma omp for schedule(dynamic)
#endif
				for (j=0;j<NumLyap;j++) {
					pCV = CV + j * Size;
					pMatC = MatC + j*(j+1)/2;
					ppStoredDX = pStoredDX;
					for (i=0;i<=j;i++) {
						for (k=0;k<Size;k++) pCV[k] += pMatC[i] * ppStoredDX[k];
						ppStoredDX += Size;
					}
				}
#if OMPSW
			}
#endif
		}
		virtual void CheckCVExp() {
			//Check exponents for covariant vectors;
			// (Caution! It changes values of vectors);
			int i,j;
			double Norm;
			double *pCV;
			//Progress of counter;
			CVExpRecCount++;
			//Do nothing if isCheckingCVExp = false;
			if (!isCheckingCVExp) return;
			//Time evolution of vectors;
			isDuringCVExp = true;
			df(pStoredX, CV);
			isDuringCVExp = false;
			//Calculate expansion rate;
			for (i=0;i<NumLyap;i++) {
				pCV = CV + i*Size;
				Norm = VectorNorm(pCV);
				LocalExpRate[i] = log(Norm);
				CVExp[i] += LocalExpRate[i];
				CVExpStd[i] += LocalExpRate[i] * LocalExpRate[i];
			}
		}
		void CheckBackCVExp() {
			//Check exponents for covariant vectors obtained by 2nd forward evolution;
			// (Caution! It changes values of vectors);
			int i,j;
			double Norm;
			double *pCV;
			//Progress of counter;
			CVExpRecCount++;
			//Do nothing if isCheckingCVExp = false;
			if (!isCheckingCVExp) return;
			//Time evolution of vectors;
			BackDF(CV);
			//Put back time and the pointer for X;
			t++;
			ReRestoreX();
			//Calculate expansion rate;
			for (i=0;i<NumLyap;i++) {
				pCV = CV + i*Size;
				Norm = VectorNorm(pCV);
				LocalExpRate[i] = log(Norm);
				CVExp[i] += LocalExpRate[i];
				CVExpStd[i] += LocalExpRate[i] * LocalExpRate[i];
			}
		}
		void CalcResults() {
			//Calculate results;
			int i;
			//Lyapunov exponents & Participation ratio;
			for (i=0;i<NumLyap;i++) {
				LyapExp[i] /= (static_cast<double>(RecordPeriod) * TimeStepLength * SignDynamics);
				CVExp[i] /= (static_cast<double>(RecordPeriod) * TimeStepLength * SignDynamics);
				GSPtcpRatio[i] /= static_cast<double>(RecordPeriod);
				CVPtcpRatio[i] /= static_cast<double>(RecordPeriod);
				LyapExpStd[i] = LyapExpStd[i] / (static_cast<double>(RecordPeriod) * TimeStepLength * TimeStepLength) - LyapExp[i] * LyapExp[i];
				CVExpStd[i] = CVExpStd[i] / (static_cast<double>(RecordPeriod) * TimeStepLength * TimeStepLength) - CVExp[i] * CVExp[i];
			}
		}
		virtual void CalcOutputQuantity() {
			//Calculate OutputQuantity;
		}
		virtual void CalcOutputQuantityBack() {
			//Calculate OutputQuantity;
		}
		virtual void AccQuantityHist() {
			//Accumulate histogram for map-specific quantity;
			//
			//This must come after subroutine for each map;
			int i;
			//
			QuantityHistCt++;
			//Save;
			if (QuantityHistCt == QuantityHistRecInt) {
			//Record and reset histogram of map-specific quantity;
				//OutFile.open((OutFileName+"-quantityhistts.dat").c_str(), ios::out | ios::app);
				OutFileQuantityHist << setprecision(9) << static_cast<double>(t) * TimeStepLength << setprecision(6) << ' ';
				for (i=0;i<QuantityHistNum * QuantityHistBinNum;i++) OutFileQuantityHist << QuantityHist[i] << ' ';
				OutFileQuantityHist << endl;
				//OutFile.close();
				//
				for (i=0;i<QuantityHistNum * QuantityHistBinNum;i++) QuantityHist[i] = 0;
				QuantityHistCt = 0;
			}
		}
		virtual void InitQuantityHist() {
			//Initialize histogram for map-specific quantity;
			int i, MaxInd;
			if (isRecordingQuantityHist) {
				InvQuantityHistBinWidth = static_cast<double>(QuantityHistBinNum) / (QuantityHistMax - QuantityHistMin);
				MaxInd = QuantityHistNum * QuantityHistBinNum;
				QuantityHist = new unsigned int[MaxInd];
				for (i=0;i<MaxInd;i++) QuantityHist[i] = 0;
				if (QuantityHistSupplNum > 0) {
					MaxInd = QuantityHistSupplNum * QuantityHistBinNum;
					QuantityHistSuppl = new double[MaxInd];
					for (i=0;i<MaxInd;i++) QuantityHistSuppl[i] = 0.0;
					QuantityHistSupplDistID = new int[QuantityHistSupplNum];
					for (i=0;i<QuantityHistSupplNum;i++) QuantityHistSupplDistID[i] = 0;
				}
			}
			QuantityHistCt = 0;
			if (QuantityHistRecInt > 0) OutFileQuantityHist.open((OutFileName+"-quantityhistts.dat").c_str(), ios::out | ios::trunc);
		}
		virtual void AccPowerSpectra() {
			//Calculate power spectra of covariant vectors;
			int i,j,k;
			double Norm;
			double *pPowerSpectra, *pCV;
			for (i=0;i<LocalDOF;i++) {
				for (j=0;j<NumLyap;j++) {
					pPowerSpectra = PowerSpectra + (i * NumLyap + j) * (RealSize / 2 + 1);
					pCV = CV + j * Size;
					for (k=0;k<RealSize;k++) PSFFTdata[k] = pCV[k * LocalDOF + i];
					//Norm = 0.0;	//not normalize here since 2009/10/12
					//for (k=0;k<RealSize;k++) Norm += PSFFTdata[k] * PSFFTdata[k];
					//Norm = sqrt(Norm);
					//for (k=0;k<RealSize;k++) PSFFTdata[k] /= Norm;
					realft2(PSFFTdata-1, static_cast<unsigned long>(RealSize), -1); 
					pPowerSpectra[0] += PSFFTdata[0] * PSFFTdata[0];
					for (k=1;k<RealSize/2;k++)
						pPowerSpectra[k] += PSFFTdata[k*2] * PSFFTdata[k*2] + PSFFTdata[k*2+1] * PSFFTdata[k*2+1];
					pPowerSpectra[RealSize/2] += PSFFTdata[1] * PSFFTdata[1];
				}
			}
		}	
		virtual void AccPowSpecDyn() {
			//Calculate power spectrum of dynamics;
			int i,k;
			double Norm;
			double *pPowSpecDyn;
			for (i=0;i<LocalDOF;i++) {
				pPowSpecDyn = PowSpecDyn + i * (RealSize / 2 + 1);
				for (k=0;k<RealSize;k++) PSFFTdata[k] = X[k * LocalDOF + i];
				realft2(PSFFTdata-1, static_cast<unsigned long>(RealSize), -1); 
				pPowSpecDyn[0] += PSFFTdata[0] * PSFFTdata[0];
				for (k=1;k<RealSize/2;k++)
					pPowSpecDyn[k] += PSFFTdata[k*2] * PSFFTdata[k*2] + PSFFTdata[k*2+1] * PSFFTdata[k*2+1];
				pPowSpecDyn[RealSize/2] += PSFFTdata[1] * PSFFTdata[1];
			}
		}	
		virtual void RecordPowerSpectra() {
			//Record power spectra;
			int i,j;
			int PowerSpectraLength;
			double *pPowerSpectra;
			ofstream OutFile;
			if (isRecordingPowerSpectra) {
				PowerSpectraLength = RealSize / 2 + 1;
				OutFile.open((OutFileName+"-ps.dat").c_str());
				for (i=0;i<LocalDOF*NumLyap;i++) {
					pPowerSpectra = PowerSpectra + i * PowerSpectraLength;
					for (j=0;j<PowerSpectraLength;j++) OutFile << pPowerSpectra[j] / static_cast<double>(RecordPeriod) << ' ';
					OutFile << endl;
				}
				OutFile.close();
			}
			//Record Power Spectrum of dynamics;
			if (isRecordingPowSpecDyn) {
				PowerSpectraLength = RealSize / 2 + 1;
				OutFile.open((OutFileName+"-psdyn.dat").c_str());
				for (i=0;i<LocalDOF;i++) {
					pPowerSpectra = PowSpecDyn + i * PowerSpectraLength;
					for (j=0;j<PowerSpectraLength;j++) OutFile << pPowerSpectra[j] / static_cast<double>(RecordPeriod) << ' ';
					OutFile << endl;
				}
				OutFile.close();
			}
		}
};

//Abstract class for GCM and other general systems;
class gcm : public map {
	public:
		gcm(int TmpMapID, string TmpMapName, int TmpParamNum, int TmpLocalDOF) :
			map (TmpMapID, TmpMapName, TmpParamNum,	TmpLocalDOF) {};
		virtual void Init() {
			//Initialize map;
			Dim = 1;
			//Go to the parent;
			map::Init();
		}
};


//23: 1d Kuramoto-Sivashinsky (spectrum method, splitting method (2nd-order Adams-Moulton, Heun));
/*
  du/dt = -d^2u/dx^2 - a d^4u/dx^4 - u du/dx
  u = \sum_j [C_j exp(i k_j x)]
  k_j = 2\pi j/L
  
  Cutoff of wavelength is calculated from system resolution;
*/
class mapKSPDEp : public gcm {
	private:
		int RepNum;			//time step / dt;
		int WNCutoff;		//Cutoff of wavenumber (= largest N that satisfies 3N+1<=Size);
		int kMaxInd;		//Index for the most unstable wavenumber;
		double Pi;			//pi;
		double Kcoeff;		// = 2\pi / L;
		double FFTcoeff, IFFTcoeff;	//Coefficients to multiply results of transform;
		double kMaxLambda;	//wave length for the most unstable wavenumber;
		double *CNcoeff;		//Coefficients to multiply for Crank-Nicholson;
		double *k1, *k2, *dk1, *dk2;	//for Runge-Kutta;
		double *F, *dF;			//Nonlinear term in spectral method;
		double *u1, *u2, *du1, *du2;
		double *dudc, *ddudc;	//term evaluation at each spatial point;
		double *FFTdata;
		vector<double> PowSpecUPO;	//power spectra of UPOs;
#if FFTWSW
		int Sizec;				// = Size/2+1;
		double *FFTr;			//Array for FFT;
		complex<double> *FFTc;	//Array for FFT;
		fftw_plan FFTplan, IFFTplan;	//plan of FFTW;
		void FFT(double *v) {
			//FFT;
			int i;
			for (i=0;i<Size;i++) FFTr[i] = v[i];
			fftw_execute(FFTplan);
			for (i=0;i<=WNCutoff;i++) {
				v[2*i] = FFTcoeff * FFTc[i].real();
				v[2*i+1] = FFTcoeff * FFTc[i].imag();
			}
			for (i=2*WNCutoff+2;i<Size;i++) v[i] = 0.0;
		}
		void IFFT(double *v) {
			//Inverse FFT;
			int i;
			for (i=0;i<=WNCutoff;i++) FFTc[i] = complex<double>(v[i*2], v[i*2+1]);
			for (i=WNCutoff+1;i<Sizec;i++) FFTc[i] = 0.0;
			fftw_execute(IFFTplan);
			for (i=0;i<Size;i++) v[i] = IFFTcoeff * FFTr[i];
		}
#else
		void FFT(double *v) {
			//FFT;
			int i;
			realft2(v-1, Size, -1);
			v[0] *= FFTcoeff;
			v[1] = 0.0;
			for (i=2;i<2*WNCutoff+2;i++) v[i] *= FFTcoeff;
			for (i=2*WNCutoff+2;i<Size;i++) v[i] = 0.0;
		}
		void IFFT(double *v) {
			//Inverse FFT;
			int i;
			realft2(v-1, Size, 1);
			for (i=0;i<Size;i++) v[i] *= IFFTcoeff;
		}
#endif
	public:
		mapKSPDEp() : gcm(23, "1d Kuramoto-Sivashinsky PDE (periodic boundary)", 4, 1) {
			//Constructor;
			ParamName[0] = "dt";
			ParamName[1] = "time step / dt";
			ParamName[2] = "a";
			ParamName[3] = "Chain length";
		}
		~mapKSPDEp() {
			//Destructor;
			//Destroy arrays;
			delete[] k1, k2, dk1, dk2;
			delete[] F, dF, u1, u2, du1, du2, dudc, ddudc, CNcoeff;
			if (isRecordingOutputQuantity) delete[] OutputQuantity;
			delete[] FFTdata;
#if FFFTWSW
			//For FFTW;
			delete[] FFTr, FFTc;
			fftw_destroy_plan(FFTplan);
			fftw_destroy_plan(IFFTplan);
#endif
		}
		void Init() {
			//Initialize map;
			int i,j;
			double k, a;
			ifstream InFileTmp;
			//Go first to the parent;
			gcm::Init();
			//Set variables;
			RepNum = static_cast<int>(Param[1]);
			TimeStepLength = Param[0] * Param[1];
			Pi = 3.141592653589793;
			WNCutoff = (Size - 1) / 3;
			Kcoeff = 2.0 * Pi / Param[3];
			kMaxInd = static_cast<int>(floor(1.0 / (sqrt(2.0) * Kcoeff) + 0.5));
			kMaxLambda = Param[3] / static_cast<double>(kMaxInd);
			//Allocate memory for arrays;
			k1 = new double[Size];
			k2 = new double[Size];
			dk1 = new double[Size * omp::NumThreads];
			dk2 = new double[Size * omp::NumThreads];
			F = new double[Size];
			u1 = new double[Size];
			u2 = new double[Size];
			du1 = new double[Size];
			du2 = new double[Size];
			dF = new double[Size * omp::NumThreads];
			dudc = new double[Size * omp::NumThreads];
			ddudc = new double[Size * omp::NumThreads];
			CNcoeff = new double[Size];
			FFTdata = new double[Size];
			//Set arrays;
			for (i=0;i<=WNCutoff;i++) {
				k = Kcoeff * static_cast<double>(i);
				k *= k;
				CNcoeff[2*i] = (2.0 + Param[0] * (k - Param[2] * k * k)) / (2.0 - Param[0] * (k - Param[2] * k * k));
				CNcoeff[2*i+1] = (2.0 + Param[0] * (k - Param[2] * k * k)) / (2.0 - Param[0] * (k - Param[2] * k * k));
			}
			for (i=2*WNCutoff+2;i<Size;i++) CNcoeff[i] = 0.0;
			//For FFT;
#if FFTWSW
			Sizec = Size / 2 + 1;
			FFTr = new double[Size];
			FFTc = new complex<double>[Sizec];
			FFTcoeff = 1.0 / static_cast<double>(Size);
			IFFTcoeff = 1.0;
			FFTplan = fftw_plan_dft_r2c_1d(Size, FFTr, reinterpret_cast<fftw_complex*>(FFTc), FFTW_MEASURE);
			IFFTplan = fftw_plan_dft_c2r_1d(Size, reinterpret_cast<fftw_complex*>(FFTc), FFTr, FFTW_MEASURE);
#else
			FFTcoeff = 1.0 / static_cast<double>(Size);
			IFFTcoeff = 2.0;
#endif
			//Set OutputQuantity;
			ExistOutputQuantity = true;
			ExistOutputQuantityBack = false;
			OutputQuantityNum = 0;
			OutputQuantityNumBack = 0;
			if(isRecordingOutputQuantity && isFileExisting((OutFileName+"-upo.dat").c_str())) {
				InFileTmp.open((OutFileName+"-upo.dat").c_str(), ios::in);
				do {
					for (i=0;i<Size;i++) InFileTmp >> FFTdata[i];
					realft2(FFTdata-1, static_cast<unsigned long>(Size), -1);
					for (j=1;j<=WNCutoff;j++) {
						PowSpecUPO.push_back(FFTdata[j*2] * FFTcoeff);
						PowSpecUPO.push_back(FFTdata[j*2+1] * FFTcoeff);
					}
					OutputQuantityNum++;
				} while (!InFileTmp.eof());
				InFileTmp.close();
				OutputQuantity = new double[OutputQuantityNum];
				cout << OutputQuantityNum << " UPOs are compared." << endl;
			}
		}
		void SetInitialX() {
			//Set initial conditions;
			int i;
			double Amp, Amp1, Amp2, Amp3;
			ifstream InFile;
			//
			InFile.open((OutFileName+"-init.dat").c_str(), ios::in);
			if (InFile) {
				//Read initial condition file;
				cout << "Initial condition read from file." << endl;
				for (i=0;i<Size;i++) InFile >> X[i];
				InFile.close();
			} else {
				//Random initial condition;
				cout << "Random initial condition created." << endl;
				Amp = 0.01 / static_cast<double>(WNCutoff + 1);
				X[0] = 0.0;
				X[1] = 0.0;			//Imaginary part of zero wavelength is zero;
				Amp1 = Amp / (Kcoeff * Kcoeff * sqrt(2.0));
				Amp2 = Amp / (pow(Kcoeff, 4.0) * sqrt(2.0));
				Amp3 = sqrt(Amp / (Kcoeff * static_cast<double>(2*WNCutoff + 1) * sqrt(2.0)));
				for (i=1;i<=WNCutoff;i++) {
					Amp = MIN(MIN(Amp1 / static_cast<double>(i*i), Amp2 / static_cast<double>(i*i*i*i)), Amp3 / sqrt(static_cast<double>(i)));
					X[2*i] = Amp * (2.0 * RNG.genrand_real1() - 1.0);
					X[2*i+1] = Amp * (2.0 * RNG.genrand_real1() - 1.0);
				}
				for (i=2*WNCutoff+2;i<Size;i++) X[i] = 0.0;
				//IFFT;
				IFFT(X);
			}
			//Back to the parent;
			map::SetInitialX();
		}
		virtual void RungeKutta(const double *pX) {
			//Do Runge-Kutta;
			int i;
			func(pX, k1, u1, du1);
			for (i=0;i<Size;i++) k1[i] = CNcoeff[i] * (pX[i] + k1[i]);
			func(k1, k2, u2, du2);
		}
		virtual void BackRungeKutta(const double *pX) {
			//Do backward Runge-Kutta;
			int i;
			func(pX, k1, u1, du1);
			for (i=0;i<Size;i++) {
				k1[i] = -k1[i];
				k1[i] = CNcoeff[i] * (pX[i] + k1[i]);
			}
			func(k1, k2, u2, du2);
			for (i=0;i<Size;i++) k2[i] = -k2[i];
		}
		virtual void func(const double *Xin, double *Xout, double *u, double *du) {
			//Calculate -dt * F_j;
			int i;
			double k;
			//Get u at each spatial point by IFFT;
			u[0] = Xin[0];
			u[1] = 0.0;
			for (i=2;i<2*WNCutoff+2;i++) u[i] = Xin[i];
			for (i=2*WNCutoff+2;i<Size;i++) u[i] = 0.0;
			IFFT(u);
//for (i=0;i<NumPoints;i++) cout << u[i] << " ";
//cout << endl;
			//Get du/dx at each spatial point by IFFT;
			du[0] = 0.0;
			du[1] = 0.0;
			for (i=1;i<=WNCutoff;i++) {
				du[2*i] = - Kcoeff * static_cast<double>(i) * Xin[2*i+1];
				du[2*i+1] = Kcoeff * static_cast<double>(i) * Xin[2*i];
			}
			for (i=2*WNCutoff+2;i<Size;i++) du[i] = 0.0;
			IFFT(du);
//for (i=0;i<NumPoints;i++) cout << du[i] << " ";
//cout << endl;
			//Calculate F at each spatial point and do FFT;
			for (i=0;i<Size;i++) F[i] = u[i] * du[i];
			FFT(F);
			//Calculate -dt * F_j;
			for (i=0;i<Size;i++) Xout[i] = -Param[0] * F[i];
		}
		void dfunc(const double *Xin, const double *dXin, double *dXout, const double *u, const double *du) {
			//Calculate -dt * dF_j;
			//Xin is argument of dF_j, dXin is its multiplier;
			int i;
			double k;
			double *pdF, *pdudc, *pddudc;
			//Set pointer;
			pdF = dF + omp::ThreadID * Size;
			pdudc = dudc + omp::ThreadID * Size;
			pddudc = ddudc + omp::ThreadID * Size;
			//Get du/dc at each spatial point by IFFT;
			pdudc[0] = dXin[0];
			pdudc[1] = 0.0;
			for (i=2;i<2*WNCutoff+2;i++) pdudc[i] = dXin[i];
			for (i=2*WNCutoff+2;i<Size;i++) pdudc[i] = 0.0;
			IFFT(pdudc);
			//Get (d/dc)du/dx at each spatial point by IFFT;
			pddudc[0] = 0.0;
			pddudc[1] = 0.0;
			for (i=1;i<=WNCutoff;i++) {
				pddudc[2*i] = - Kcoeff * static_cast<double>(i) * dXin[2*i+1];
				pddudc[2*i+1] = Kcoeff * static_cast<double>(i) * dXin[2*i];
			}
			for (i=2*WNCutoff+2;i<Size;i++) pddudc[i] = 0.0;
			IFFT(pddudc);
			//Calculate dF at each spatial point and do FFT;
			for (i=0;i<Size;i++) pdF[i] = pdudc[i] * du[i] + u[i] * pddudc[i];
			FFT(pdF);
			//Calculate -dt * dF_j;
			for (i=0;i<Size;i++) dXout[i] = -Param[0] * pdF[i];
		}
		virtual void SingleDF(const double *pX, double *ppdX) {
			//Time evolution of a single vector;
			int i;
			double *pdk1, *pdk2;
			//Set pointer;
			pdk1 = dk1 + omp::ThreadID * Size;
			pdk2 = dk2 + omp::ThreadID * Size;
			//Calculate Jacobian of Runge-Kutta;
			dfunc(pX, ppdX, pdk1, u1, du1);
			for (i=0;i<Size;i++) pdk1[i] = CNcoeff[i] * (ppdX[i] + pdk1[i]);
			dfunc(k1, pdk1, pdk2, u2, du2);
			//Update dx;
			for (i=0;i<Size;i++) ppdX[i] = 0.5 * (CNcoeff[i] * ppdX[i] + pdk1[i] + pdk2[i]);
		}
		virtual void SingleBackDF(const double *pX, double *ppdX) {
			//Backward time evolution of a single vector;
			int i;
			double *pdk1, *pdk2;
			//Set pointer;
			pdk1 = dk1 + omp::ThreadID * Size;
			pdk2 = dk2 + omp::ThreadID * Size;
			//Calculate Jacobian of Runge-Kutta;
			dfunc(pX, ppdX, pdk1, u1, du1);
			for (i=0;i<Size;i++) {
				pdk1[i] = -pdk1[i];
				pdk1[i] = CNcoeff[i] * (ppdX[i] + pdk1[i]);
			}
			dfunc(k1, pdk1, pdk2, u2, du2);
			for (i=0;i<Size;i++) pdk2[i] = -pdk2[i];
			//Update dx;
			for (i=0;i<Size;i++) ppdX[i] = 0.5 * (CNcoeff[i] * ppdX[i] + pdk1[i] + pdk2[i]);
		}
		virtual void df(double *pX, double *pdX) {
			//Time evolution of vectors;
			int i,k;
			//FFT;
			FFT(pX);
			for (i=0;i<NumLyap;i++) FFT(pdX + i * Size);
			//Repeat small steps;
			for (k=0;k<RepNum;k++) {
				//Evolve finite-amplitude perturbation;
				//Runge-Kutta;
				RungeKutta(pX);
				//Compute dx;
#if OMPSW
				#pragma omp parallel private(i)
				#pragma omp for
#endif
				for (i=0;i<NumLyap;i++)
					SingleDF(pX, pdX + i * Size);
				//Update X;
				for (i=0;i<Size;i++) pX[i] = 0.5 * (CNcoeff[i] * pX[i] + k1[i] + k2[i]);
			}
			//Keep Fourier-space values for power spectra;
			if (isRecordingPowSpecDyn || isRecordingOutputQuantity)
				for (i=0;i<Size;i++) FFTdata[i] = pX[i];
			//IFFT;
			IFFT(pX);
			for (i=0;i<NumLyap;i++) IFFT(pdX + i * Size);
		}
		virtual void BackDF(double *pdX) {
			//Backward time evolution of vectors;
			//Values of X are restored here;
			int i,k;
			//Reset X;
			for (i=0;i<Size;i++) X[i] = pStoredX[i];
			//FFT;
			FFT(X);
			for (i=0;i<NumLyap;i++) FFT(pdX + i * Size);
			//Repeat small steps;
			for (k=0;k<RepNum;k++) {
				//Runge-Kutta;
				BackRungeKutta(X);
				//Compute dx;
#if OMPSW
				#pragma omp parallel private(i)
				#pragma omp for
#endif
				for (i=0;i<NumLyap;i++)
					SingleBackDF(X, pdX + i * Size);
				//Update X;
				for (i=0;i<Size;i++) X[i] = 0.5 * (CNcoeff[i] * X[i] + k1[i] + k2[i]);
			}
			//IFFT;
			IFFT(X);
			for (i=0;i<NumLyap;i++) IFFT(pdX + i * Size);
			//Restore X;
			RestoreX();
			//Back to the parent;
			map::BackDF(pdX);
		}
		virtual void f() {
			//Time evolution;
			int i;
			//df does everything;
			df(X, dX);
			//Back to the parent;
			map::f();
		}
		void AccPowSpecDyn() {
			//Calculate power spectrum of dynamics;
			int k;
			double Norm;
			realft2(FFTdata-1, static_cast<unsigned long>(RealSize), -1); 
			PowSpecDyn[0] += FFTdata[0] * FFTdata[0];
			for (k=1;k<RealSize/2;k++)
				PowSpecDyn[k] += FFTdata[k*2] * FFTdata[k*2] + FFTdata[k*2+1] * FFTdata[k*2+1];
			PowSpecDyn[RealSize/2] += FFTdata[1] * FFTdata[1];
		}	
		void CalcOutputQuantity() {
			//Calculate OutputQuantity;
			//	"distance" from UPOs = \sum_k |u(k) - u_{UPO}(k)|^2 = (1/L) \int dx |u(x) - u_{UPO}(x)|^2;
			//
			int i,j,i0,k;
			double x0,x0Tmp;
			double Diff, DblTmp;
			double c, s, c0, s0;
			//
			for (i=0;i<OutputQuantityNum;i++) {
				//estimate shift in space;
				/*
				i0 = (kMaxInd-1)*2;
				x0 = atan2(PowSpecUPO[i*2*WNCutoff+i0]*FFTdata[3+i0] - PowSpecUPO[i*2*WNCutoff+1+i0]*FFTdata[2+i0],
							 PowSpecUPO[i*2*WNCutoff+i0]*FFTdata[2+i0] + PowSpecUPO[i*2*WNCutoff+1+i0]*FFTdata[3+i0]) / (Kcoeff * static_cast<double>(kMaxInd));
				x0Tmp = x0;
				OutputQuantity[i] = DBL_MAX;
				for (k=0;k<kMaxInd;k++) {
					x0 = x0Tmp + static_cast<double>(k) * kMaxLambda;
					Diff = 0.0;
					c = 1.0; s = 0.0;
					c0 = cos(Kcoeff * x0);
					s0 = sin(Kcoeff * x0);
					for (j=1;j<=WNCutoff;j++) {
						DblTmp = c * c0 - s * s0;
						s = s * c0 + c * s0;
						c = DblTmp;
						//
						DblTmp = (c*FFTdata[j*2]+s*FFTdata[j*2+1] - PowSpecUPO[i*2*WNCutoff+(j-1)*2]);
						Diff += DblTmp * DblTmp;
						DblTmp = (-s*FFTdata[j*2]+c*FFTdata[j*2+1] - PowSpecUPO[i*2*WNCutoff+(j-1)*2+1]);
						Diff += DblTmp * DblTmp;
					}
					OutputQuantity[i] = MIN(OutputQuantity[i], Diff);
				}
				*/
				x0Tmp = atan2(PowSpecUPO[i*2*WNCutoff]*FFTdata[3] - PowSpecUPO[i*2*WNCutoff+1]*FFTdata[2],
								 PowSpecUPO[i*2*WNCutoff]*FFTdata[2] + PowSpecUPO[i*2*WNCutoff+1]*FFTdata[3]) / Kcoeff;
				i0 = (kMaxInd-1)*2;
				x0 = atan2(PowSpecUPO[i*2*WNCutoff+i0]*FFTdata[3+i0] - PowSpecUPO[i*2*WNCutoff+1+i0]*FFTdata[2+i0],
							 PowSpecUPO[i*2*WNCutoff+i0]*FFTdata[2+i0] + PowSpecUPO[i*2*WNCutoff+1+i0]*FFTdata[3+i0]) / (Kcoeff * static_cast<double>(kMaxInd));
				x0 += floor((x0Tmp - x0) / kMaxLambda + 0.5) * kMaxLambda;
				c = 1.0; s = 0.0;
				c0 = cos(Kcoeff * x0);
				s0 = sin(Kcoeff * x0);
				Diff = 0.0;
				for (j=1;j<=WNCutoff;j++) {
					DblTmp = c * c0 - s * s0;
					s = s * c0 + c * s0;
					c = DblTmp;
					//
					DblTmp = (c*FFTdata[j*2]+s*FFTdata[j*2+1] - PowSpecUPO[i*2*WNCutoff+(j-1)*2]);
					Diff += DblTmp * DblTmp;
					DblTmp = (-s*FFTdata[j*2]+c*FFTdata[j*2+1] - PowSpecUPO[i*2*WNCutoff+(j-1)*2+1]);
					Diff += DblTmp * DblTmp;
				}
				OutputQuantity[i] = Diff;
			}
			//for (j=1;j<=WNCutoff;j++) F[j-1] = FFTdata[j*2] * FFTdata[j*2] + FFTdata[j*2+1] * FFTdata[j*2+1];
			//for (i=0;i<OutputQuantityNum;i++) {
			//	OutputQuantity[i] = 0.0;
			//	for (j=0;j<WNCutoff;j++) OutputQuantity[i] += fabs(F[j] - PowSpecUPO[i*WNCutoff+j]);
			//}
		}
};


//76: 1d Kuramoto-Sivashinsky for UPO (spectrum method, splitting method (2nd-order Adams-Moulton, Heun));
/*
  du/dt = -d^2u/dx^2 - a d^4u/dx^4 - u du/dx
  u = \sum_j [C_j exp(i k_j x)]
  k_j = 2\pi j/L
  
  Cutoff of wavelength is calculated from system resolution;
*/
class mapKSPDEpUPO : public gcm {
	private:
		bool UseSymm;		//use reflection symmetry to close the trajectory or not;
		bool SymmCt;		//true if next return to the initial condition is with reflection of coordinates;
		int RepNum;			//time step / dt;
		int PeriodTimeStep;	//time step for period;
		int PeriodCt;		//counter for reseting to initial condition;
		int WNCutoff;		//Cutoff of wavenumber (= largest N that satisfies 3N+1<=Size);
		double Pi;			//pi;
		double Kcoeff;		// = 2\pi / L;
		double FFTcoeff, IFFTcoeff;	//Coefficients to multiply results of transform;
		double *CNcoeff;		//Coefficients to multiply for Crank-Nicholson;
		double *k1, *k2, *dk1, *dk2;	//for Runge-Kutta;
		double *F, *dF;			//Nonlinear term in spectral method;
		double *u1, *u2, *du1, *du2;
		double *dudc, *ddudc;	//term evaluation at each spatial point;
		double *XInit;			//initial condition for UPO;
#if FFTWSW
		int Sizec;				// = Size/2+1;
		double *FFTr;			//Array for FFT;
		complex<double> *FFTc;	//Array for FFT;
		fftw_plan FFTplan, IFFTplan;	//plan of FFTW;
		void FFT(double *v) {
			//FFT;
			int i;
			for (i=0;i<Size;i++) FFTr[i] = v[i];
			fftw_execute(FFTplan);
			for (i=0;i<=WNCutoff;i++) {
				v[2*i] = FFTcoeff * FFTc[i].real();
				v[2*i+1] = FFTcoeff * FFTc[i].imag();
			}
			for (i=2*WNCutoff+2;i<Size;i++) v[i] = 0.0;
		}
		void IFFT(double *v) {
			//Inverse FFT;
			int i;
			for (i=0;i<=WNCutoff;i++) FFTc[i] = complex<double>(v[i*2], v[i*2+1]);
			for (i=WNCutoff+1;i<Sizec;i++) FFTc[i] = 0.0;
			fftw_execute(IFFTplan);
			for (i=0;i<Size;i++) v[i] = IFFTcoeff * FFTr[i];
		}
#else
		void FFT(double *v) {
			//FFT;
			int i;
			realft2(v-1, Size, -1);
			v[0] *= FFTcoeff;
			v[1] = 0.0;
			for (i=2;i<2*WNCutoff+2;i++) v[i] *= FFTcoeff;
			for (i=2*WNCutoff+2;i<Size;i++) v[i] = 0.0;
		}
		void IFFT(double *v) {
			//Inverse FFT;
			int i;
			realft2(v-1, Size, 1);
			for (i=0;i<Size;i++) v[i] *= IFFTcoeff;
		}
#endif
	public:
		mapKSPDEpUPO() : gcm(76, "1d Kuramoto-Sivashinsky PDE for UPO (periodic boundary)", 5, 1) {
			//Constructor;
			ParamName[0] = "dt";
			ParamName[1] = "time step / dt";
			ParamName[2] = "Chain length";
			ParamName[3] = "Period of UPO";
			ParamName[4] = "Use reflection symmetry u'=-u, x'=-x";
		}
		~mapKSPDEpUPO() {
			//Destructor;
			//Destroy arrays;
			delete[] k1, k2, dk1, dk2;
			delete[] F, dF, u1, u2, du1, du2, dudc, ddudc, CNcoeff;
			delete[] XInit;
			delete[] OutputQuantity;
#if FFFTWSW
			//For FFTW;
			delete[] FFTr, FFTc;
			fftw_destroy_plan(FFTplan);
			fftw_destroy_plan(IFFTplan);
#endif
		}
		void Init() {
			//Initialize map;
			int i;
			double k;
			ifstream InFile;
			//Go first to the parent;
			gcm::Init();
			//Set variables;
			PeriodTimeStep = floor(Param[3] / (Param[0]*Param[1]) + 0.5);
			Param[0] = Param[3] / static_cast<double>(PeriodTimeStep) / Param[1];
			cout << "dt changed to " << Param[0] << endl;
			cout << "1 period = " << PeriodTimeStep << " time steps." << endl;
			RepNum = static_cast<int>(Param[1]);
			TimeStepLength = Param[0] * Param[1];
			Pi = 3.141592653589793;
			WNCutoff = (Size - 1) / 3;
			Kcoeff = 2.0 * Pi / Param[2];
			UseSymm = static_cast<bool>(Param[4]);
			//Allocate memory for arrays;
			k1 = new double[Size];
			k2 = new double[Size];
			dk1 = new double[Size * omp::NumThreads];
			dk2 = new double[Size * omp::NumThreads];
			F = new double[Size];
			u1 = new double[Size];
			u2 = new double[Size];
			du1 = new double[Size];
			du2 = new double[Size];
			dF = new double[Size * omp::NumThreads];
			dudc = new double[Size * omp::NumThreads];
			ddudc = new double[Size * omp::NumThreads];
			CNcoeff = new double[Size];
			XInit = new double[Size*2];
			//Set arrays;
			for (i=0;i<=WNCutoff;i++) {
				k = Kcoeff * static_cast<double>(i);
				k *= k;
				CNcoeff[2*i] = (2.0 + Param[0] * (k - k * k)) / (2.0 - Param[0] * (k - k * k));
				CNcoeff[2*i+1] = (2.0 + Param[0] * (k - k * k)) / (2.0 - Param[0] * (k - k * k));
			}
			for (i=2*WNCutoff+2;i<Size;i++) CNcoeff[i] = 0.0;
			//For FFT;
#if FFTWSW
			Sizec = Size / 2 + 1;
			FFTr = new double[Size];
			FFTc = new complex<double>[Sizec];
			FFTcoeff = 1.0 / static_cast<double>(Size);
			IFFTcoeff = 1.0;
			FFTplan = fftw_plan_dft_r2c_1d(Size, FFTr, reinterpret_cast<fftw_complex*>(FFTc), FFTW_MEASURE);
			IFFTplan = fftw_plan_dft_c2r_1d(Size, reinterpret_cast<fftw_complex*>(FFTc), FFTr, FFTW_MEASURE);
#else
			FFTcoeff = 1.0 / static_cast<double>(Size);
			IFFTcoeff = 2.0;
#endif
			//Read initial condition file;
			InFile.open((OutFileName+"-init.dat").c_str(), ios::in);
			if (!InFile) ErrorExit("Can't find initial condition file.");
			for (i=0;i<Size;i++) InFile >> XInit[i];
			InFile.close();
			XInit[Size] = -XInit[0];
			for (i=1;i<Size;i++) XInit[Size+i] = -XInit[Size-i];
			//Set OutputQuantity;
			// OutputQuantity[] = inverse of matric C;
			ExistOutputQuantity = false;
			ExistOutputQuantityBack = true;
			OutputQuantityNumBack = SizeMatR;
			OutputQuantity = new double[OutputQuantityNumBack];
			//Revise LoopByte;
			LoopByte += sizeof(int) + sizeof(bool);
		}
		void SetInitialX() {
			//Set initial conditions;
			int i;
			for (i=0;i<Size;i++) X[i] = XInit[i];
			PeriodCt = 0;
			SymmCt = true;
			//Back to the parent;
			map::SetInitialX();
		}
		virtual void RungeKutta(const double *pX) {
			//Do Runge-Kutta;
			int i;
			func(pX, k1, u1, du1);
			for (i=0;i<Size;i++) k1[i] = CNcoeff[i] * (pX[i] + k1[i]);
			func(k1, k2, u2, du2);
		}
		virtual void BackRungeKutta(const double *pX) {
			//Do backward Runge-Kutta;
			int i;
			func(pX, k1, u1, du1);
			for (i=0;i<Size;i++) {
				k1[i] = -k1[i];
				k1[i] = CNcoeff[i] * (pX[i] + k1[i]);
			}
			func(k1, k2, u2, du2);
			for (i=0;i<Size;i++) k2[i] = -k2[i];
		}
		virtual void func(const double *Xin, double *Xout, double *u, double *du) {
			//Calculate -dt * F_j;
			int i;
			double k;
			//Get u at each spatial point by IFFT;
			u[0] = Xin[0];
			u[1] = 0.0;
			for (i=2;i<2*WNCutoff+2;i++) u[i] = Xin[i];
			for (i=2*WNCutoff+2;i<Size;i++) u[i] = 0.0;
			IFFT(u);
//for (i=0;i<NumPoints;i++) cout << u[i] << " ";
//cout << endl;
			//Get du/dx at each spatial point by IFFT;
			du[0] = 0.0;
			du[1] = 0.0;
			for (i=1;i<=WNCutoff;i++) {
				du[2*i] = - Kcoeff * static_cast<double>(i) * Xin[2*i+1];
				du[2*i+1] = Kcoeff * static_cast<double>(i) * Xin[2*i];
			}
			for (i=2*WNCutoff+2;i<Size;i++) du[i] = 0.0;
			IFFT(du);
//for (i=0;i<NumPoints;i++) cout << du[i] << " ";
//cout << endl;
			//Calculate F at each spatial point and do FFT;
			for (i=0;i<Size;i++) F[i] = u[i] * du[i];
			FFT(F);
			//Calculate -dt * F_j;
			for (i=0;i<Size;i++) Xout[i] = -Param[0] * F[i];
		}
		void dfunc(const double *Xin, const double *dXin, double *dXout, const double *u, const double *du) {
			//Calculate -dt * dF_j;
			//Xin is argument of dF_j, dXin is its multiplier;
			int i;
			double k;
			double *pdF, *pdudc, *pddudc;
			//Set pointer;
			pdF = dF + omp::ThreadID * Size;
			pdudc = dudc + omp::ThreadID * Size;
			pddudc = ddudc + omp::ThreadID * Size;
			//Get du/dc at each spatial point by IFFT;
			pdudc[0] = dXin[0];
			pdudc[1] = 0.0;
			for (i=2;i<2*WNCutoff+2;i++) pdudc[i] = dXin[i];
			for (i=2*WNCutoff+2;i<Size;i++) pdudc[i] = 0.0;
			IFFT(pdudc);
			//Get (d/dc)du/dx at each spatial point by IFFT;
			pddudc[0] = 0.0;
			pddudc[1] = 0.0;
			for (i=1;i<=WNCutoff;i++) {
				pddudc[2*i] = - Kcoeff * static_cast<double>(i) * dXin[2*i+1];
				pddudc[2*i+1] = Kcoeff * static_cast<double>(i) * dXin[2*i];
			}
			for (i=2*WNCutoff+2;i<Size;i++) pddudc[i] = 0.0;
			IFFT(pddudc);
			//Calculate dF at each spatial point and do FFT;
			for (i=0;i<Size;i++) pdF[i] = pdudc[i] * du[i] + u[i] * pddudc[i];
			FFT(pdF);
			//Calculate -dt * dF_j;
			for (i=0;i<Size;i++) dXout[i] = -Param[0] * pdF[i];
		}
		virtual void SingleDF(const double *pX, double *ppdX) {
			//Time evolution of a single vector;
			int i;
			double *pdk1, *pdk2;
			//Set pointer;
			pdk1 = dk1 + omp::ThreadID * Size;
			pdk2 = dk2 + omp::ThreadID * Size;
			//Calculate Jacobian of Runge-Kutta;
			dfunc(pX, ppdX, pdk1, u1, du1);
			for (i=0;i<Size;i++) pdk1[i] = CNcoeff[i] * (ppdX[i] + pdk1[i]);
			dfunc(k1, pdk1, pdk2, u2, du2);
			//Update dx;
			for (i=0;i<Size;i++) ppdX[i] = 0.5 * (CNcoeff[i] * ppdX[i] + pdk1[i] + pdk2[i]);
		}
		virtual void SingleBackDF(const double *pX, double *ppdX) {
			//Backward time evolution of a single vector;
			int i;
			double *pdk1, *pdk2;
			//Set pointer;
			pdk1 = dk1 + omp::ThreadID * Size;
			pdk2 = dk2 + omp::ThreadID * Size;
			//Calculate Jacobian of Runge-Kutta;
			dfunc(pX, ppdX, pdk1, u1, du1);
			for (i=0;i<Size;i++) {
				pdk1[i] = -pdk1[i];
				pdk1[i] = CNcoeff[i] * (ppdX[i] + pdk1[i]);
			}
			dfunc(k1, pdk1, pdk2, u2, du2);
			for (i=0;i<Size;i++) pdk2[i] = -pdk2[i];
			//Update dx;
			for (i=0;i<Size;i++) ppdX[i] = 0.5 * (CNcoeff[i] * ppdX[i] + pdk1[i] + pdk2[i]);
		}
		virtual void df(double *pX, double *pdX) {
			//Time evolution of vectors;
			int i,k;
			//FFT;
			FFT(pX);
			for (i=0;i<NumLyap;i++) FFT(pdX + i * Size);
			//Repeat small steps;
			for (k=0;k<RepNum;k++) {
				//Evolve finite-amplitude perturbation;
				//Runge-Kutta;
				RungeKutta(pX);
				//Compute dx;
#if OMPSW
				#pragma omp parallel private(i)
				#pragma omp for
#endif
				for (i=0;i<NumLyap;i++)
					SingleDF(pX, pdX + i * Size);
				//Update X;
				for (i=0;i<Size;i++) pX[i] = 0.5 * (CNcoeff[i] * pX[i] + k1[i] + k2[i]);
			}
			//IFFT;
			IFFT(pX);
			for (i=0;i<NumLyap;i++) IFFT(pdX + i * Size);
		}
		virtual void BackDF(double *pdX) {
			//Backward time evolution of vectors;
			//Values of X are restored here;
			int i,k;
			//Reset X;
			for (i=0;i<Size;i++) X[i] = pStoredX[i];
			//FFT;
			FFT(X);
			for (i=0;i<NumLyap;i++) FFT(pdX + i * Size);
			//Repeat small steps;
			for (k=0;k<RepNum;k++) {
				//Runge-Kutta;
				BackRungeKutta(X);
				//Compute dx;
#if OMPSW
				#pragma omp parallel private(i)
				#pragma omp for
#endif
				for (i=0;i<NumLyap;i++)
					SingleBackDF(X, pdX + i * Size);
				//Update X;
				for (i=0;i<Size;i++) X[i] = 0.5 * (CNcoeff[i] * X[i] + k1[i] + k2[i]);
			}
			//IFFT;
			IFFT(X);
			for (i=0;i<NumLyap;i++) IFFT(pdX + i * Size);
			//Restore X;
			RestoreX();
			//Back to the parent;
			map::BackDF(pdX);
		}
		virtual void f() {
			//Time evolution;
			int i;
			//df does everything;
			df(X, dX);
			//Reset to initial condition;
			PeriodCt++;
			if (PeriodCt == PeriodTimeStep) {
				if (UseSymm && SymmCt)
					for (i=0;i<Size;i++) X[i] = XInit[Size+i];
				else
					for (i=0;i<Size;i++) X[i] = XInit[i];
				PeriodCt = 0;
				SymmCt = !SymmCt;
			}
			//Back to the parent;
			map::f();
		}
		void WriteVar() {
			//Write variables at junctions of forward-backward loops;
			//Go to the parent;
			map::WriteVar();
			//
			OutFileTmp.write(reinterpret_cast<char*>(&PeriodCt),sizeof(int));
			OutFileTmp.write(reinterpret_cast<char*>(&SymmCt),sizeof(bool));
			OutFileTmp.flush();
		}
		void ReadVar() {
			//Read variables at junctions of forward-backward loops;
			//Go to the parent;
			map::ReadVar();
			//Recover the seed;
			OutFileTmp.read(reinterpret_cast<char*>(&PeriodCt),sizeof(int));
			OutFileTmp.read(reinterpret_cast<char*>(&SymmCt),sizeof(bool));
		}
		void WriteVarX() {
			//Write variables at junctions of forward-backward loops;
			//Go to the parent;
			map::WriteVarX();
			//
			OutFileTmp.write(reinterpret_cast<char*>(&PeriodCt),sizeof(int));
			OutFileTmp.write(reinterpret_cast<char*>(&SymmCt),sizeof(bool));
			OutFileTmp.flush();
		}
		void ReadVarX() {
			//Read variables at junctions of forward-backward loops;
			//Go to the parent;
			map::ReadVarX();
			//Recover the seed;
			OutFileTmp.read(reinterpret_cast<char*>(&PeriodCt),sizeof(int));
			OutFileTmp.read(reinterpret_cast<char*>(&SymmCt),sizeof(bool));
		}
		void ReadVarNoDX() {
			//Read variables except dX at junctions of forward-backward loops;
			//Go to the parent;
			map::ReadVarNoDX();
			//
			OutFileTmp.read(reinterpret_cast<char*>(&PeriodCt),sizeof(int));
			OutFileTmp.read(reinterpret_cast<char*>(&SymmCt),sizeof(bool));
		}
		virtual void CalcOutputQuantityBack() {
			//Calculate OutputQuantity;
			int i;
			//Just output matrix C;
			for (i=0;i<SizeMatR;i++) OutputQuantity[i] = MatC[i];
			/*
			//Compute inverse of matrix C (same algorithm as BackwardEvolution());
			//Bug;
			int i,j,k;
			double Sum;
			double *pCinv;	//pointer for OutputQuantity;
			//Solve I = C C';
#if OMPSW
			#pragma omp parallel private(i,j,k,pCinv,Sum)
			#pragma omp for schedule(dynamic)
#endif
			for (k=0;k<NumLyap;k++) {
				//Set pointer;
				pCinv = OutputQuantity + k*(k+1)/2;
				//Main part;
				for (i=k;i>=0;i--) {
					//Sum = 0.0;
					pCinv[i] = (i==k ? 1.0 : 0.0);
					for (j=i+1;j<=k;j++) pCinv[i] -= MatC[i*NumLyap + j - i*(i+1)/2] * pCinv[j];
					pCinv[i] /= MatC[i*NumLyap + i - i*(i+1)/2];
				}
				//Normalization;
				//Sum = 0.0;
				//for (i=0;i<=k;i++) Sum += pCinv[i] * pCinv[i];
				//Sum = 1.0 / sqrt(Sum);
				//for (i=0;i<=k;i++) pMatC[i] *= Sum;
			}
			*/
		}
};

//77: 1d Kuramoto-Sivashinsky with UPO (spectrum method, splitting method (2nd-order Adams-Moulton, Heun));
/*
  du/dt = -d^2u/dx^2 - a d^4u/dx^4 - u du/dx
  u = \sum_j [C_j exp(i k_j x)]
  k_j = 2\pi j/L
  
  Cutoff of wavelength is calculated from system resolution;
*/
class mapKSPDEpUPO2 : public gcm {
	private:
		bool OutputQuantityMode;	//false for correspondence to UPO, true for distance from UPO
									//need time-series data for the latter (1st col = period, 2nd col = sign);
		bool AdjustDt;		//If dt is adjusted to match UPO;
		int RepNum;			//time step / dt;
		int WNCutoff;		//Cutoff of wavenumber (= largest N that satisfies 3N+1<=Size);
		int kMaxInd;		//Index for the most unstable wavenumber;
		int UPOPeriod;		//Period of UPO (unit of time step);
		double Pi;			//pi;
		double Kcoeff;		// = 2\pi / L;
		double FFTcoeff, IFFTcoeff;	//Coefficients to multiply results of transform;
		double kMaxLambda;	//wave length for the most unstable wavenumber;
		double *CNcoeff;		//Coefficients to multiply for Crank-Nicholson;
		double *k1, *k2, *dk1, *dk2;	//for Runge-Kutta;
		double *F, *dF;			//Nonlinear term in spectral method;
		double *u1, *u2, *du1, *du2;
		double *dudc, *ddudc;	//term evaluation at each spatial point;
		double *FFTdata;
		double *PowSpecUPO;		//power spectra of UPOs;
		double *DiffVec;		//Difference vector;
		ifstream InFileUPO;		//Input file for time series of UPO correspondence;
#if FFTWSW
		int Sizec;				// = Size/2+1;
		double *FFTr;			//Array for FFT;
		complex<double> *FFTc;	//Array for FFT;
		fftw_plan FFTplan, IFFTplan;	//plan of FFTW;
		void FFT(double *v) {
			//FFT;
			int i;
			for (i=0;i<Size;i++) FFTr[i] = v[i];
			fftw_execute(FFTplan);
			for (i=0;i<=WNCutoff;i++) {
				v[2*i] = FFTcoeff * FFTc[i].real();
				v[2*i+1] = FFTcoeff * FFTc[i].imag();
			}
			for (i=2*WNCutoff+2;i<Size;i++) v[i] = 0.0;
		}
		void IFFT(double *v) {
			//Inverse FFT;
			int i;
			for (i=0;i<=WNCutoff;i++) FFTc[i] = complex<double>(v[i*2], v[i*2+1]);
			for (i=WNCutoff+1;i<Sizec;i++) FFTc[i] = 0.0;
			fftw_execute(IFFTplan);
			for (i=0;i<Size;i++) v[i] = IFFTcoeff * FFTr[i];
		}
#else
		void FFT(double *v) {
			//FFT;
			int i;
			realft2(v-1, Size, -1);
			v[0] *= FFTcoeff;
			v[1] = 0.0;
			for (i=2;i<2*WNCutoff+2;i++) v[i] *= FFTcoeff;
			for (i=2*WNCutoff+2;i<Size;i++) v[i] = 0.0;
		}
		void IFFT(double *v) {
			//Inverse FFT;
			int i;
			realft2(v-1, Size, 1);
			for (i=0;i<Size;i++) v[i] *= IFFTcoeff;
		}
#endif
		double CalcUPODistance(int Ct, int Sign, double *x0out) {
			//Calculate "distance" from UPOs = \sum_k |u(k) - u_{UPO}(k)|^2 = (1/L) \int dx |u(x) - u_{UPO}(x)|^2;
			int i, i0;
			double Res;
			double x0, x0Tmp, DblTmp;
			double c, s, c0, s0;
			//estimate shift in space;
			x0Tmp = atan2(Sign*PowSpecUPO[Ct*3*WNCutoff]*FFTdata[3] - PowSpecUPO[Ct*3*WNCutoff+1]*FFTdata[2],
							 Sign*PowSpecUPO[Ct*3*WNCutoff]*FFTdata[2] + PowSpecUPO[Ct*3*WNCutoff+1]*FFTdata[3]) / Kcoeff;
			i0 = (kMaxInd-1)*2;
			x0 = atan2(Sign*PowSpecUPO[Ct*3*WNCutoff+i0]*FFTdata[3+i0] - PowSpecUPO[Ct*3*WNCutoff+1+i0]*FFTdata[2+i0],
						 Sign*PowSpecUPO[Ct*3*WNCutoff+i0]*FFTdata[2+i0] + PowSpecUPO[Ct*3*WNCutoff+1+i0]*FFTdata[3+i0]) / (Kcoeff * static_cast<double>(kMaxInd));
			x0 += floor((x0Tmp - x0) / kMaxLambda + 0.5) * kMaxLambda;
			*x0out = x0;
			c = 1.0; s = 0.0;
			c0 = cos(Kcoeff * x0);
			s0 = sin(Kcoeff * x0);
			Res = 0.0;
			for (i=1;i<=WNCutoff;i++) {
				DblTmp = c * c0 - s * s0;
				s = s * c0 + c * s0;
				c = DblTmp;
				//
				DblTmp = (c*FFTdata[i*2]+s*FFTdata[i*2+1] - Sign*PowSpecUPO[Ct*3*WNCutoff+(i-1)*2]);
				Res += DblTmp * DblTmp;
				DblTmp = (-s*FFTdata[i*2]+c*FFTdata[i*2+1] - PowSpecUPO[Ct*3*WNCutoff+(i-1)*2+1]);
				Res += DblTmp * DblTmp;
			}
			Res = sqrt(2.0 * Res);
			return Res;
		}
	public:
		mapKSPDEpUPO2() : gcm(77, "1d Kuramoto-Sivashinsky PDE with UPO (periodic boundary)", 4, 1) {
			//Constructor;
			ParamName[0] = "dt";
			ParamName[1] = "time step / dt";
			ParamName[2] = "Chain length";
			ParamName[3] = "Adjust dt?";
		}
		~mapKSPDEpUPO2() {
			//Destructor;
			//Destroy arrays;
			delete[] k1, k2, dk1, dk2;
			delete[] F, dF, u1, u2, du1, du2, dudc, ddudc, CNcoeff;
			if (isRecordingOutputQuantity) delete[] OutputQuantity;
			if (isRecordingOutputQuantity && OutputQuantityMode) InFileUPO.close();
			delete[] FFTdata;
			delete[] PowSpecUPO;
			delete[] DiffVec;
#if FFFTWSW
			//For FFTW;
			delete[] FFTr, FFTc;
			fftw_destroy_plan(FFTplan);
			fftw_destroy_plan(IFFTplan);
#endif
		}
		void Init() {
			//Initialize map;
			int i,j;
			double k, a;
			ifstream InFileTmp;
			//Go first to the parent;
			gcm::Init();
			//Set variables;
			Pi = 3.141592653589793;
			WNCutoff = (Size - 1) / 3;
			Kcoeff = 2.0 * Pi / Param[2];
			kMaxInd = static_cast<int>(floor(1.0 / (sqrt(2.0) * Kcoeff) + 0.5));
			kMaxLambda = Param[2] / static_cast<double>(kMaxInd);
			AdjustDt = static_cast<bool>(Param[3]);
			//For FFT;
#if FFTWSW
			Sizec = Size / 2 + 1;
			FFTr = new double[Size];
			FFTc = new complex<double>[Sizec];
			FFTcoeff = 1.0 / static_cast<double>(Size);
			IFFTcoeff = 1.0;
			FFTplan = fftw_plan_dft_r2c_1d(Size, FFTr, reinterpret_cast<fftw_complex*>(FFTc), FFTW_MEASURE);
			IFFTplan = fftw_plan_dft_c2r_1d(Size, reinterpret_cast<fftw_complex*>(FFTc), FFTr, FFTW_MEASURE);
#else
			FFTcoeff = 1.0 / static_cast<double>(Size);
			IFFTcoeff = 2.0;
#endif
			//Allocate memory for arrays;
			k1 = new double[Size];
			k2 = new double[Size];
			dk1 = new double[Size * omp::NumThreads];
			dk2 = new double[Size * omp::NumThreads];
			F = new double[Size];
			u1 = new double[Size];
			u2 = new double[Size];
			du1 = new double[Size];
			du2 = new double[Size];
			dF = new double[Size * omp::NumThreads];
			dudc = new double[Size * omp::NumThreads];
			ddudc = new double[Size * omp::NumThreads];
			CNcoeff = new double[Size];
			FFTdata = new double[Size];
			DiffVec = new double[Size];
			//Read UPO data (period + # of time steps + snapshot data);
			InFileTmp.open((OutFileName+"-upo.dat").c_str(), ios::in);
			if (!InFileTmp) ErrorExit("Can't find UPO file.");
			InFileTmp >> a;				//Period;
			InFileTmp >> UPOPeriod;		//Number of time steps per period;
			if (AdjustDt) {
				Param[0] = a / (static_cast<double>(UPOPeriod) * Param[1]);
				cout << "dt changed to " << Param[0] << endl;
			}
			RepNum = static_cast<int>(Param[1]);
			TimeStepLength = Param[0] * Param[1];
			PowSpecUPO = new double[UPOPeriod * WNCutoff*3];
			for (i=0;i<UPOPeriod;i++) {
				InFileTmp >> a;
				for (j=0;j<Size;j++) InFileTmp >> FFTdata[j];
				FFT(FFTdata);
				for (j=1;j<=WNCutoff;j++) {
					PowSpecUPO[i*WNCutoff*3+(j-1)*2] = FFTdata[j*2];
					PowSpecUPO[i*WNCutoff*3+(j-1)*2+1] = FFTdata[j*2+1];
					PowSpecUPO[(i*3+2)*WNCutoff+(j-1)] = FFTdata[j*2] * FFTdata[j*2] + FFTdata[j*2+1] * FFTdata[j*2+1];
				}
			}
			InFileTmp.close();
			//Set arrays;
			for (i=0;i<=WNCutoff;i++) {
				k = Kcoeff * static_cast<double>(i);
				k *= k;
				CNcoeff[2*i] = (2.0 + Param[0] * (k - k * k)) / (2.0 - Param[0] * (k - k * k));
				CNcoeff[2*i+1] = (2.0 + Param[0] * (k - k * k)) / (2.0 - Param[0] * (k - k * k));
			}
			for (i=2*WNCutoff+2;i<Size;i++) CNcoeff[i] = 0.0;
			//Set OutputQuantity;
			if (isRecordingOutputQuantity) {
				InFileUPO.open((OutFileName+"-upots.dat").c_str(), ios::in);
				OutputQuantityMode = static_cast<bool>(InFileUPO);
				cout << "OutputQuantityMode = " << OutputQuantityMode << endl;
			} else {
				OutputQuantityMode = false;
			}
			ExistOutputQuantity = true;
			if (OutputQuantityMode) OutputQuantityNum = 3;
				else OutputQuantityNum = 5;
			ExistOutputQuantityBack = true;
			OutputQuantityNumBack = NumLyap+5;
			OutputQuantity = new double[MAX(OutputQuantityNum, OutputQuantityNumBack)];
		}
		void SetInitialX() {
			//Set initial conditions;
			int i;
			double Amp, Amp1, Amp2, Amp3;
			ifstream InFile;
			//
			InFile.open((OutFileName+"-init.dat").c_str(), ios::in);
			if (InFile) {
				//Read initial condition file;
				cout << "Initial condition read from file." << endl;
				for (i=0;i<Size;i++) InFile >> X[i];
				InFile.close();
			} else {
				//Random initial condition;
				cout << "Random initial condition created." << endl;
				Amp = 0.01 / static_cast<double>(WNCutoff + 1);
				X[0] = 0.0;
				X[1] = 0.0;			//Imaginary part of zero wavelength is zero;
				Amp1 = Amp / (Kcoeff * Kcoeff * sqrt(2.0));
				Amp2 = Amp / (pow(Kcoeff, 4.0) * sqrt(2.0));
				Amp3 = sqrt(Amp / (Kcoeff * static_cast<double>(2*WNCutoff + 1) * sqrt(2.0)));
				for (i=1;i<=WNCutoff;i++) {
					Amp = MIN(MIN(Amp1 / static_cast<double>(i*i), Amp2 / static_cast<double>(i*i*i*i)), Amp3 / sqrt(static_cast<double>(i)));
					X[2*i] = Amp * (2.0 * RNG.genrand_real1() - 1.0);
					X[2*i+1] = Amp * (2.0 * RNG.genrand_real1() - 1.0);
				}
				for (i=2*WNCutoff+2;i<Size;i++) X[i] = 0.0;
				//IFFT;
				IFFT(X);
			}
			//Back to the parent;
			map::SetInitialX();
		}
		virtual void RungeKutta(const double *pX) {
			//Do Runge-Kutta;
			int i;
			func(pX, k1, u1, du1);
			for (i=0;i<Size;i++) k1[i] = CNcoeff[i] * (pX[i] + k1[i]);
			func(k1, k2, u2, du2);
		}
		virtual void BackRungeKutta(const double *pX) {
			//Do backward Runge-Kutta;
			int i;
			func(pX, k1, u1, du1);
			for (i=0;i<Size;i++) {
				k1[i] = -k1[i];
				k1[i] = CNcoeff[i] * (pX[i] + k1[i]);
			}
			func(k1, k2, u2, du2);
			for (i=0;i<Size;i++) k2[i] = -k2[i];
		}
		virtual void func(const double *Xin, double *Xout, double *u, double *du) {
			//Calculate -dt * F_j;
			int i;
			double k;
			//Get u at each spatial point by IFFT;
			u[0] = Xin[0];
			u[1] = 0.0;
			for (i=2;i<2*WNCutoff+2;i++) u[i] = Xin[i];
			for (i=2*WNCutoff+2;i<Size;i++) u[i] = 0.0;
			IFFT(u);
//for (i=0;i<NumPoints;i++) cout << u[i] << " ";
//cout << endl;
			//Get du/dx at each spatial point by IFFT;
			du[0] = 0.0;
			du[1] = 0.0;
			for (i=1;i<=WNCutoff;i++) {
				du[2*i] = - Kcoeff * static_cast<double>(i) * Xin[2*i+1];
				du[2*i+1] = Kcoeff * static_cast<double>(i) * Xin[2*i];
			}
			for (i=2*WNCutoff+2;i<Size;i++) du[i] = 0.0;
			IFFT(du);
//for (i=0;i<NumPoints;i++) cout << du[i] << " ";
//cout << endl;
			//Calculate F at each spatial point and do FFT;
			for (i=0;i<Size;i++) F[i] = u[i] * du[i];
			FFT(F);
			//Calculate -dt * F_j;
			for (i=0;i<Size;i++) Xout[i] = -Param[0] * F[i];
		}
		void dfunc(const double *Xin, const double *dXin, double *dXout, const double *u, const double *du) {
			//Calculate -dt * dF_j;
			//Xin is argument of dF_j, dXin is its multiplier;
			int i;
			double k;
			double *pdF, *pdudc, *pddudc;
			//Set pointer;
			pdF = dF + omp::ThreadID * Size;
			pdudc = dudc + omp::ThreadID * Size;
			pddudc = ddudc + omp::ThreadID * Size;
			//Get du/dc at each spatial point by IFFT;
			pdudc[0] = dXin[0];
			pdudc[1] = 0.0;
			for (i=2;i<2*WNCutoff+2;i++) pdudc[i] = dXin[i];
			for (i=2*WNCutoff+2;i<Size;i++) pdudc[i] = 0.0;
			IFFT(pdudc);
			//Get (d/dc)du/dx at each spatial point by IFFT;
			pddudc[0] = 0.0;
			pddudc[1] = 0.0;
			for (i=1;i<=WNCutoff;i++) {
				pddudc[2*i] = - Kcoeff * static_cast<double>(i) * dXin[2*i+1];
				pddudc[2*i+1] = Kcoeff * static_cast<double>(i) * dXin[2*i];
			}
			for (i=2*WNCutoff+2;i<Size;i++) pddudc[i] = 0.0;
			IFFT(pddudc);
			//Calculate dF at each spatial point and do FFT;
			for (i=0;i<Size;i++) pdF[i] = pdudc[i] * du[i] + u[i] * pddudc[i];
			FFT(pdF);
			//Calculate -dt * dF_j;
			for (i=0;i<Size;i++) dXout[i] = -Param[0] * pdF[i];
		}
		virtual void SingleDF(const double *pX, double *ppdX) {
			//Time evolution of a single vector;
			int i;
			double *pdk1, *pdk2;
			//Set pointer;
			pdk1 = dk1 + omp::ThreadID * Size;
			pdk2 = dk2 + omp::ThreadID * Size;
			//Calculate Jacobian of Runge-Kutta;
			dfunc(pX, ppdX, pdk1, u1, du1);
			for (i=0;i<Size;i++) pdk1[i] = CNcoeff[i] * (ppdX[i] + pdk1[i]);
			dfunc(k1, pdk1, pdk2, u2, du2);
			//Update dx;
			for (i=0;i<Size;i++) ppdX[i] = 0.5 * (CNcoeff[i] * ppdX[i] + pdk1[i] + pdk2[i]);
		}
		virtual void SingleBackDF(const double *pX, double *ppdX) {
			//Backward time evolution of a single vector;
			int i;
			double *pdk1, *pdk2;
			//Set pointer;
			pdk1 = dk1 + omp::ThreadID * Size;
			pdk2 = dk2 + omp::ThreadID * Size;
			//Calculate Jacobian of Runge-Kutta;
			dfunc(pX, ppdX, pdk1, u1, du1);
			for (i=0;i<Size;i++) {
				pdk1[i] = -pdk1[i];
				pdk1[i] = CNcoeff[i] * (ppdX[i] + pdk1[i]);
			}
			dfunc(k1, pdk1, pdk2, u2, du2);
			for (i=0;i<Size;i++) pdk2[i] = -pdk2[i];
			//Update dx;
			for (i=0;i<Size;i++) ppdX[i] = 0.5 * (CNcoeff[i] * ppdX[i] + pdk1[i] + pdk2[i]);
		}
		virtual void df(double *pX, double *pdX) {
			//Time evolution of vectors;
			int i,k;
			//FFT;
			FFT(pX);
			for (i=0;i<NumLyap;i++) FFT(pdX + i * Size);
			//Repeat small steps;
			for (k=0;k<RepNum;k++) {
				//Evolve finite-amplitude perturbation;
				//Runge-Kutta;
				RungeKutta(pX);
				//Compute dx;
#if OMPSW
				#pragma omp parallel private(i)
				#pragma omp for
#endif
				for (i=0;i<NumLyap;i++)
					SingleDF(pX, pdX + i * Size);
				//Update X;
				for (i=0;i<Size;i++) pX[i] = 0.5 * (CNcoeff[i] * pX[i] + k1[i] + k2[i]);
			}
			//Keep Fourier-space values;
			if (isRecordingPowSpecDyn || isRecordingOutputQuantity)
				for (i=0;i<Size;i++) FFTdata[i] = pX[i];
			//IFFT;
			IFFT(pX);
			for (i=0;i<NumLyap;i++) IFFT(pdX + i * Size);
		}
		virtual void BackDF(double *pdX) {
			//Backward time evolution of vectors;
			//Values of X are restored here;
			int i,k;
			//Reset X;
			for (i=0;i<Size;i++) X[i] = pStoredX[i];
			//FFT;
			FFT(X);
			for (i=0;i<NumLyap;i++) FFT(pdX + i * Size);
			//Repeat small steps;
			for (k=0;k<RepNum;k++) {
				//Runge-Kutta;
				BackRungeKutta(X);
				//Compute dx;
#if OMPSW
				#pragma omp parallel private(i)
				#pragma omp for
#endif
				for (i=0;i<NumLyap;i++)
					SingleBackDF(X, pdX + i * Size);
				//Update X;
				for (i=0;i<Size;i++) X[i] = 0.5 * (CNcoeff[i] * X[i] + k1[i] + k2[i]);
			}
			//IFFT;
			IFFT(X);
			for (i=0;i<NumLyap;i++) IFFT(pdX + i * Size);
			//Restore X;
			RestoreX();
			//Back to the parent;
			map::BackDF(pdX);
		}
		virtual void f() {
			//Time evolution;
			int i;
			//df does everything;
			df(X, dX);
			//Back to the parent;
			map::f();
		}
		void AccPowSpecDyn() {
			//Calculate power spectrum of dynamics;
			int k;
			double Norm;
			PowSpecDyn[0] += FFTdata[0] * FFTdata[0];
			for (k=1;k<RealSize/2;k++)
				PowSpecDyn[k] += FFTdata[k*2] * FFTdata[k*2] + FFTdata[k*2+1] * FFTdata[k*2+1];
			PowSpecDyn[RealSize/2] += FFTdata[1] * FFTdata[1];
		}
		void CalcOutputQuantity() {
			//Calculate OutputQuantity;
			//
			int i,j,i0,k;
			int UPOTime;		//Time counter for UPO (-1 when trajectory is far from UPO);
			int UPOSign;		//sign of UPO;
			double DblTmp, DblTmp2, MinDev, MinDevX0;
			double x0, x0p;
			//
			if (OutputQuantityMode) {
				InFileUPO >> UPOTime;
				InFileUPO >> UPOSign;
				OutputQuantity[0] = CalcUPODistance(0, +1, &x0);
				OutputQuantity[1] = (UPOTime < 0 ? -1.0 : CalcUPODistance(UPOTime, UPOSign, &x0));
				OutputQuantity[2] = (UPOTime < 0 ? -1.0 : x0);
			} else {
				//look for the minimum deviation in power spectrum;
				MinDev = DBL_MAX;
				for (i=0;i<UPOPeriod;i++) {
					/*
					DblTmp = 0.0;
					for (j=1;j<=WNCutoff;j++)
						DblTmp += fabs(FFTdata[j*2] * FFTdata[j*2] + FFTdata[j*2+1] * FFTdata[j*2+1] - PowSpecUPO[(i*3+2)*WNCutoff+(j-1)]);
					if (DblTmp < MinDev) {
						MinDev = DblTmp;
						UPOTime = i;
					}
					*/
					DblTmp = CalcUPODistance(i, +1, &x0);
					DblTmp2 = CalcUPODistance(i, -1, &x0p);
					if (MIN(DblTmp, DblTmp2) < MinDev) {
						if (DblTmp < DblTmp2) {
							MinDev = DblTmp; MinDevX0 = x0;
						} else {
							MinDev = DblTmp2; MinDevX0 = x0p;
						}
						UPOSign = (DblTmp < DblTmp2 ? +1 : -1);
						UPOTime = i;
					}
				}
				//compute distance;
				/*
				DblTmp = CalcUPODistance(UPOTime, +1);
				DblTmp2 = CalcUPODistance(UPOTime, -1);
				//MinDev = MIN(DblTmp, DblTmp2);
				//Set sign and time counter of UPO;
				UPOSign = (DblTmp < DblTmp2 ? +1 : -1);
				*/
				OutputQuantity[0] = CalcUPODistance(0, +1, &x0);
				OutputQuantity[1] = MinDev;
				OutputQuantity[2] = static_cast<double>(UPOTime);
				OutputQuantity[3] = static_cast<double>(UPOSign);
				OutputQuantity[4] = MinDevX0;
			}
		}
#if ARMASW
		mat SubspaceMatrix(int StartInd, int NumInd) {
			//Construct (Size)d-orthogonal-bases matrix for subspace spanned by CV #(StartInd)-#(StartInd + NumInd - 1);
			//See TeX note for subspace-angle;
			int i,j,k,n;
			int Ind;
			double Norm;
			double *pdX;
			mat MatTmp, MatR, MatQ;
			MatTmp.zeros(Size, NumInd);
			//Calculate orthogonal-bases matrix;
			for (j=0;j<NumInd;j++) {
				Ind = StartInd + j;
				//Orthogonalize;
				pdX = CV + Ind * Size;
				for (i=0;i<Size;i++) MatTmp.at(i,j) = pdX[i];
			}
			qr(MatQ, MatR, MatTmp);
			return MatQ.cols(0,NumInd-1);
		}
#endif
		void CalcOutputQuantityBack() {
			//Calculate OutputQuantity;
			//
			// OutputQuantity[0] = || u - u_p || (distance from UPO);
			// OutputQuantity[1-NumLyap] = angle between (u-u_p) and subspace 1-i ;
			// OutputQuantity[(NumLyap+1)-(NumLyap+4) = angle between (u-u_p) and subspace 2-4, 5-9, 2-9, 10-NumLyap;
			//
			int i,j,i0,k;
			int UPOTime;		//Time counter for UPO (-1 when trajectory is far from UPO);
			int UPOSign;		//sign of UPO;
			double DblTmp, DblTmp2, MinDev, MinDevX0;
			double x0, x0p;
			double c,s,c0,s0;
#if ARMASW
			mat MatDiffVec;
			colvec SingVal;
#endif
			//
			//Reset X;
			for (i=0;i<Size;i++) {X[i] = pStoredX[i]; FFTdata[i] = pStoredX[i];}
			FFT(FFTdata);
			//look for the minimum deviation in power spectrum;
			MinDev = DBL_MAX;
			for (i=0;i<UPOPeriod;i++) {
				DblTmp = CalcUPODistance(i, +1, &x0);
				DblTmp2 = CalcUPODistance(i, -1, &x0p);
				if (MIN(DblTmp, DblTmp2) < MinDev) {
					if (DblTmp < DblTmp2) {
						MinDev = DblTmp; MinDevX0 = x0;
					} else {
						MinDev = DblTmp2; MinDevX0 = x0p;
					}
					UPOSign = (DblTmp < DblTmp2 ? +1 : -1);
					UPOTime = i;
				}
			}
			//compute difference vector;
			for (i=0;i<Size;i++) DiffVec[i] = 0.0;
			c = 1.0; s = 0.0;
			c0 = cos(-Kcoeff * MinDevX0);
			s0 = sin(-Kcoeff * MinDevX0);
			for (i=1;i<=WNCutoff;i++) {
				DblTmp = c * c0 - s * s0;
				s = s * c0 + c * s0;
				c = DblTmp;
				//
				DiffVec[i*2] = c * UPOSign*PowSpecUPO[UPOTime*3*WNCutoff+(i-1)*2] + s * PowSpecUPO[UPOTime*3*WNCutoff+(i-1)*2+1];
				DiffVec[i*2+1] = -s * UPOSign*PowSpecUPO[UPOTime*3*WNCutoff+(i-1)*2] + c * PowSpecUPO[UPOTime*3*WNCutoff+(i-1)*2+1];
			}
			IFFT(DiffVec);
			for (i=0;i<Size;i++) DiffVec[i] -= X[i];
			DblTmp = 0.0;
			for (i=0;i<Size;i++) DblTmp += DiffVec[i]*DiffVec[i];
			OutputQuantity[0] = sqrt(DblTmp);
			DblTmp = 1.0 / OutputQuantity[0];
			for (i=0;i<Size;i++) DiffVec[i] *= DblTmp;
			//compute angle;
#if ARMASW
			MatDiffVec.zeros(Size, 1);
			for (i=0;i<Size;i++) MatDiffVec.at(i,0) = DiffVec[i];
			for (i=0;i<NumLyap;i++) {
				SingVal = svd(trans(MatDiffVec) * SubspaceMatrix(0, i+1));
				OutputQuantity[i+1] = acos(MIN(SingVal[0], 1.0));
			}
			SingVal = svd(trans(MatDiffVec) * SubspaceMatrix(1, 3));
			OutputQuantity[NumLyap+1] = acos(MIN(SingVal[0], 1.0));
			SingVal = svd(trans(MatDiffVec) * SubspaceMatrix(4, 5));
			OutputQuantity[NumLyap+2] = acos(MIN(SingVal[0], 1.0));
			SingVal = svd(trans(MatDiffVec) * SubspaceMatrix(1, 8));
			OutputQuantity[NumLyap+3] = acos(MIN(SingVal[0], 1.0));
			SingVal = svd(trans(MatDiffVec) * SubspaceMatrix(9, NumLyap-9));
			OutputQuantity[NumLyap+4] = acos(MIN(SingVal[0], 1.0));
#else
			for (i=0;i<NumLyap;i++) {
				DblTmp = 0.0;
				for (j=0;j<Size;j++) DblTmp += DiffVec[j] * CV[i*Size+j];
				//OutputQuantity[i+1] = acos(MIN(DblTmp, 1.0));
				OutputQuantity[i+1] = 0.0;
			}
			for (i=0;i<4;i++) OutputQuantity[NumLyap+1+i] = 0.0;
#endif
		}
};

#
//84: 1d Kuramoto-Sivashinsky for RPO (spectrum method, splitting method (2nd-order Adams-Moulton, Heun));
/*
  du/dt = -d^2u/dx^2 - a d^4u/dx^4 - u du/dx
  u = \sum_j [C_j exp(i k_j x)]
  k_j = 2\pi j/L
  
  Cutoff of wavelength is calculated from system resolution;
*/
class mapKSPDEpRPO : public gcm {
	private:
		int RepNum;			//time step / dt;
		int PeriodTimeStep;	//time step for period;
		int PeriodCt;		//counter for reseting to initial condition;
		int WNCutoff;		//Cutoff of wavenumber (= largest N that satisfies 3N+1<=Size);
		double Pi;			//pi;
		double Kcoeff;		// = 2\pi / L;
		double FFTcoeff, IFFTcoeff;	//Coefficients to multiply results of transform;
		double *CNcoeff;		//Coefficients to multiply for Crank-Nicholson;
		double *k1, *k2, *dk1, *dk2;	//for Runge-Kutta;
		double *F, *dF;			//Nonlinear term in spectral method;
		double *u1, *u2, *du1, *du2;
		double *dudc, *ddudc;	//term evaluation at each spatial point;
		double *XInit;			//initial condition for UPO;
		double *ShiftCoeff;		//Fourier multipliers for shift;
#if FFTWSW
		int Sizec;				// = Size/2+1;
		double *FFTr;			//Array for FFT;
		complex<double> *FFTc;	//Array for FFT;
		fftw_plan FFTplan, IFFTplan;	//plan of FFTW;
		void FFT(double *v) {
			//FFT;
			int i;
			for (i=0;i<Size;i++) FFTr[i] = v[i];
			fftw_execute(FFTplan);
			for (i=0;i<=WNCutoff;i++) {
				v[2*i] = FFTcoeff * FFTc[i].real();
				v[2*i+1] = FFTcoeff * FFTc[i].imag();
			}
			for (i=2*WNCutoff+2;i<Size;i++) v[i] = 0.0;
		}
		void IFFT(double *v) {
			//Inverse FFT;
			int i;
			for (i=0;i<=WNCutoff;i++) FFTc[i] = complex<double>(v[i*2], v[i*2+1]);
			for (i=WNCutoff+1;i<Sizec;i++) FFTc[i] = 0.0;
			fftw_execute(IFFTplan);
			for (i=0;i<Size;i++) v[i] = IFFTcoeff * FFTr[i];
		}
#else
		void FFT(double *v) {
			//FFT;
			int i;
			realft2(v-1, Size, -1);
			v[0] *= FFTcoeff;
			v[1] = 0.0;
			for (i=2;i<2*WNCutoff+2;i++) v[i] *= FFTcoeff;
			for (i=2*WNCutoff+2;i<Size;i++) v[i] = 0.0;
		}
		void IFFT(double *v) {
			//Inverse FFT;
			int i;
			realft2(v-1, Size, 1);
			for (i=0;i<Size;i++) v[i] *= IFFTcoeff;
		}
#endif
	public:
		mapKSPDEpRPO() : gcm(84, "1d Kuramoto-Sivashinsky PDE for RPO (periodic boundary)", 5, 1) {
			//Constructor;
			ParamName[0] = "dt";
			ParamName[1] = "time step / dt";
			ParamName[2] = "Chain length";
			ParamName[3] = "Period of UPO";
			ParamName[4] = "Spatial shift";
		}
		~mapKSPDEpRPO() {
			//Destructor;
			//Destroy arrays;
			delete[] k1, k2, dk1, dk2;
			delete[] F, dF, u1, u2, du1, du2, dudc, ddudc, CNcoeff;
			delete[] XInit, ShiftCoeff;
			delete[] OutputQuantity;
#if FFFTWSW
			//For FFTW;
			delete[] FFTr, FFTc;
			fftw_destroy_plan(FFTplan);
			fftw_destroy_plan(IFFTplan);
#endif
		}
		void Init() {
			//Initialize map;
			int i;
			double k;
			ifstream InFile;
			//Go first to the parent;
			gcm::Init();
			//Set variables;
			PeriodTimeStep = floor(Param[3] / (Param[0]*Param[1]) + 0.5);
			Param[0] = Param[3] / static_cast<double>(PeriodTimeStep) / Param[1];
			cout << "dt changed to " << Param[0] << endl;
			cout << "1 period = " << PeriodTimeStep << " time steps." << endl;
			RepNum = static_cast<int>(Param[1]);
			TimeStepLength = Param[0] * Param[1];
			Pi = 3.141592653589793;
			WNCutoff = (Size - 1) / 3;
			Kcoeff = 2.0 * Pi / Param[2];
			//Allocate memory for arrays;
			k1 = new double[Size];
			k2 = new double[Size];
			dk1 = new double[Size * omp::NumThreads];
			dk2 = new double[Size * omp::NumThreads];
			F = new double[Size];
			u1 = new double[Size];
			u2 = new double[Size];
			du1 = new double[Size];
			du2 = new double[Size];
			dF = new double[Size * omp::NumThreads];
			dudc = new double[Size * omp::NumThreads];
			ddudc = new double[Size * omp::NumThreads];
			CNcoeff = new double[Size];
			XInit = new double[Size];
			ShiftCoeff = new double[Size];
			//Set arrays;
			for (i=0;i<=WNCutoff;i++) {
				k = Kcoeff * static_cast<double>(i);
				k *= k;
				CNcoeff[2*i] = (2.0 + Param[0] * (k - k * k)) / (2.0 - Param[0] * (k - k * k));
				CNcoeff[2*i+1] = (2.0 + Param[0] * (k - k * k)) / (2.0 - Param[0] * (k - k * k));
			}
			for (i=2*WNCutoff+2;i<Size;i++) CNcoeff[i] = 0.0;
			//For FFT;
#if FFTWSW
			Sizec = Size / 2 + 1;
			FFTr = new double[Size];
			FFTc = new complex<double>[Sizec];
			FFTcoeff = 1.0 / static_cast<double>(Size);
			IFFTcoeff = 1.0;
			FFTplan = fftw_plan_dft_r2c_1d(Size, FFTr, reinterpret_cast<fftw_complex*>(FFTc), FFTW_MEASURE);
			IFFTplan = fftw_plan_dft_c2r_1d(Size, reinterpret_cast<fftw_complex*>(FFTc), FFTr, FFTW_MEASURE);
#else
			FFTcoeff = 1.0 / static_cast<double>(Size);
			IFFTcoeff = 2.0;
#endif
			//Read initial condition file;
			InFile.open((OutFileName+"-init.dat").c_str(), ios::in);
			if (!InFile) ErrorExit("Can't find initial condition file.");
			for (i=0;i<Size;i++) InFile >> XInit[i];
			InFile.close();
			for (i=0;i<=WNCutoff;i++) {
				ShiftCoeff[2*i] = cos(Kcoeff * static_cast<double>(i) * Param[4]);
				ShiftCoeff[2*i+1] = sin(Kcoeff * static_cast<double>(i) * Param[4]);
			}
			for (i=2*WNCutoff+2;i<Size;i++) ShiftCoeff[i] = 0.0;			
			//Set OutputQuantity;
			// OutputQuantity[] = inverse of matric C;
			ExistOutputQuantity = false;
			ExistOutputQuantityBack = true;
			OutputQuantityNumBack = SizeMatR;
			OutputQuantity = new double[OutputQuantityNumBack];
			//Revise LoopByte;
			LoopByte += sizeof(int);
		}
		void SetInitialX() {
			//Set initial conditions;
			int i;
			for (i=0;i<Size;i++) X[i] = XInit[i];
			IFFT(X);
			PeriodCt = 0;
			//Back to the parent;
			map::SetInitialX();
		}
		virtual void RungeKutta(const double *pX) {
			//Do Runge-Kutta;
			int i;
			func(pX, k1, u1, du1);
			for (i=0;i<Size;i++) k1[i] = CNcoeff[i] * (pX[i] + k1[i]);
			func(k1, k2, u2, du2);
		}
		virtual void BackRungeKutta(const double *pX) {
			//Do backward Runge-Kutta;
			int i;
			func(pX, k1, u1, du1);
			for (i=0;i<Size;i++) {
				k1[i] = -k1[i];
				k1[i] = CNcoeff[i] * (pX[i] + k1[i]);
			}
			func(k1, k2, u2, du2);
			for (i=0;i<Size;i++) k2[i] = -k2[i];
		}
		virtual void func(const double *Xin, double *Xout, double *u, double *du) {
			//Calculate -dt * F_j;
			int i;
			double k;
			//Get u at each spatial point by IFFT;
			u[0] = Xin[0];
			u[1] = 0.0;
			for (i=2;i<2*WNCutoff+2;i++) u[i] = Xin[i];
			for (i=2*WNCutoff+2;i<Size;i++) u[i] = 0.0;
			IFFT(u);
//for (i=0;i<NumPoints;i++) cout << u[i] << " ";
//cout << endl;
			//Get du/dx at each spatial point by IFFT;
			du[0] = 0.0;
			du[1] = 0.0;
			for (i=1;i<=WNCutoff;i++) {
				du[2*i] = - Kcoeff * static_cast<double>(i) * Xin[2*i+1];
				du[2*i+1] = Kcoeff * static_cast<double>(i) * Xin[2*i];
			}
			for (i=2*WNCutoff+2;i<Size;i++) du[i] = 0.0;
			IFFT(du);
//for (i=0;i<NumPoints;i++) cout << du[i] << " ";
//cout << endl;
			//Calculate F at each spatial point and do FFT;
			for (i=0;i<Size;i++) F[i] = u[i] * du[i];
			FFT(F);
			//Calculate -dt * F_j;
			for (i=0;i<Size;i++) Xout[i] = -Param[0] * F[i];
		}
		void dfunc(const double *Xin, const double *dXin, double *dXout, const double *u, const double *du) {
			//Calculate -dt * dF_j;
			//Xin is argument of dF_j, dXin is its multiplier;
			int i;
			double k;
			double *pdF, *pdudc, *pddudc;
			//Set pointer;
			pdF = dF + omp::ThreadID * Size;
			pdudc = dudc + omp::ThreadID * Size;
			pddudc = ddudc + omp::ThreadID * Size;
			//Get du/dc at each spatial point by IFFT;
			pdudc[0] = dXin[0];
			pdudc[1] = 0.0;
			for (i=2;i<2*WNCutoff+2;i++) pdudc[i] = dXin[i];
			for (i=2*WNCutoff+2;i<Size;i++) pdudc[i] = 0.0;
			IFFT(pdudc);
			//Get (d/dc)du/dx at each spatial point by IFFT;
			pddudc[0] = 0.0;
			pddudc[1] = 0.0;
			for (i=1;i<=WNCutoff;i++) {
				pddudc[2*i] = - Kcoeff * static_cast<double>(i) * dXin[2*i+1];
				pddudc[2*i+1] = Kcoeff * static_cast<double>(i) * dXin[2*i];
			}
			for (i=2*WNCutoff+2;i<Size;i++) pddudc[i] = 0.0;
			IFFT(pddudc);
			//Calculate dF at each spatial point and do FFT;
			for (i=0;i<Size;i++) pdF[i] = pdudc[i] * du[i] + u[i] * pddudc[i];
			FFT(pdF);
			//Calculate -dt * dF_j;
			for (i=0;i<Size;i++) dXout[i] = -Param[0] * pdF[i];
		}
		virtual void SingleDF(const double *pX, double *ppdX) {
			//Time evolution of a single vector;
			int i;
			double *pdk1, *pdk2;
			//Set pointer;
			pdk1 = dk1 + omp::ThreadID * Size;
			pdk2 = dk2 + omp::ThreadID * Size;
			//Calculate Jacobian of Runge-Kutta;
			dfunc(pX, ppdX, pdk1, u1, du1);
			for (i=0;i<Size;i++) pdk1[i] = CNcoeff[i] * (ppdX[i] + pdk1[i]);
			dfunc(k1, pdk1, pdk2, u2, du2);
			//Update dx;
			for (i=0;i<Size;i++) ppdX[i] = 0.5 * (CNcoeff[i] * ppdX[i] + pdk1[i] + pdk2[i]);
		}
		virtual void SingleBackDF(const double *pX, double *ppdX) {
			//Backward time evolution of a single vector;
			int i;
			double *pdk1, *pdk2;
			//Set pointer;
			pdk1 = dk1 + omp::ThreadID * Size;
			pdk2 = dk2 + omp::ThreadID * Size;
			//Calculate Jacobian of Runge-Kutta;
			dfunc(pX, ppdX, pdk1, u1, du1);
			for (i=0;i<Size;i++) {
				pdk1[i] = -pdk1[i];
				pdk1[i] = CNcoeff[i] * (ppdX[i] + pdk1[i]);
			}
			dfunc(k1, pdk1, pdk2, u2, du2);
			for (i=0;i<Size;i++) pdk2[i] = -pdk2[i];
			//Update dx;
			for (i=0;i<Size;i++) ppdX[i] = 0.5 * (CNcoeff[i] * ppdX[i] + pdk1[i] + pdk2[i]);
		}
		virtual void df(double *pX, double *pdX) {
			//Time evolution of vectors;
			int i,k;
			//FFT;
			FFT(pX);
			for (i=0;i<NumLyap;i++) FFT(pdX + i * Size);
			//Repeat small steps;
			for (k=0;k<RepNum;k++) {
				//Evolve finite-amplitude perturbation;
				//Runge-Kutta;
				RungeKutta(pX);
				//Compute dx;
#if OMPSW
				#pragma omp parallel private(i)
				#pragma omp for
#endif
				for (i=0;i<NumLyap;i++)
					SingleDF(pX, pdX + i * Size);
				//Update X;
				for (i=0;i<Size;i++) pX[i] = 0.5 * (CNcoeff[i] * pX[i] + k1[i] + k2[i]);
			}
			//IFFT;
			IFFT(pX);
			for (i=0;i<NumLyap;i++) IFFT(pdX + i * Size);
		}
		virtual void BackDF(double *pdX) {
			//Backward time evolution of vectors;
			//Values of X are restored here;
			int i,k;
			//Reset X;
			for (i=0;i<Size;i++) X[i] = pStoredX[i];
			//FFT;
			FFT(X);
			for (i=0;i<NumLyap;i++) FFT(pdX + i * Size);
			//Repeat small steps;
			for (k=0;k<RepNum;k++) {
				//Runge-Kutta;
				BackRungeKutta(X);
				//Compute dx;
#if OMPSW
				#pragma omp parallel private(i)
				#pragma omp for
#endif
				for (i=0;i<NumLyap;i++)
					SingleBackDF(X, pdX + i * Size);
				//Update X;
				for (i=0;i<Size;i++) X[i] = 0.5 * (CNcoeff[i] * X[i] + k1[i] + k2[i]);
			}
			//IFFT;
			IFFT(X);
			for (i=0;i<NumLyap;i++) IFFT(pdX + i * Size);
			//Restore X;
			RestoreX();
			//Back to the parent;
			map::BackDF(pdX);
		}
		virtual void f() {
			//Time evolution;
			int i,j;
			double *pdX;
			double TmpDbl;
			//df does everything;
			df(X, dX);
			//Reset to initial condition & shift vectors by -lambda;
			PeriodCt++;
			if (PeriodCt == PeriodTimeStep) {
				for (i=0;i<Size;i++) X[i] = XInit[i];
				for (i=0;i<NumLyap;i++) {
					pdX = dX + i * Size;
					FFT(pdX);
					for (i=0;i<=WNCutoff;i++) {
						TmpDbl = pdX[2*i];
						pdX[2*i] = ShiftCoeff[2*i] * pdX[2*i] - ShiftCoeff[2*i+1] * pdX[2*i+1];
						pdX[2*i+1] = ShiftCoeff[2*i] * pdX[2*i+1] + ShiftCoeff[2*i+1] * TmpDbl;
					}
					IFFT(pdX);
				}
				PeriodCt = 0;
			}
			//Back to the parent;
			map::f();
		}
		void WriteVar() {
			//Write variables at junctions of forward-backward loops;
			//Go to the parent;
			map::WriteVar();
			//
			OutFileTmp.write(reinterpret_cast<char*>(&PeriodCt),sizeof(int));
			OutFileTmp.flush();
		}
		void ReadVar() {
			//Read variables at junctions of forward-backward loops;
			//Go to the parent;
			map::ReadVar();
			//Recover the seed;
			OutFileTmp.read(reinterpret_cast<char*>(&PeriodCt),sizeof(int));
		}
		void WriteVarX() {
			//Write variables at junctions of forward-backward loops;
			//Go to the parent;
			map::WriteVarX();
			//
			OutFileTmp.write(reinterpret_cast<char*>(&PeriodCt),sizeof(int));
			OutFileTmp.flush();
		}
		void ReadVarX() {
			//Read variables at junctions of forward-backward loops;
			//Go to the parent;
			map::ReadVarX();
			//Recover the seed;
			OutFileTmp.read(reinterpret_cast<char*>(&PeriodCt),sizeof(int));
		}
		void ReadVarNoDX() {
			//Read variables except dX at junctions of forward-backward loops;
			//Go to the parent;
			map::ReadVarNoDX();
			//
			OutFileTmp.read(reinterpret_cast<char*>(&PeriodCt),sizeof(int));
		}
		virtual void CalcOutputQuantityBack() {
			//Calculate OutputQuantity;
			int i;
			//Just output matrix C;
			for (i=0;i<SizeMatR;i++) OutputQuantity[i] = MatC[i];
			/*
			//Compute inverse of matrix C (same algorithm as BackwardEvolution());
			//Bug;
			int i,j,k;
			double Sum;
			double *pCinv;	//pointer for OutputQuantity;
			//Solve I = C C';
#if OMPSW
			#pragma omp parallel private(i,j,k,pCinv,Sum)
			#pragma omp for schedule(dynamic)
#endif
			for (k=0;k<NumLyap;k++) {
				//Set pointer;
				pCinv = OutputQuantity + k*(k+1)/2;
				//Main part;
				for (i=k;i>=0;i--) {
					//Sum = 0.0;
					pCinv[i] = (i==k ? 1.0 : 0.0);
					for (j=i+1;j<=k;j++) pCinv[i] -= MatC[i*NumLyap + j - i*(i+1)/2] * pCinv[j];
					pCinv[i] /= MatC[i*NumLyap + i - i*(i+1)/2];
				}
				//Normalization;
				//Sum = 0.0;
				//for (i=0;i<=k;i++) Sum += pCinv[i] * pCinv[i];
				//Sum = 1.0 / sqrt(Sum);
				//for (i=0;i<=k;i++) pMatC[i] *= Sum;
			}
			*/
		}
};


//Class for data recording;
class record {
	private:
		int ProgBarCount;		//Counter for progress bar;
		int ProgBarStep;		//One time step for progress bar;
	public:
		bool isRecordingTimeSeries, isRecordingLyap;	//What to record;
		bool isRecordingCVExp;
		bool isRecordingSnapshot, isRecordingSnapshotGS, isRecordingSnapshotCV;
		bool isRecordingSmallestAngle, isRecordingAllAngles;
		bool isRecordingLyapTimeSeriesAlways, isRecordingAngleTimeSeries;
		bool isRecordingDOSViolation, isRecordingNeighAngle;
		bool isRecordingAllDOSViolation, isRecordingStepOrderParam, isRecordingFTLEFluc;
		bool AngleHistSeparate;
		bool isRecordingSnapHist, isRecordingPtcpHist, isRecordingLyapHist;
		bool isRecordingOutputQuantity;
		bool OutputQuantityFirstLastSw;
		bool isRecordingDensityTimeSeries;
		bool SnapshotFirstLastSw;
		bool isDOSFirstTurn;
		bool isRecordingSupplPtcpRatio;
		bool isRecordingSubspaceAngle, isRecordingSubspaceAngleTimeSeries;
		bool LyapHistGSCV;									//if using CV exponent;
		bool RecordTimeSeriesSw;							//true at the moment to record time series;
		bool *SubspaceComplSw;								//if true, indicating complementary set;
		int SnapshotGSNum, SnapshotCVNum;
		int *SnapshotGSInd, *SnapshotCVInd;
		int ExpTimeSeriesNum;
		int *ExpTimeSeriesInd;
		int MaxExpInd, MinContInd;							//Max/Min index of expanding/contracting directions;
		int SmallestAngleBinNum, AllAngleBinNum;			//Number of bins for angle histogram;
		int NeighAngleBinNum, NeighAngleDataNum, AllAngleDataNum;
		int NeighborDistance;
		int AngleTimeSeriesPairNum, AngleTimeSeriesNeighborSw, *AngleTimeSeriesInd;
		int SnapHistBinNum, SnapHistInterval, SnapHistCt;
		int RecordTimeSeriesInt;							//time interval for recording time series
		int NumSnapshot, NumOutputQuantity;
		int DOSPeriodNum, *DOSPeriod, DOSPeriodMax;			//List of periods over which finite-time LEs are overaged
		int DOSLocalExpRateCt, AllDOSDataNum;
		int PtcpHistBinNum, LyapHistBinNum;
		int DensityTimeSeriesBinNum, DensityTimeSeriesTotalBinNum;
		int NumSubspace, NumSubspacePair, SubspaceAngleBinNum;
		int *SubspaceVecInd, *SubspacePairInd;				//Index of vectors spanning subspace and of subspace itself;
		int *SubspaceDim;									//Dimension of subspace;
		unsigned int *SmallestAngleHist, *AllAngleHist;		//Angle histogram;
		unsigned int *NeighAngleHist;
		unsigned int *DOSViolationCount;					//Counting DOS violation;
		unsigned int *AllDOSViolationCount;
		unsigned int *SnapHist;								//Histogram of snapshot;
		unsigned int *PtcpHist, *LyapHist;					//Histogram of participation ratio;
		unsigned int *DensityTimeSeriesCount;
		unsigned int *SubspaceAngleHist;
		double Pi, InvSmallestAngleBinWidth, InvAllAngleBinWidth, InvNeighAngleBinWidth;
		double *NeighAngleMoment, *AllAngleMoment;			//Moments of neighbors-angle;
		double *DOSLocalExpRate, *DOSFiniteTimeLE;			//Stored local expansion rate;
		double *DOSMoment, *FTLECov;						//Moments of step order parameter, covariance of finite-time LE;
		double SnapHistMinX, SnapHistMaxX, InvSnapHistBinWidth;
		double LyapHistMin, LyapHistMax, InvLyapHistBinWidth;
		double PtcpHistMin, PtcpHistMax, InvPtcpHistBinWidth;
		double *LyapPtcpHist;								//Lyapunov exponent for each bin of Y2 histogram;
		double *DensityTimeSeriesMin, *DensityTimeSeriesMax, *InvDensityTimeSeriesBinWidth;
		double *SupplPtcpRatio, *SupplPtcpRatioTmp;
		double InvSubspaceAngleBinWidth;
		ofstream OutFile, OutFileTimeSeries, OutFileLyap;	//Output files;
		ofstream OutFileCVExp, OutFileSnapshot;
		ofstream *OutFileSnapshotGS, *OutFileSnapshotCV;
		ofstream OutFileAngleTimeSeries, OutFileSubspaceAngleTimeSeries;
		ofstream OutFileOutputQuantity, OutFileOutputQuantityBack;
		ofstream OutFileDensityTimeSeries;
		string OutFileName;		//Output file name;
#if ARMASW
		mat *MatSubspace;								//Orthogonal-bases matrix for subspace;
#endif
		~record() {
			//Destructor;
			int i;
			//Close time series files;
			if (isRecordingTimeSeries) OutFileTimeSeries.close();
			if (isRecordingLyap) OutFileLyap.close();
			if (isRecordingCVExp) OutFileCVExp.close();
			if (isRecordingSnapshot) OutFileSnapshot.close();
			if (isRecordingSnapshotGS)
				for (i=0;i<SnapshotGSNum;i++) OutFileSnapshotGS[i].close();
			if (isRecordingSnapshotCV)
				for (i=0;i<SnapshotCVNum;i++) OutFileSnapshotCV[i].close();
			if (isRecordingAngleTimeSeries) OutFileAngleTimeSeries.close();
			if (isRecordingSubspaceAngle && isRecordingSubspaceAngleTimeSeries) OutFileSubspaceAngleTimeSeries.close();
			if (isRecordingOutputQuantity && OutFileOutputQuantity.good()) OutFileOutputQuantity.close();
			if (isRecordingOutputQuantity && OutFileOutputQuantityBack.good()) OutFileOutputQuantityBack.close();
			if (isRecordingDensityTimeSeries) OutFileDensityTimeSeries.close();
			//Free memory;
			if (isRecordingSnapshotGS) delete[] SnapshotGSInd, OutFileSnapshotGS;
			if (isRecordingSnapshotCV) delete[] SnapshotCVInd, OutFileSnapshotCV;
			if (isRecordingLyapTimeSeriesAlways && ExpTimeSeriesNum > 0) delete[] ExpTimeSeriesInd;
			if (isRecordingSmallestAngle) delete[] SmallestAngleHist;
			if (isRecordingAllAngles) delete[] AllAngleHist, AllAngleMoment;
			if (isRecordingAngleTimeSeries && AngleTimeSeriesNeighborSw == 0) delete[] AngleTimeSeriesInd;
			if (isRecordingDOSViolation) {
				delete[] DOSViolationCount, DOSPeriod, DOSLocalExpRate, DOSFiniteTimeLE;
				if (isRecordingAllDOSViolation) delete[] AllDOSViolationCount;
				if (isRecordingStepOrderParam) delete[] DOSMoment;
				if (isRecordingFTLEFluc) delete[] FTLECov;
			}
			if (isRecordingNeighAngle) delete[] NeighAngleHist, NeighAngleMoment;
			if (isRecordingSnapHist) delete[] SnapHist;
			if (isRecordingPtcpHist) delete[] PtcpHist, LyapPtcpHist;
			if (isRecordingLyapHist) delete[] LyapHist;
			if (isRecordingDensityTimeSeries) delete[] DensityTimeSeriesMin, DensityTimeSeriesMax, InvDensityTimeSeriesBinWidth, DensityTimeSeriesCount;
			if (isRecordingSupplPtcpRatio) delete[] SupplPtcpRatio, SupplPtcpRatioTmp;
			if (isRecordingSubspaceAngle) delete[] SubspaceVecInd, SubspacePairInd, SubspaceDim, SubspaceComplSw, SubspaceAngleHist;
#if ARMASW
			if (isRecordingSubspaceAngle) delete[] MatSubspace;
#endif
		}
		void Init(map *Map) {
			//Initialize recording stuff;
			int i,j;
			int i1,i2;
			int MaxInd;
			char buf[80];
			//Set variables for progress bar;
			ProgBarCount = 0;
			switch (Map->BackwardMode) {
				case 0:
				case 2:
				case 3:
					ProgBarStep = Map->RecordEndTime;
					break;
				case 1:
				case 4:
					if (Map->isDividingLoop) ProgBarStep = Map->ThermalizationTime + Map->TransientTime * 3 + Map->RecordPeriod * 3;
						else ProgBarStep = Map->ThermalizationTime + Map->TransientTime * 2 + Map->RecordPeriod * 2;
					break;
				case 5:
					ProgBarStep = Map->RecordPeriod * 2;
					break;
				default:
					ErrorExit("Undefined mode.");
					break;
			}
			if (Map->BackwardMode == 1 && Map->isContinued) {
				if (!(Map->isFromBeginning)) ProgBarStep -= Map->ThermalizationTime + Map->TransientTime + Map->SavedRecordPeriod;
					else ProgBarStep -= Map->ThermalizationTime + Map->TransientTime;
			}
			ProgBarStep /= 100;
			//Initialization for angle distribution;
			Pi = 3.1415926536;
			if (isRecordingSmallestAngle) {
				MaxInd = SmallestAngleBinNum;
				SmallestAngleHist = new unsigned int[MaxInd];
				for (i=0;i<MaxInd;i++) SmallestAngleHist[i] = 0;
				InvSmallestAngleBinWidth = static_cast<double>(SmallestAngleBinNum) / (0.5 * Pi);
			}
			if (isRecordingAllAngles) {
				//Variables for histogram;
				AllAngleDataNum = (MaxExpInd + 1) * (Map->NumLyap - MinContInd) - (MaxExpInd > MinContInd ? (MaxExpInd - MinContInd + 1) * (MaxExpInd - MinContInd + 2) / 2 : 0);
				//AllAngleDataNum = Map->NumLyap * (Map->NumLyap - 1) / 2;
				MaxInd = AllAngleBinNum * (AngleHistSeparate ? AllAngleDataNum : 1);
				AllAngleHist = new unsigned int[MaxInd];
				for (i=0;i<MaxInd;i++) AllAngleHist[i] = 0;
				InvAllAngleBinWidth = static_cast<double>(AllAngleBinNum) / Pi;
				//Variables for moments;
				MaxInd = AllAngleDataNum * 4;
				AllAngleMoment = new double[MaxInd];
				for (i=0;i<MaxInd;i++) AllAngleMoment[i] = 0.0;
			}
			if (isRecordingNeighAngle) {
				//Variables for histogram;
				NeighAngleDataNum = Map->NumLyap - NeighborDistance;
				MaxInd = NeighAngleDataNum * NeighAngleBinNum;
				NeighAngleHist = new unsigned int[MaxInd];
				for (i=0;i<MaxInd;i++) NeighAngleHist[i] = 0;
				InvNeighAngleBinWidth = static_cast<double>(NeighAngleBinNum) / Pi;
				//Variables for moments;
				MaxInd = NeighAngleDataNum * 4;
				NeighAngleMoment = new double[MaxInd];
				for (i=0;i<MaxInd;i++) NeighAngleMoment[i] = 0.0;
			}
			if (isRecordingSubspaceAngle) {
				//Set subspace variables;
				for (i=0;i<NumSubspace;i++) {
					if (SubspaceVecInd[i*2] * SubspaceVecInd[i*2+1] < 0) ErrorExit("error in subspace index");
					SubspaceComplSw[i] = (SubspaceVecInd[i*2] < 0);
					if (SubspaceComplSw[i]) {
						SubspaceVecInd[i*2] = -SubspaceVecInd[i*2]; SubspaceVecInd[i*2+1] = -SubspaceVecInd[i*2+1];
						SubspaceDim[i] = Map->NumLyap - (SubspaceVecInd[i*2+1] - SubspaceVecInd[i*2] + 1);
					} else {
						SubspaceDim[i] = (SubspaceVecInd[i*2+1] - SubspaceVecInd[i*2] + 1);
					}
				}
				//Set pair variables;
				for (i=0;i<NumSubspacePair;i++) {
					i1 = SubspacePairInd[i*2]; i2 = SubspacePairInd[i*2+1];
					if (SubspaceDim[i2] > SubspaceDim[i1]) {
						SubspacePairInd[i*2] = i2;
						SubspacePairInd[i*2+1] = i1;
					}
				}
				//Initialize histogram;
				MaxInd = NumSubspacePair * SubspaceAngleBinNum;
				SubspaceAngleHist = new unsigned int[MaxInd];
				for (i=0;i<MaxInd;i++) SubspaceAngleHist[i] = 0;
				InvSubspaceAngleBinWidth = static_cast<double>(SubspaceAngleBinNum) / (0.5 * Pi);
			}
			//Initialization for snapshot distribution;
			if (isRecordingSnapHist) {
				InvSnapHistBinWidth = static_cast<double>(SnapHistBinNum) / (SnapHistMaxX - SnapHistMinX);
				MaxInd = SnapHistBinNum * SnapHistInterval;
				SnapHist = new unsigned int[MaxInd];
				for (i=0;i<MaxInd;i++) SnapHist[i] = 0;
				SnapHistCt = 0;
			}
			//Initialization for inv. ptcp. ratio distribution;
			if (isRecordingPtcpHist) {
				PtcpHistMin = log(1.0 / static_cast<double>(Map->Size));
				PtcpHistMax = 0.0;	// = log(1.0);
				InvPtcpHistBinWidth = static_cast<double>(PtcpHistBinNum) / (PtcpHistMax - PtcpHistMin);
				MaxInd = PtcpHistBinNum * Map->NumLyap;
				PtcpHist = new unsigned int[MaxInd];
				LyapPtcpHist = new double[MaxInd];
				for (i=0;i<MaxInd;i++) {
					PtcpHist[i] = 0;
					LyapPtcpHist[i] = 0.0;
				}
			}
			//Initialization for Lyap exp distribution;
			if (isRecordingLyapHist) {
				InvLyapHistBinWidth = static_cast<double>(LyapHistBinNum) / (LyapHistMax - LyapHistMin);
				MaxInd = LyapHistBinNum * Map->NumLyap;
				LyapHist = new unsigned int[MaxInd];
				for (i=0;i<MaxInd;i++) LyapHist[i] = 0;
			}
			//Initialization for time series of dynamics density;
			if (isRecordingDensityTimeSeries) {
				//Initialize variables;
				for (i=0;i<Map->LocalDOF;i++)
					InvDensityTimeSeriesBinWidth[i] = static_cast<double>(DensityTimeSeriesBinNum) / (DensityTimeSeriesMax[i] - DensityTimeSeriesMin[i]);
				DensityTimeSeriesTotalBinNum = 1;
				for (i=0;i<Map->LocalDOF;i++) DensityTimeSeriesTotalBinNum *= DensityTimeSeriesBinNum;
				double DensityTimeSeriesFactor = 1.0 / static_cast<double>(Map->RealSize);
				for (i=0;i<Map->LocalDOF;i++) DensityTimeSeriesFactor *= InvDensityTimeSeriesBinWidth[i];
				DensityTimeSeriesCount = new unsigned int[DensityTimeSeriesTotalBinNum];
				//Open file;
				OutFileDensityTimeSeries.open((OutFileName+"-densts.dat").c_str());
				OutFileDensityTimeSeries << setprecision(15);
				OutFileDensityTimeSeries << DensityTimeSeriesFactor;
				for (i=0;i<Map->LocalDOF;i++)
					for (j=0;j<DensityTimeSeriesBinNum;j++)
						OutFileDensityTimeSeries << " " << DensityTimeSeriesMin[i] + (static_cast<double>(j) + 0.5) / InvDensityTimeSeriesBinWidth[i];
				for (i=0;i<DensityTimeSeriesTotalBinNum - Map->LocalDOF * DensityTimeSeriesBinNum;i++)
					OutFileDensityTimeSeries << " " << 0;
				OutFileDensityTimeSeries << endl;
			}
			//Allocate memory;
			if (isRecordingDOSViolation) {
				MaxInd = Map->NumLyap * NUMDOSV * DOSPeriodNum;
				DOSViolationCount = new unsigned int[MaxInd];
				for (i=0;i<MaxInd;i++) DOSViolationCount[i] = 0;
				DOSPeriodMax = DOSPeriod[DOSPeriodNum-1];
				MaxInd = Map->NumLyap * DOSPeriodMax;
				DOSLocalExpRate = new double[MaxInd];
				for (i=0;i<MaxInd;i++) DOSLocalExpRate[i] = 0.0;
				DOSFiniteTimeLE = new double[Map->NumLyap];
				for (i=0;i<Map->NumLyap;i++) DOSFiniteTimeLE[i] = 0.0;
				isDOSFirstTurn = true;
				DOSLocalExpRateCt = 0;
				if (isRecordingAllDOSViolation) {
					AllDOSDataNum = Map->NumLyap * (Map->NumLyap - 1) / 2;
					MaxInd = AllDOSDataNum * DOSPeriodNum;
					AllDOSViolationCount = new unsigned int[MaxInd];
					for (i=0;i<MaxInd;i++) AllDOSViolationCount[i] = 0;
				}
				if (isRecordingStepOrderParam) {
					MaxInd = (Map->NumLyap-2) * 4;
					DOSMoment = new double[MaxInd];
					for (i=0;i<MaxInd;i++) DOSMoment[i] = 0.0; 
				}
				if (isRecordingFTLEFluc) {
					MaxInd = Map->NumLyap * (Map->NumLyap + 1) / 2 * DOSPeriodNum;
					FTLECov = new double[MaxInd];
					for (i=0;i<MaxInd;i++) FTLECov[i] = 0.0;
				}
			}
			if (isRecordingSupplPtcpRatio) {
				MaxInd = Map->NumLyap * NUMSUPPLPTCP;
				SupplPtcpRatio = new double[MaxInd];
				for (i=0;i<MaxInd;i++) SupplPtcpRatio[i] = 0.0;
				SupplPtcpRatioTmp = new double[NUMSUPPLPTCP*3];
			}
			if (Map->isRecordingTsCumulants) {
				MaxInd = 4 * Map->LocalDOF;
				Map->TsCumulants = new double[MaxInd];
				for (i=0;i<MaxInd;i++) Map->TsCumulants[i] = 0.0;
			}
			//Make parameter file;
			OutFile.open((OutFileName+"-param.dat").c_str());
			OutFile << "% Map : " << Map->MapName << endl;
			OutFile << "% Dim : " << Map->Dim << endl;
			OutFile << "% Size : " << Map->RealSize << endl;
			OutFile << "% Degree of freedom : " << Map->Size << endl;
			for (i=0;i<Map->ParamNum;i++) {
				OutFile << "% " << Map->ParamName[i] << " : " << Map->Param[i] << endl;
			}
			OutFile << "% Random seed : " << Map->Seed << endl;
			OutFile << "% Number of calculated Lyapunov vectors : " << Map->NumLyap << endl;
			OutFile << "% Thermalization time : " << Map->ThermalizationTime << endl;
			OutFile << "% Transient time : " << Map->TransientTime << endl;
			OutFile << "% Recording period : " << Map->RecordPeriod << endl;
			OutFile << "% Time interval between orthonormalizing vectors : " << Map->GSInterval << endl;
			OutFile << "% Maximum occupied memory : " << Map->MaxMem << endl; 
			OutFile.close();
			//Make parameter binary file;
			if (Map->BackwardMode == 2) {
				OutFile.open((OutFileName+"-param.bin").c_str(), ios::out | ios::binary | ios::trunc);
				i = VERSION;
				OutFile.write(reinterpret_cast<char*>(&i),sizeof(int));				
				i = Map->MapID;
				OutFile.write(reinterpret_cast<char*>(&i),sizeof(int));
				OutFile.write(reinterpret_cast<char*>(&(Map->Dim)),sizeof(int));
				OutFile.write(reinterpret_cast<char*>(&(Map->RealSize)),sizeof(int));
				OutFile.write(reinterpret_cast<char*>(Map->Param),Map->ParamNum * sizeof(double));
				OutFile.write(reinterpret_cast<char*>(&(Map->Seed)),sizeof(unsigned long));
				OutFile.write(reinterpret_cast<char*>(&(Map->NumLyap)),sizeof(int));
				OutFile.write(reinterpret_cast<char*>(&(Map->ThermalizationTime)),sizeof(int));
				OutFile.write(reinterpret_cast<char*>(&(Map->TransientTime)),sizeof(int));
				OutFile.write(reinterpret_cast<char*>(&(Map->RecordPeriod)),sizeof(int));
				OutFile.write(reinterpret_cast<char*>(&(Map->GSInterval)),sizeof(int));
				OutFile.write(reinterpret_cast<char*>(&(Map->MaxMem)),sizeof(double));
				OutFile.close();
			}
			//Open record files;
			if (isRecordingTimeSeries)
				if (Map->isContinued && !Map->isFromBeginning && isFileExisting((OutFileName+"-ts.dat").c_str())) OutFileTimeSeries.open((OutFileName+"-ts.dat").c_str(), ios::out | ios::app);
					else OutFileTimeSeries.open((OutFileName+"-ts.dat").c_str(), ios::out | ios::trunc);
			if (isRecordingLyap) {
				if (Map->isContinued && !Map->isFromBeginning && isFileExisting((OutFileName+"-lyapexp.dat").c_str())) OutFileLyap.open((OutFileName+"-lyapexp.dat").c_str(), ios::out | ios::app);
					else OutFileLyap.open((OutFileName+"-lyapexp.dat").c_str(), ios::out | ios::trunc);
				OutFileLyap << setprecision(15);
			}
			if (isRecordingCVExp) {
				OutFileCVExp.open((OutFileName+"-cvexp.dat").c_str());
				OutFileCVExp << setprecision(15);
			}
			if (isRecordingSnapshot) OutFileSnapshot.open((OutFileName+"-snapshot.dat").c_str());
			if (isRecordingSnapshotGS) for (i=0;i<SnapshotGSNum;i++) {
				sprintf(buf,"%d",SnapshotGSInd[i]);
				OutFileSnapshotGS[i].open((OutFileName+"-snapshotgs"+buf+".dat").c_str());
			}
			if (isRecordingSnapshotCV) for (i=0;i<SnapshotCVNum;i++) {
				sprintf(buf,"%d",SnapshotCVInd[i]);
				OutFileSnapshotCV[i].open((OutFileName+"-snapshotcv"+buf+".dat").c_str());
			}
			if (isRecordingAngleTimeSeries) OutFileAngleTimeSeries.open((OutFileName+"-anglets.dat").c_str());
			if (isRecordingSubspaceAngle && isRecordingSubspaceAngleTimeSeries) OutFileSubspaceAngleTimeSeries.open((OutFileName+"-subspaceangts.dat").c_str());
			if (Map->isDividingLoop) {
				Map->OutFileTmp.open((OutFileName + ".var").c_str(),ios::in | ios::out | ios::binary | (Map->isContinued ? ios::ate : ios::trunc));
				Map->CurrTmpFileID = 0;
			}
			if (isRecordingOutputQuantity && Map->ExistOutputQuantity) OutFileOutputQuantity.open((OutFileName+"-quantity.dat").c_str());
			if (isRecordingOutputQuantity && Map->ExistOutputQuantityBack) OutFileOutputQuantityBack.open((OutFileName+"-quantityback.dat").c_str());
		}
		void InitDensityTimeSeries(map *Map) {
			DensityTimeSeriesMin = new double[Map->LocalDOF];
			DensityTimeSeriesMax = new double[Map->LocalDOF];
			InvDensityTimeSeriesBinWidth = new double[Map->LocalDOF];
		}
		void WriteResults(map *Map) {
			//Add results to binary parameter file;
			OutFile.open((OutFileName+"-param.bin").c_str(), ios::out | ios::binary | ios::app);
			OutFile.write(reinterpret_cast<char*>(&(Map->LoopInd)), sizeof(int));
			OutFile.write(reinterpret_cast<char*>(Map->LyapExp), Map->NumLyap * sizeof(double));
			OutFile.write(reinterpret_cast<char*>(Map->GSPtcpRatio), Map->NumLyap * sizeof(double));
			OutFile.write(reinterpret_cast<char*>(Map->LyapExpStd), Map->NumLyap * sizeof(double));
			OutFile.close();
		}
		void WriteAfterBackwardTransient(map *Map) {
			//Output data after backward transient;
			OutFile.open((OutFileName+"-abtrans.bin").c_str(), ios::out | ios::binary);
			OutFile.write(reinterpret_cast<char*>(Map->MatC), Map->SizeMatR * sizeof(double));
			OutFile.close();
		}
		void ReadAfterBackwardTransient(map *Map) {
			//Input data after backward transient;
			ifstream InFile;
			InFile.open((OutFileName+"-abtrans.bin").c_str(), ios::in | ios::binary);
			if (!InFile) ErrorExit("Can't continue.");
			InFile.read(reinterpret_cast<char*>(Map->MatC), Map->SizeMatR * sizeof(double));
			InFile.close();
		}
		void RecordTimeSeries(map *Map) {
			//Record time series;
			int i;
			if (!RecordTimeSeriesSw) return;
			OutFileTimeSeries << setprecision(9) << static_cast<double>(Map->t) * Map->TimeStepLength << setprecision(6);
			for (i=0;i<Map->LocalDOF;i++) OutFileTimeSeries << ' ' << Map->MeanX[i];
			OutFileTimeSeries << endl;
		}
		void RecordLyap(map *Map) {
			//Record Lyapunov exponents;
			int i;
			if ((Map->BackwardMode != 2 || isRecordingLyapTimeSeriesAlways) && (!RecordTimeSeriesSw)) return;
			OutFileLyap << static_cast<double>(Map->t) * Map->TimeStepLength;
			if (isRecordingLyapTimeSeriesAlways && ExpTimeSeriesNum > 0) {
				for (i=0;i<ExpTimeSeriesNum;i++)
					OutFileLyap << ' ' << Map->LocalExpRate[ExpTimeSeriesInd[i]] / (Map->TimeStepLength * Map->SignDynamics);
				for (i=0;i<ExpTimeSeriesNum;i++)
					OutFileLyap << ' ' << Map->LocalPtcpRatio[ExpTimeSeriesInd[i]];
			} else {
				for (i=0;i<Map->NumLyap;i++)
					OutFileLyap << ' ' << Map->LyapExp[i] / (static_cast<double>(Map->LyapRecCount) * Map->TimeStepLength * Map->SignDynamics);
				for (i=0;i<Map->NumLyap;i++)
					OutFileLyap << ' ' << Map->GSPtcpRatio[i] / static_cast<double>(Map->LyapRecCount);
				for (i=0;i<Map->NumLyap;i++)
					OutFileLyap << ' ' << Map->LocalExpRate[i] / (Map->TimeStepLength * Map->SignDynamics);
				for (i=0;i<Map->NumLyap;i++)
					OutFileLyap << ' ' << Map->LocalPtcpRatio[i];
			}
			OutFileLyap << endl;
		}
		void RecordCVExp(map *Map) {
			//Record exponents for covariant vectors;
			int i;
			if ((Map->BackwardMode != 2 || isRecordingLyapTimeSeriesAlways) && (!RecordTimeSeriesSw)) return;
			OutFileCVExp << static_cast<double>(Map->t + 1) * Map->TimeStepLength;	//Because of 1 time evolution;
			if (isRecordingLyapTimeSeriesAlways && ExpTimeSeriesNum > 0) {
				if (Map->isCheckingCVExp) for (i=0;i<ExpTimeSeriesNum;i++)
					OutFileCVExp << ' ' << Map->LocalExpRate[ExpTimeSeriesInd[i]] / (Map->TimeStepLength * Map->SignDynamics);
				for (i=0;i<ExpTimeSeriesNum;i++)
					OutFileCVExp << ' ' << Map->LocalPtcpRatio[ExpTimeSeriesInd[i]];
			} else {
				if (Map->isCheckingCVExp) for (i=0;i<Map->NumLyap;i++)
					OutFileCVExp << ' ' << Map->CVExp[i] / (static_cast<double>(Map->CVExpRecCount) * Map->TimeStepLength * Map->SignDynamics);
				for (i=0;i<Map->NumLyap;i++)
					OutFileCVExp << ' ' << Map->CVPtcpRatio[i] / static_cast<double>(Map->CVExpRecCount);
				if (Map->isCheckingCVExp) for (i=0;i<Map->NumLyap;i++)
					OutFileCVExp << ' ' << Map->LocalExpRate[i] / (Map->TimeStepLength * Map->SignDynamics);
				for (i=0;i<Map->NumLyap;i++)
					OutFileCVExp << ' ' << Map->LocalPtcpRatio[i];
			}
			OutFileCVExp << endl;
		}
		void RecordBackCVExp(map *Map) {
			//Record exponents for covariant vectors for inversed dynamics case;
			int i;
			if ((Map->BackwardMode != 2 || isRecordingLyapTimeSeriesAlways) && (!RecordTimeSeriesSw)) return;
			OutFileCVExp << static_cast<double>(Map->t - 1) * Map->TimeStepLength;	//Because of 1 time evolution;
			if (Map->isCheckingCVExp) for (i=0;i<Map->NumLyap;i++)
				OutFileCVExp << ' ' << Map->CVExp[i] / (static_cast<double>(Map->CVExpRecCount) * Map->TimeStepLength * Map->SignDynamics);
			for (i=0;i<Map->NumLyap;i++)
				OutFileCVExp << ' ' << Map->CVPtcpRatio[i] / static_cast<double>(Map->CVExpRecCount);
			if (Map->isCheckingCVExp) for (i=0;i<Map->NumLyap;i++)
				OutFileCVExp << ' ' << Map->LocalExpRate[i] / (Map->TimeStepLength * Map->SignDynamics);
			for (i=0;i<Map->NumLyap;i++)
				OutFileCVExp << ' ' << Map->LocalPtcpRatio[i];
			OutFileCVExp << endl;
		}
		void RecordAngleTimeSeries(map *Map) {
			//Record time-series of angles;
			int i,j;
			int MaxInd1, MaxInd2;
			double *pdX1, *pdX2;
			double DotProd, DotProdMax;
			if (!RecordTimeSeriesSw) return;
			OutFileAngleTimeSeries << setprecision(9) << static_cast<double>(Map->t) * Map->TimeStepLength << setprecision(6);
			switch (AngleTimeSeriesNeighborSw) {
				case 0:
					for (i=0;i<AngleTimeSeriesPairNum;i++) {
						pdX1 = Map->CV + AngleTimeSeriesInd[i*2] * Map->Size;
						pdX2 = Map->CV + AngleTimeSeriesInd[i*2+1] * Map->Size;
						DotProd = Map->VectorDotProduct(pdX1, pdX2);
						DotProd = MIN(MAX(DotProd, -1.0 + DBL_EPSILON), 1.0 - DBL_EPSILON);
						OutFileAngleTimeSeries << ' ' << acos(DotProd);
					}
					break;
				case 1:
					for (i=0;i<Map->NumLyap-1;i++) {
						pdX1 = Map->CV + i * Map->Size;
						pdX2 = Map->CV + (i+1) * Map->Size;
						DotProd = Map->VectorDotProduct(pdX1, pdX2);
						DotProd = MIN(MAX(DotProd, -1.0 + DBL_EPSILON), 1.0 - DBL_EPSILON);
						OutFileAngleTimeSeries << ' ' << acos(DotProd);
					}
					break;
				case 2:
					DotProdMax = 0.0;
					for (i=0;i<Map->NumLyap-1;i++) {
						pdX1 = Map->CV + i * Map->Size;
						for (j=i+1;j<Map->NumLyap;j++) {
							pdX2 = Map->CV + j * Map->Size;
							DotProd = fabs(Map->VectorDotProduct(pdX1, pdX2));
							DotProd = MIN(MAX(DotProd, -1.0 + DBL_EPSILON), 1.0 - DBL_EPSILON);
							if (DotProd > DotProdMax) {
								DotProdMax = DotProd;
								MaxInd1 = i; MaxInd2 = j;
							}
						}
					}
					OutFileAngleTimeSeries << ' ' << acos(DotProdMax) << ' ' << MaxInd1 << ' ' << MaxInd2;
					break;
			}
			OutFileAngleTimeSeries << endl;
		}
		void ProgBar() {
			//Display progress bar;
			ProgBarCount++;
			if (ProgBarCount >= ProgBarStep) {
				ProgBarCount = 0;
				cout << ".";
				cout.flush();
			}
		}
		void RecordResults(map *Map) {
			//Record results of calculation;
			int i,j,k,ct;
			int Sum;
			unsigned int *pNeighAngleHist, *pDOSViolationCount;
			double TmpDbl;
			double InvSum, InvBinWidth;
			double *pMdX;
			double Std, ForthMom;
			double *pMoment, *pFTLECov;
			OutFile.open((OutFileName+"-result.dat").c_str());
			OutFile << setprecision(12);
			for (i=0;i<Map->NumLyap;i++) {
				OutFile << Map->GSPtcpRatio[i] << ' ';
				if (Map->isDoingBackward) OutFile << Map->CVPtcpRatio[i] << ' ';
				OutFile << Map->LyapExp[i] << ' ';
				if (Map->isCheckingCVExp) OutFile << Map->CVExp[i] << ' ';
				if (Map->Ver == -1 && Map->isContinued && !(Map->isFromBeginning)) OutFile << "-1 ";
					else OutFile << Map->LyapExpStd[i] << ' ';
				if (Map->isCheckingCVExp) OutFile << Map->CVExpStd[i] << ' ';
				OutFile << endl;
			}
			OutFile.close();
			//Record angle distributions;
			if (isRecordingSmallestAngle) {
				OutFile.open((OutFileName+"-minang.dat").c_str());
				Sum = 0;
				for (i=0;i<SmallestAngleBinNum;i++) Sum += SmallestAngleHist[i];
				InvSum = 1.0 / static_cast<double>(Sum);
				for (i=0;i<SmallestAngleBinNum;i++) {
					OutFile << (static_cast<double>(i)+0.5) / InvSmallestAngleBinWidth << ' ';
					for (j=0;j<(AngleHistSeparate ? Map->NumLyap : 1);j++) {
						OutFile << SmallestAngleHist[i + j * SmallestAngleBinNum] << ' ';
						OutFile << static_cast<double>(SmallestAngleHist[i + j * SmallestAngleBinNum]) * InvSum * InvSmallestAngleBinWidth << ' ';
					}
					OutFile << endl;
				}
				OutFile.close();
			}
			if (isRecordingAllAngles) {
				//Histogram;
				OutFile.open((OutFileName+"-allang.dat").c_str());
				Sum = 0;
				for (i=0;i<AllAngleBinNum;i++) Sum += AllAngleHist[i];
				InvSum = 1.0 / static_cast<double>(Sum);
				for (i=0;i<AllAngleBinNum;i++) {
					OutFile << (static_cast<double>(i)+0.5) / InvAllAngleBinWidth << ' ';
					for (j=0;j<(AngleHistSeparate ? AllAngleDataNum : 1);j++)
						OutFile << AllAngleHist[i + j * AllAngleBinNum] << ' ';
					for (j=0;j<(AngleHistSeparate ? AllAngleDataNum : 1);j++)
						OutFile << static_cast<double>(AllAngleHist[i + j * AllAngleBinNum]) * InvSum * InvAllAngleBinWidth << ' ';
					OutFile << endl;
				}
				OutFile.close();
				//Moments;
				for (i=0;i<4*AllAngleDataNum;i++) AllAngleMoment[i] /= static_cast<double>(Map->RecordPeriod);
				OutFile.open((OutFileName+"-allangmom.dat").c_str());
				//Column: 1st-4th moments, Standard deviation, 4th moment around center, kurtosis excess; 
				for (i=0;i<AllAngleDataNum;i++) {
					pMoment = AllAngleMoment + i * 4;
					for (j=0;j<4;j++) OutFile << pMoment[j] << ' ';
					Std = pMoment[1] - pMoment[0] * pMoment[0];
					OutFile << Std << ' ';
					ForthMom = pMoment[3] - 4.0 * pMoment[2] * pMoment[0] + 6.0 * pMoment[1] * pMoment[0] * pMoment[0] - 3.0 * pMoment[0] * pMoment[0] * pMoment[0] * pMoment[0];
					OutFile << ForthMom << ' ' << ForthMom / (Std * Std) - 3.0 << endl;
				}
				OutFile.close();
			}
			if (Map->isRecordingTsCumulants) {
				//Cumulants for mean-field time series;
				for (i=0;i<4*(Map->LocalDOF);i++) Map->TsCumulants[i] /= static_cast<double>(Map->RecordPeriod);
				OutFile.open((OutFileName+"-tscum.dat").c_str());
				//Column: 1st-4th moments, Standard deviation, 4th cumulant, Binder's cumulant (kurtosis); 
				for (i=0;i<(Map->LocalDOF);i++) {
					pMoment = Map->TsCumulants + i * 4;
					for (j=0;j<4;j++) OutFile << pMoment[j] << ' ';
					Std = pMoment[1] - pMoment[0] * pMoment[0];
					OutFile << Std << ' ';
					ForthMom = pMoment[3] - 4.0 * pMoment[2] * pMoment[0] + 6.0 * pMoment[1] * pMoment[0] * pMoment[0] - 3.0 * pMoment[0] * pMoment[0] * pMoment[0] * pMoment[0];
					OutFile << ForthMom - 3.0 * Std * Std << ' ';
					OutFile << ForthMom / (Std * Std) - 3.0 << endl;
				}
				OutFile.close();
			}
			if (isRecordingNeighAngle) {
				//Histogram;
				OutFile.open((OutFileName+"-neighang.dat").c_str());
				Sum = 0;
				for (i=0;i<NeighAngleBinNum;i++) Sum += NeighAngleHist[i];
				InvSum = 1.0 / static_cast<double>(Sum);
				for (i=0;i<NeighAngleBinNum;i++) {
					OutFile << (static_cast<double>(i)+0.5) / InvNeighAngleBinWidth << ' ';
					for (j=0;j<NeighAngleDataNum;j++)
						OutFile << NeighAngleHist[j * NeighAngleBinNum + i] << ' ';
					for (j=0;j<NeighAngleDataNum;j++)
						OutFile << static_cast<double>(NeighAngleHist[j * NeighAngleBinNum + i]) * InvSum * InvNeighAngleBinWidth << ' ';
					OutFile << endl;
				}
				OutFile.close();
				//Moments;
				for (i=0;i<4*NeighAngleDataNum;i++) NeighAngleMoment[i] /= static_cast<double>(Map->RecordPeriod);
				OutFile.open((OutFileName+"-neighangmom.dat").c_str());
				//Column: 1st-4th moments, Standard deviation, 4th moment around center, kurtosis excess; 
				for (i=0;i<NeighAngleDataNum;i++) {
					pMoment = NeighAngleMoment + i * 4;
					for (j=0;j<4;j++) OutFile << pMoment[j] << ' ';
					Std = pMoment[1] - pMoment[0] * pMoment[0];
					OutFile << Std << ' ';
					ForthMom = pMoment[3] - 4.0 * pMoment[2] * pMoment[0] + 6.0 * pMoment[1] * pMoment[0] * pMoment[0] - 3.0 * pMoment[0] * pMoment[0] * pMoment[0] * pMoment[0];
					OutFile << ForthMom << ' ' << ForthMom / (Std * Std) - 3.0 << endl;
				}
				OutFile.close();
			}
			if (isRecordingSubspaceAngle) {
				//Histogram;
				OutFile.open((OutFileName+"-subspaceang.dat").c_str());
				Sum = 0;
				for (i=0;i<SubspaceAngleBinNum;i++) Sum += SubspaceAngleHist[i];
				InvSum = 1.0 / static_cast<double>(Sum);
				for (i=0;i<SubspaceAngleBinNum;i++) {
					OutFile << (static_cast<double>(i)+0.5) / InvSubspaceAngleBinWidth << ' ';
					for (j=0;j<NumSubspacePair;j++)
						OutFile << static_cast<double>(SubspaceAngleHist[j * SubspaceAngleBinNum + i]) * InvSum * InvSubspaceAngleBinWidth << ' ';
					for (j=0;j<NumSubspacePair;j++)
						OutFile << SubspaceAngleHist[j * SubspaceAngleBinNum + i] << ' ';
					OutFile << endl;
				}
				OutFile.close();
			}
			//Record DOS Violations;
			if (isRecordingDOSViolation) {
				//v;
				OutFile.open((OutFileName+"-dosv.dat").c_str());
				for (i=0;i<Map->NumLyap;i++) {
					for (k=0;k<DOSPeriodNum;k++)
						for (j=0;j<NUMDOSV;j++)
							OutFile << static_cast<double>(DOSViolationCount[(k*Map->NumLyap+i)*NUMDOSV+j]) / static_cast<double>(Map->RecordPeriod - DOSPeriod[k] + 1) << ' ';
					OutFile << endl;
				}
				OutFile.close();
				//v all pairs;
				if (isRecordingAllDOSViolation) {
					OutFile.open((OutFileName+"-alldosv.dat").c_str());
					for (k=0;k<DOSPeriodNum;k++) {
						pDOSViolationCount = AllDOSViolationCount + k * AllDOSDataNum;
						for (i=0;i<AllDOSDataNum;i++)
							OutFile << static_cast<double>(pDOSViolationCount[i]) / static_cast<double>(Map->RecordPeriod - DOSPeriod[k]) << ' ';
						OutFile << endl;
					}
					OutFile.close();
				}
				//Moments and Binder's cumulant;
				if (isRecordingStepOrderParam) {
					for (i=0;i<(Map->NumLyap-2)*4;i++) DOSMoment[i] /= static_cast<double>(Map->RecordPeriod);
					OutFile.open((OutFileName+"-dosmom.dat").c_str());
					for (i=0;i<Map->NumLyap-2;i++) {
						pMoment = DOSMoment + i * 4;
						for (j=0;j<4;j++) OutFile << pMoment[j] << ' ';
						OutFile << pMoment[3] / (pMoment[1] * pMoment[1]) - 3.0 << ' ';
						Std = pMoment[1] - pMoment[0] * pMoment[0];
						OutFile << Std << ' ';
						ForthMom = pMoment[3] - 4.0 * pMoment[2] * pMoment[0] + 6.0 * pMoment[1] * pMoment[0] * pMoment[0] - 3.0 * pMoment[0] * pMoment[0] * pMoment[0] * pMoment[0];
						OutFile << ForthMom << ' ' << ForthMom / (Std * Std) - 3.0 << endl;
					}
					OutFile.close();
				}
				//FTLE fluctuations;
				if (isRecordingFTLEFluc) {
					OutFile.open((OutFileName+"-ftlefluc.dat").c_str());
					for (k=0;k<DOSPeriodNum;k++) {
						pFTLECov = FTLECov + Map->NumLyap * (Map->NumLyap + 1) / 2 * k;
						TmpDbl = static_cast<double>(DOSPeriod[k]);
						ct=0;
						for (i=0;i<Map->NumLyap;i++) {
							for (j=i;j<Map->NumLyap;j++) {
								OutFile << pFTLECov[ct] / static_cast<double>(Map->RecordPeriod - DOSPeriod[k]) / TmpDbl - Map->CVExp[i] * Map->CVExp[j] * TmpDbl << ' ';
								ct++;
							}
						}
						OutFile << endl;
					}
					OutFile.close();
				}
			}
			//Record Power Spectra;
			Map->RecordPowerSpectra();
			//Record snapshot histogram;
			if (isRecordingSnapHist) {
				OutFile.open((OutFileName+"-snaphist.dat").c_str());
				for (i=0;i<SnapHistBinNum;i++) {
					OutFile << SnapHistMinX + (static_cast<double>(i)+0.5) / InvSnapHistBinWidth << ' ';
					for (j=0;j<SnapHistInterval;j++) {
						OutFile << SnapHist[i + j * SnapHistBinNum] << ' ';
					}
					OutFile << endl;
				}
				OutFile.close();
			}
			//Record inv. ptcp. ratio histogram;
			if (isRecordingPtcpHist) {
				OutFile.open((OutFileName+"-ptcphist.dat").c_str());
				InvSum = 1.0 / static_cast<double>(Map->RecordPeriod);
				for (i=0;i<PtcpHistBinNum;i++) {
					InvBinWidth = InvSum / ((exp(1.0 / InvPtcpHistBinWidth) - 1.0) * exp(PtcpHistMin + static_cast<double>(i) / InvPtcpHistBinWidth));
					OutFile << exp(PtcpHistMin + (static_cast<double>(i) + 0.5) / InvPtcpHistBinWidth) << ' ';
					for (j=0;j<Map->NumLyap;j++)
						OutFile << PtcpHist[i + j * PtcpHistBinNum] << ' ';
					for (j=0;j<Map->NumLyap;j++)
						OutFile << static_cast<double>(PtcpHist[i + j * PtcpHistBinNum]) * InvBinWidth << ' ';
					if (Map->isCheckingCVExp) for (j=0;j<Map->NumLyap;j++)
						OutFile << (PtcpHist[i + j * PtcpHistBinNum] == 0 ? 0 : LyapPtcpHist[i + j * PtcpHistBinNum] / static_cast<double>(PtcpHist[i + j * PtcpHistBinNum]) / (Map->TimeStepLength * Map->SignDynamics)) << ' ';
					OutFile << endl;
				}
				OutFile.close();
			}
			//Record Lyap exp histogram;
			if (isRecordingLyapHist) {
				OutFile.open((OutFileName+"-"+(LyapHistGSCV ? "cv" : "gs")+"lyaphist.dat").c_str());
				InvSum = 1.0 / static_cast<double>(Map->RecordPeriod);
				for (i=0;i<LyapHistBinNum;i++) {
					OutFile << LyapHistMin + (static_cast<double>(i)+0.5) / InvLyapHistBinWidth << ' ';
					for (j=0;j<Map->NumLyap;j++)
						OutFile << LyapHist[i + j * LyapHistBinNum] << ' ';
					for (j=0;j<Map->NumLyap;j++)
						OutFile << static_cast<double>(LyapHist[i + j * LyapHistBinNum]) * InvSum * InvLyapHistBinWidth << ' ';
					OutFile << endl;
				}
				OutFile.close();
			}
			//Record supplementary ptcp. ratio;
			if (isRecordingSupplPtcpRatio) {
				OutFile.open((OutFileName+"-supplptcp.dat").c_str());
				OutFile << setprecision(12);
				for (i=0;i<Map->NumLyap;i++) {
					OutFile << Map->LyapExp[i] << ' ';
					if (Map->isCheckingCVExp) OutFile << Map->CVExp[i] << ' ';
					for (j=0;j<NUMSUPPLPTCP;j++) OutFile << SupplPtcpRatio[i*NUMSUPPLPTCP+j] / static_cast<double>(Map->RecordPeriod) << ' ';
					OutFile << endl;
				}
				OutFile << setprecision(6);
				OutFile.close();
			}
			//Record histogram of map-specific quantity;
			if (Map->isRecordingQuantityHist) {
				OutFile.open((OutFileName+"-quantityhist.dat").c_str());
				InvSum = 1.0 / static_cast<double>(Map->RecordPeriod);
				for (i=0;i<(Map->QuantityHistBinNum);i++) {
					OutFile << Map->QuantityHistMin + (static_cast<double>(i)+0.5) / Map->InvQuantityHistBinWidth << ' ';
					for (j=0;j<(Map->QuantityHistNum);j++)
						OutFile << static_cast<double>(Map->QuantityHist[j*(Map->QuantityHistBinNum)+i]) * InvSum * Map->InvQuantityHistBinWidth << ' ';
					for (j=0;j<(Map->QuantityHistNum);j++)
						OutFile << Map->QuantityHist[j*(Map->QuantityHistBinNum)+i] << ' ';
					for (j=0;j<(Map->QuantityHistSupplNum);j++) {
						ct = Map->QuantityHist[(Map->QuantityHistSupplDistID[j])*(Map->QuantityHistBinNum)+i];
						if (ct > 0)
							OutFile << Map->QuantityHistSuppl[j*(Map->QuantityHistBinNum)+i] / static_cast<double>(ct) << ' ';
						else
							OutFile << "nan ";
					}
					OutFile << endl;
				}
				OutFile.close();
			}
		}
		void RecordSnapshot(map *Map) {
			//Record snapshots of local variables;
			int i;
			if (!SnapshotFirstLastSw && !RecordTimeSeriesSw) return;
			OutFileSnapshot << setprecision(9) << static_cast<double>(Map->t) * Map->TimeStepLength << setprecision(15);
			for (i=0;i<Map->Size;i++) OutFileSnapshot << ' ' << Map->X[i];
			OutFileSnapshot << endl;
		}
		void RecordSnapshotGS(map *Map) {
			//Record snapshots of GS vectors;
			int i,j;
			double *pdX;
			if (!SnapshotFirstLastSw && !RecordTimeSeriesSw) return;
			for (i=0;i<SnapshotGSNum;i++) {
				OutFileSnapshotGS[i] << setprecision(9) << static_cast<double>(Map->t) * Map->TimeStepLength << setprecision(15);
				pdX = Map->dX + SnapshotGSInd[i] * Map->Size;
				for (j=0;j<Map->Size;j++) OutFileSnapshotGS[i] << ' ' << pdX[j];
				OutFileSnapshotGS[i] << endl;
			}
		}
		void RecordSnapshotCV(map *Map) {
			//Record snapshots of Cov. vectors;
			int i,j;
			double *pdX;
			if (!SnapshotFirstLastSw && !RecordTimeSeriesSw) return;
			for (i=0;i<SnapshotCVNum;i++) {
				OutFileSnapshotCV[i] << setprecision(9) << static_cast<double>(Map->t) * Map->TimeStepLength << setprecision(15);
				pdX = Map->CV + SnapshotCVInd[i] * Map->Size;
				for (j=0;j<Map->Size;j++) OutFileSnapshotCV[i] << ' ' << pdX[j];
				OutFileSnapshotCV[i] << endl;
			}
		}
		void RecordOutputQuantity(map *Map) {
			//Record map-specific quantity;
			int i;
			if (!OutputQuantityFirstLastSw && !RecordTimeSeriesSw) return;
			if (!(Map->isNowRecordingOutputQuantity)) return;
			OutFileOutputQuantity << setprecision(9) << static_cast<double>(Map->t) * Map->TimeStepLength << setprecision(6);
			for (i=0;i<Map->OutputQuantityNum;i++) OutFileOutputQuantity << ' ' << Map->OutputQuantity[i];
			OutFileOutputQuantity << endl;
		}
		void RecordOutputQuantityBack(map *Map) {
			//Record map-specific quantity;
			int i;
			if (!OutputQuantityFirstLastSw && !RecordTimeSeriesSw) return;
			if (!(Map->isNowRecordingOutputQuantity)) return;
			OutFileOutputQuantityBack << setprecision(9) << static_cast<double>(Map->t) * Map->TimeStepLength << setprecision(6);
			for (i=0;i<Map->OutputQuantityNumBack;i++) OutFileOutputQuantityBack << ' ' << Map->OutputQuantity[i];
			OutFileOutputQuantityBack << endl;
		}
		void InitSnapshotGS() {
			//Initialization for snapshots of GS vectors;
			SnapshotGSInd = new int[SnapshotGSNum];
			OutFileSnapshotGS = new ofstream[SnapshotGSNum];
		} 
		void InitSnapshotCV() {
			//Initialization for snapshots of Cov. vectors;
			SnapshotCVInd = new int[SnapshotCVNum];
			OutFileSnapshotCV = new ofstream[SnapshotCVNum];
		}
		void InitSubspaceAngle() {
			//Initialization for subspace angle;
#if ARMASW
			SubspaceVecInd = new int[NumSubspace * 2];
			SubspacePairInd = new int[NumSubspacePair * 2];
			SubspaceDim = new int[NumSubspace];
			SubspaceComplSw = new bool[NumSubspace];
			MatSubspace = new mat[NumSubspace];
#endif
		} 
		void InitDOS() {
			//Initialization for DOS;
			DOSPeriod = new int[DOSPeriodNum];
		}
		void InitAngleTimeSeries() {
			//Initialization for angle time series;
			AngleTimeSeriesInd = new int[AngleTimeSeriesPairNum*2];
		}
		void AccAngleHist(map *Map) {
			//Accumulation for angle histograms;
			int i,j,k;
			double *pdX1, *pdX2;
			double DotProd, Angle, AnglePow;
			double SmallestAngle = 0.5 * Pi;
			unsigned int *pAllAngleHist;
			double *pMoment;
			//Set pointers;
			pAllAngleHist = AllAngleHist;
			pMoment = AllAngleMoment;
			//Calculate angles;
			for (i=0;i<=MaxExpInd;i++) {
				pdX1 = Map->CV + i * Map->Size;
				for (j=MAX(i+1,MinContInd);j<Map->NumLyap;j++) {
					pdX2 = Map->CV + j * Map->Size;
					DotProd = Map->VectorDotProduct(pdX1, pdX2);
					DotProd = MIN(MAX(DotProd, -1.0 + DBL_EPSILON), 1.0 - DBL_EPSILON);
					//Put this angle to histogram;
					Angle = acos(DotProd);
					if (isRecordingAllAngles) {
						pAllAngleHist[static_cast<int>(Angle * InvAllAngleBinWidth)]++;
						//Count for moment;
						AnglePow = 1.0;
						for (k=0;k<4;k++) {
							AnglePow *= Angle;
							pMoment[k] += AnglePow;
						}
						pMoment += 4;
					}
					Angle = acos(fabs(DotProd));
					if (Angle < SmallestAngle) SmallestAngle = Angle;
					if (AngleHistSeparate) pAllAngleHist += AllAngleBinNum;
				}
			}
			//Put smallest angle to histogram;
			if (isRecordingSmallestAngle) SmallestAngleHist[static_cast<int>(SmallestAngle * InvSmallestAngleBinWidth)]++;
		}
		void AccNeighAngleHist(map *Map) {
			//Accumulation for neighbors-angle histograms;
			int i,j,k;
			double *pdX1, *pdX2;
			double DotProd, Angle;
			double AnglePow;
			//Calculate angles;
			for (i=0;i<NeighAngleDataNum;i++) {
				pdX1 = Map->CV + i * Map->Size;
				pdX2 = Map->CV + (i + NeighborDistance) * Map->Size;
				DotProd = Map->VectorDotProduct(pdX1, pdX2);
				DotProd = MAX(MIN(DotProd, 1.0 - DBL_EPSILON), -1.0 + DBL_EPSILON);
				Angle = acos(DotProd);
				//Put this angle to histogram;
				NeighAngleHist[static_cast<int>(Angle * InvNeighAngleBinWidth) + i * NeighAngleBinNum]++;
				//Count for moment;
				AnglePow = 1.0;
				for (j=0;j<4;j++) {
					AnglePow *= Angle;
					NeighAngleMoment[i * 4 + j] += AnglePow;
				}
			}
		}
		void CountDOSViolation(map *Map) {
			//Calculate if DOS is violated or not;
			//Considered ratios;
			//local(i;i+1), local(i;i+d), local(i;i-1,i+1), local(i;i-d,i+d);
			//semi-global(i;i+2,i+4,...,i-2,i-4,...);
			//global(i;all except one in the same pair), global(i;all);
			int i,j,k,ip,ct;
			int NowPeriod, NowInd;
			unsigned int *pDOSViolationCount;
			double m, m2;
			double *pDOSLocalExpRate;
			double *pDOSMoment, *pFTLECov;
			bool DOSViolated;
			//Store local expansion rate & clear finite-time LEs;
			pDOSLocalExpRate = DOSLocalExpRate + DOSLocalExpRateCt * Map->NumLyap;
			for (i=0;i<Map->NumLyap;i++) {
				pDOSLocalExpRate[i] = Map->LocalExpRate[i];
				DOSFiniteTimeLE[i] = 0.0;
			}
			//Check if DOS is violated;
			NowPeriod = 0;
#if OMPSW
			#pragma omp parallel private(i,j,k,pDOSViolationCount,DOSViolated)
			for (k=0;k<DOSPeriodNum;k++) {
				if (isDOSFirstTurn && DOSLocalExpRateCt + 1 < DOSPeriod[k]) continue; 
				//Calculate finite-time LEs;
				#pragma omp single
				{
				for (;NowPeriod < DOSPeriod[k];NowPeriod++) {
					NowInd = DOSLocalExpRateCt - NowPeriod + (DOSLocalExpRateCt < NowPeriod ? DOSPeriodMax : 0);
					pDOSLocalExpRate = DOSLocalExpRate + NowInd * Map->NumLyap;
					for (i=0;i<Map->NumLyap;i++) DOSFiniteTimeLE[i] += pDOSLocalExpRate[i];
				}
				}
				#pragma omp for
				for (i=0;i<Map->NumLyap;i++) {
					pDOSViolationCount = DOSViolationCount + (k * Map->NumLyap + i) * NUMDOSV;
					//local(i;i+1), local(i;i+d), local(i;i-1,i+1), local(i;i-d,i+d);
					if (i<Map->NumLyap-1 && DOSFiniteTimeLE[i] < DOSFiniteTimeLE[i+1]) pDOSViolationCount[0]++;
					if (i<Map->NumLyap-NeighborDistance && DOSFiniteTimeLE[i] < DOSFiniteTimeLE[i+NeighborDistance]) pDOSViolationCount[1]++;
					if (i<Map->NumLyap-2 && (DOSFiniteTimeLE[i+1] < DOSFiniteTimeLE[i+2] || DOSFiniteTimeLE[i+1] > DOSFiniteTimeLE[i])) pDOSViolationCount[2]++;
					if (i<Map->NumLyap-2*NeighborDistance && (DOSFiniteTimeLE[i+NeighborDistance] < DOSFiniteTimeLE[i+2*NeighborDistance] || DOSFiniteTimeLE[i+NeighborDistance] > DOSFiniteTimeLE[i])) pDOSViolationCount[3]++;
					//semi-global(i;i+2,i+4,...,i-2,i-4,...);
					DOSViolated = false;
					for (j=(i%NeighborDistance);j<i;j+=NeighborDistance) DOSViolated = DOSViolated || (DOSFiniteTimeLE[i] > DOSFiniteTimeLE[j]);
					for (j=i+NeighborDistance;j<Map->NumLyap;j+=NeighborDistance) DOSViolated = DOSViolated || (DOSFiniteTimeLE[i] < DOSFiniteTimeLE[j]);
					if (DOSViolated) pDOSViolationCount[4]++;
					//global(i;all except one in the same pair), global(i;all);
					DOSViolated = false;
					for (j=0;j<i-1;j++) DOSViolated = DOSViolated || (DOSFiniteTimeLE[i] > DOSFiniteTimeLE[j]);
					for (j=i+2;j<Map->NumLyap;j++) DOSViolated = DOSViolated || (DOSFiniteTimeLE[i] < DOSFiniteTimeLE[j]);
					if (DOSViolated || (i<Map->NumLyap-1 ? DOSFiniteTimeLE[i] < DOSFiniteTimeLE[i+1] : false)) pDOSViolationCount[5+(i%2)]++;
					if (DOSViolated || (i>0 ? DOSFiniteTimeLE[i] > DOSFiniteTimeLE[i-1] : false)) pDOSViolationCount[6-(i%2)]++;
					//if (i!=0 && i!=Map->NumLyap-1) DOSViolated = DOSViolated || (i>ip ? Map->LocalExpRate[i] > Map->LocalExpRate[ip] : Map->LocalExpRate[i] < Map->LocalExpRate[ip]);
					DOSViolated = DOSViolated || (i<Map->NumLyap-1 ? DOSFiniteTimeLE[i] < DOSFiniteTimeLE[i+1] : false) || (i>0 ? DOSFiniteTimeLE[i] > DOSFiniteTimeLE[i-1] : false);
					if (DOSViolated) pDOSViolationCount[7]++;
					//DOS for all pairs;
					if (isRecordingAllDOSViolation) {
						pDOSViolationCount = AllDOSViolationCount + k * AllDOSDataNum + i * Map->NumLyap - i*(i+1)/2;
						for (j=i+1;j<Map->NumLyap;j++)
							if (DOSFiniteTimeLE[i] < DOSFiniteTimeLE[j]) pDOSViolationCount[j-(i+1)]++;
					}
				}
				#pragma omp single
				{
				if (isRecordingFTLEFluc) {
					pFTLECov = FTLECov + Map->NumLyap * (Map->NumLyap + 1) / 2 * k;
					ct=0;
					for (i=0;i<Map->NumLyap;i++) {
						for (j=i;j<Map->NumLyap;j++) {
							pFTLECov[ct] += DOSFiniteTimeLE[i] * DOSFiniteTimeLE[j];
							ct++;
						}
					}
				}
				}
			}
#else
			for (k=0;k<DOSPeriodNum;k++) {
				if (isDOSFirstTurn && DOSLocalExpRateCt + 1 < DOSPeriod[k]) continue; 
				//Calculate finite-time LEs;
				for (;NowPeriod < DOSPeriod[k];NowPeriod++) {
					NowInd = DOSLocalExpRateCt - NowPeriod + (DOSLocalExpRateCt < NowPeriod ? DOSPeriodMax : 0);
					pDOSLocalExpRate = DOSLocalExpRate + NowInd * Map->NumLyap;
					for (i=0;i<Map->NumLyap;i++) DOSFiniteTimeLE[i] += pDOSLocalExpRate[i];
				}
				for (i=0;i<Map->NumLyap;i++) {
					pDOSViolationCount = DOSViolationCount + (k * Map->NumLyap + i) * NUMDOSV;
					//local(i;i+1), local(i;i+d), local(i;i-1,i+1), local(i;i-d,i+d);
					if (i<Map->NumLyap-1 && DOSFiniteTimeLE[i] < DOSFiniteTimeLE[i+1]) pDOSViolationCount[0]++;
					if (i<Map->NumLyap-NeighborDistance && DOSFiniteTimeLE[i] < DOSFiniteTimeLE[i+NeighborDistance]) pDOSViolationCount[1]++;
					if (i<Map->NumLyap-2 && (DOSFiniteTimeLE[i+1] < DOSFiniteTimeLE[i+2] || DOSFiniteTimeLE[i+1] > DOSFiniteTimeLE[i])) pDOSViolationCount[2]++;
					if (i<Map->NumLyap-2*NeighborDistance && (DOSFiniteTimeLE[i+NeighborDistance] < DOSFiniteTimeLE[i+2*NeighborDistance] || DOSFiniteTimeLE[i+NeighborDistance] > DOSFiniteTimeLE[i])) pDOSViolationCount[3]++;
					//semi-global(i;i+2,i+4,...,i-2,i-4,...);
					DOSViolated = false;
					for (j=(i%NeighborDistance);j<i;j+=NeighborDistance) DOSViolated = DOSViolated || (DOSFiniteTimeLE[i] > DOSFiniteTimeLE[j]);
					for (j=i+NeighborDistance;j<Map->NumLyap;j+=NeighborDistance) DOSViolated = DOSViolated || (DOSFiniteTimeLE[i] < DOSFiniteTimeLE[j]);
					if (DOSViolated) pDOSViolationCount[4]++;
					//global(i;all except one in the same pair), global(i;all);
					DOSViolated = false;
					for (j=0;j<i-1;j++) DOSViolated = DOSViolated || (DOSFiniteTimeLE[i] > DOSFiniteTimeLE[j]);
					for (j=i+2;j<Map->NumLyap;j++) DOSViolated = DOSViolated || (DOSFiniteTimeLE[i] < DOSFiniteTimeLE[j]);
					if (DOSViolated || (i<Map->NumLyap-1 ? DOSFiniteTimeLE[i] < DOSFiniteTimeLE[i+1] : false)) pDOSViolationCount[5+(i%2)]++;
					if (DOSViolated || (i>0 ? DOSFiniteTimeLE[i] > DOSFiniteTimeLE[i-1] : false)) pDOSViolationCount[6-(i%2)]++;
					//if (i!=0 && i!=Map->NumLyap-1) DOSViolated = DOSViolated || (i>ip ? Map->LocalExpRate[i] > Map->LocalExpRate[ip] : Map->LocalExpRate[i] < Map->LocalExpRate[ip]);
					DOSViolated = DOSViolated || (i<Map->NumLyap-1 ? DOSFiniteTimeLE[i] < DOSFiniteTimeLE[i+1] : false) || (i>0 ? DOSFiniteTimeLE[i] > DOSFiniteTimeLE[i-1] : false);
					if (DOSViolated) pDOSViolationCount[7]++;
					//DOS for all pairs;
					if (isRecordingAllDOSViolation) {
						pDOSViolationCount = AllDOSViolationCount + k * AllDOSDataNum + i * Map->NumLyap - i*(i+1)/2;
						for (j=i+1;j<Map->NumLyap;j++)
							if (DOSFiniteTimeLE[i] < DOSFiniteTimeLE[j]) pDOSViolationCount[j-(i+1)]++;
					}
				}
				if (isRecordingFTLEFluc) {
					pFTLECov = FTLECov + Map->NumLyap * (Map->NumLyap + 1) / 2 * k;
					ct=0;
					for (i=0;i<Map->NumLyap;i++) {
						for (j=i;j<Map->NumLyap;j++) {
							pFTLECov[ct] += DOSFiniteTimeLE[i] * DOSFiniteTimeLE[j];
							ct++;
						}
					}
				}
			}
#endif
			//aftertreatment;
			DOSLocalExpRateCt++;
			if (DOSLocalExpRateCt >= DOSPeriodMax) {
				isDOSFirstTurn = false;
				DOSLocalExpRateCt = 0;
			}
			//Calculate moments of DOS order-parameter;
			// Order-parameter = N|\lambda_{i+1} - 2\lambda_i + \lambda_{i-1}|;
			if (isRecordingStepOrderParam) {
				for (i=0;i<Map->NumLyap-2;i++) {
					pDOSMoment = DOSMoment + i * 4;
					m = fabs(Map->LocalExpRate[i] - 2.0 * Map->LocalExpRate[i+1] + Map->LocalExpRate[i+2]) * static_cast<double>(Map->RealSize);
					m2 = 1.0;
					for (j=0;j<3;j++) {
						m2 *= m;
						pDOSMoment[j] += m2;
					}
				}
			}
		}
		void AccSnapHist(map *Map) {
			//Accumulation for snapshot histograms;
			int i,j,k;
			unsigned int *pSnapHist;
			int BinId;
			//Accumulate histogram;
			pSnapHist = SnapHist + SnapHistCt * SnapHistBinNum;
			for (i=0;i<Map->Size;i++) {
				BinId = static_cast<int>((Map->X[i]-SnapHistMinX)*InvSnapHistBinWidth);
				if (BinId < 0 || BinId >= SnapHistBinNum)
					cout << "Out of range in accumulating snapshot histogram." << endl;
				pSnapHist[BinId]++;
			}
			//Increment counter;
			SnapHistCt++;
			if (SnapHistCt == SnapHistInterval) SnapHistCt = 0;
		}
		void AccPtcpHist(map *Map) {
			//Accumulation for inv. ptcp. ratio histograms;
			int i,j,k;
			unsigned int *pPtcpHist;
			double *pLyapPtcpHist;
			//Accumulate histogram;
			for (i=0;i<Map->NumLyap;i++) {
				pPtcpHist = PtcpHist + i * PtcpHistBinNum;
				pLyapPtcpHist = LyapPtcpHist + i * PtcpHistBinNum;
				pPtcpHist[static_cast<int>((log(Map->LocalPtcpRatio[i])-PtcpHistMin)*InvPtcpHistBinWidth)]++;
				if (Map->isCheckingCVExp) pLyapPtcpHist[static_cast<int>((log(Map->LocalPtcpRatio[i])-PtcpHistMin)*InvPtcpHistBinWidth)] += Map->LocalExpRate[i];
			}
		}
		void AccLyapHist(map *Map) {
			//Accumulation for Lyap exp histograms;
			int i,j,k;
			unsigned int *pLyapHist;
			//Accumulate histogram;
			for (i=0;i<Map->NumLyap;i++) {
				pLyapHist = LyapHist + i * LyapHistBinNum;
				pLyapHist[static_cast<int>((Map->LocalExpRate[i] / (Map->TimeStepLength * Map->SignDynamics) - LyapHistMin) * InvLyapHistBinWidth)]++;
			}
		}
		void RecordDensityTimeSeries(map *Map) {
			//Record time series of density of dynamical variables;
			int i,j,ind;
			//Count;
			for (i=0;i<DensityTimeSeriesTotalBinNum;i++)
				DensityTimeSeriesCount[i] = 0;
			for (i=0;i<Map->Size;i+=Map->LocalDOF) {
				ind = 0;
				for (j=0;j<Map->LocalDOF;j++) ind = ind * DensityTimeSeriesBinNum + static_cast<int>((Map->X[i+j] - DensityTimeSeriesMin[j]) * InvDensityTimeSeriesBinWidth[j]);
				DensityTimeSeriesCount[ind]++;				
			}
			//Record it;
			OutFileDensityTimeSeries << static_cast<double>(Map->t) * Map->TimeStepLength;
			for (i=0;i<DensityTimeSeriesTotalBinNum;i++)
				OutFileDensityTimeSeries << " " << DensityTimeSeriesCount[i];
			OutFileDensityTimeSeries << endl;
			//Attention! Dimension is aligned inversely (e.g., y,x);
		}
		void CalcSupplPtcpRatio(map *Map) {
			//Calculate supplementary participation ratio;
			//1-norm, 4-norm, Morriss entropy, log-averaged 2-norm;
			int i,j,k;
			double a;
			double *ppdX;
			double SupplPtcpRatio1[3];
			double SupplPtcpRatio2[3];
			double SupplPtcpRatioLoc[3];
			for (i=0;i<Map->NumLyap;i++) {
				for (j=0;j<3;j++) {SupplPtcpRatio1[j] = 0.0; SupplPtcpRatio2[j] = 0.0;}
				ppdX = Map->CV + i * (Map->Size);
				for (j=0;j<Map->Size;j+=Map->LocalDOF) {
					for (k=0;k<3;k++) SupplPtcpRatioLoc[k] = 0.0;
					for (k=0;k<Map->LocalDOF;k++) {
						a = fabs(ppdX[j+k]);
						SupplPtcpRatioLoc[0] += a;
						a = a*a;
						SupplPtcpRatioLoc[1] += a;
						a = a*a;
						SupplPtcpRatioLoc[2] += a;
					}
					for (k=0;k<3;k++) {
						SupplPtcpRatio1[k] += (k != 1 ? SupplPtcpRatioLoc[k] : -SupplPtcpRatioLoc[k] * log(SupplPtcpRatioLoc[k]));
						SupplPtcpRatio2[k] += SupplPtcpRatioLoc[k] * SupplPtcpRatioLoc[k];
					};
				}
				SupplPtcpRatio[i*NUMSUPPLPTCP] += SupplPtcpRatio2[0] / (SupplPtcpRatio1[0] * SupplPtcpRatio1[0]);
				SupplPtcpRatio[i*NUMSUPPLPTCP+1] += SupplPtcpRatio2[2] / (SupplPtcpRatio1[2] * SupplPtcpRatio1[2]);
				SupplPtcpRatio[i*NUMSUPPLPTCP+2] += SupplPtcpRatio1[1];
				SupplPtcpRatio[i*NUMSUPPLPTCP+3] += log(SupplPtcpRatio2[1]);
			}
		}
#if ARMASW
		mat SubspaceMatrix(map *Map, int StartInd, int EndInd, int NumInd, bool ComplSw) {
			//Construct (Dim)d-orthogonal-bases matrix for subspace spanned by CV #(StartInd)-#(EndInd);
			//See TeX note for subspace-angle;
			int i,j,k,n;
			int Ind;
			double Norm;
			double *pMatC1;
			mat MatTmp, MatR, MatQ;
			MatTmp.zeros(Map->NumLyap, NumInd);
			//Calculate orthogonal-bases matrix;
			if (ComplSw) {
				for (j=0;j<StartInd;j++) {
					pMatC1 = Map->MatC + j*(j+1)/2;
					for (i=0;i<=j;i++) MatTmp.at(i,j) = pMatC1[i];
				}
				for (Ind=EndInd+1;Ind<(Map->NumLyap);Ind++) {
					j = StartInd + (Ind-EndInd-1);
					pMatC1 = Map->MatC + Ind*(Ind+1)/2;
					for (i=0;i<=Ind;i++) MatTmp.at(i,j) = pMatC1[i];
				}
			} else {
				for (j=0;j<NumInd;j++) {
					Ind = StartInd + j;
					//Orthogonalize;
					pMatC1 = Map->MatC + Ind*(Ind+1)/2;
					for (i=0;i<=Ind;i++) MatTmp.at(i,j) = pMatC1[i];
				}
			}
			qr(MatQ, MatR, MatTmp);
			return MatQ.cols(0,NumInd-1);
		}
#endif
		void AccSubspaceAngleHist(map *Map) {
			//Accumulate histogram for subspace angles;
#if ARMASW
			int i,j;
			double Ang;
			mat MatTmp;
			colvec SingVal;
			//Record time in time series;
			if (isRecordingSubspaceAngleTimeSeries)
				OutFileSubspaceAngleTimeSeries << setprecision(9) << static_cast<double>(Map->t) * Map->TimeStepLength << setprecision(6);
			//construct orthogonal-bases matrix for each subspace;
			for (i=0;i<NumSubspace;i++)
				MatSubspace[i] = SubspaceMatrix(Map, SubspaceVecInd[i*2], SubspaceVecInd[i*2+1], SubspaceDim[i], SubspaceComplSw[i]);
			//calculate angles and accumulate them;
			for (i=0;i<NumSubspacePair;i++) {
				MatTmp = trans(MatSubspace[SubspacePairInd[i*2]]) * MatSubspace[SubspacePairInd[i*2+1]];
				SingVal = svd(MatTmp);
				Ang = acos(MIN(SingVal[0], 1.0));
				//Put this angle to histogram;
				SubspaceAngleHist[static_cast<int>(Ang * InvSubspaceAngleBinWidth) + i * SubspaceAngleBinNum]++;
				if (isRecordingSubspaceAngleTimeSeries) OutFileSubspaceAngleTimeSeries << ' ' << Ang;
			}
			if (isRecordingSubspaceAngleTimeSeries) OutFileSubspaceAngleTimeSeries << endl;
#endif
		}
		void SetRecordTimeSeriesSw(map *Map) {
			//Set if the moment to record time series or not;
			RecordTimeSeriesSw = (Map->t % RecordTimeSeriesInt == 0);
		}
		void SetRecordTimeSeriesSw(map *Map, int dt) {
			//Set if the moment to record time series or not;
			RecordTimeSeriesSw = ((Map->t + dt) % RecordTimeSeriesInt == 0);
		}
};

//Set map;
void SetMap(map **Map, int MapID) {
	switch (MapID) {
		case 23:
			//1D Kuramoto-Sivashinsky PDE (periodic boundary conditions);
			*Map = new mapKSPDEp;
			break;
		case 76:
			//1D Kuramoto-Sivashinsky PDE for UPO (periodic boundary conditions);
			*Map = new mapKSPDEpUPO;
			break;
		case 77:
			//1D Kuramoto-Sivashinsky PDE with UPO (periodic boundary conditions);
			*Map = new mapKSPDEpUPO2;
			break;
		case 84:
			//1D Kuramoto-Sivashinsky PDE for RPO (periodic boundary conditions);
			*Map = new mapKSPDEpRPO;
			break;
		default:
			//Error;
			ErrorExit("Undefined map.");
	}
}


//Read parameters from saved data;
void ReadParameters(int *MapID, map **Map, record *Record) {
	//Continue from saved data;
	int Ver;
	ifstream *InFileTmp = new ifstream;
	InFileTmp->open((Record->OutFileName+"-param.bin").c_str(), ios::in | ios::binary);
	InFileTmp->read(reinterpret_cast<char*>(&Ver), sizeof(int));
	if (Ver >= 0) {
		//old version;
		Ver = 0;
		InFileTmp->seekg(0, ios::beg);
	}
	InFileTmp->read(reinterpret_cast<char*>(MapID), sizeof(int));
	SetMap(Map, *MapID);
	if (Ver <= -1) InFileTmp->read(reinterpret_cast<char*>(&((*Map)->Dim)), sizeof(int));
		else (*Map)->Dim = 1;
	InFileTmp->read(reinterpret_cast<char*>(&((*Map)->RealSize)), sizeof(int));
	InFileTmp->read(reinterpret_cast<char*>((*Map)->Param), (*Map)->ParamNum * sizeof(double));
	InFileTmp->read(reinterpret_cast<char*>(&((*Map)->Seed)), sizeof(unsigned long));
	InFileTmp->read(reinterpret_cast<char*>(&((*Map)->NumLyap)), sizeof(int));
	if (Ver <= -4) 	InFileTmp->read(reinterpret_cast<char*>(&((*Map)->ThermalizationTime)), sizeof(int));
		else (*Map)->ThermalizationTime = 0;
	InFileTmp->read(reinterpret_cast<char*>(&((*Map)->TransientTime)), sizeof(int));
	InFileTmp->read(reinterpret_cast<char*>(&((*Map)->SavedRecordPeriod)), sizeof(int));
	InFileTmp->read(reinterpret_cast<char*>(&((*Map)->GSInterval)),sizeof(int));
	InFileTmp->read(reinterpret_cast<char*>(&((*Map)->MaxMem)),sizeof(double));
	InFileTmp->peek(); (*Map)->isFromBeginning = InFileTmp->eof();
	if (!((*Map)->isFromBeginning)) {
		InFileTmp->read(reinterpret_cast<char*>(&((*Map)->SavedLoopInd)),sizeof(int));
		(*Map)->InFile = InFileTmp;
	} else {
		if ((*Map)->BackwardMode == 5) ErrorExit("Can't continue.");
		cout << "Can't continue from the end of data. Starting from the beginning..." << endl;
	}
	(*Map)->Ver = Ver;
}

//Input parameters;
void InputParameters(int argc, char *argv[], map **Map, record *Record) {
	//Define variables;
	bool isContinued;
	int i;
	int MapID, Ind, TmpInt;
	unsigned int OutputMode;
	//Input parameters;
	if (argc == 1) {
		//Input from keyboard;
		cout << endl;
		cout << "0 = Start from beginning, 1 = Continue from saved data : ";
		cin >> isContinued;
		if (!isContinued) {
			//Start from beginning;
			cout << "Which map? (int)" << endl;
			cout << "   23 = 1D Kuramoto-Sivashinsky PDE (periodic boundary)" << endl;
			cout << "   76 = 1D Kuramoto-Sivashinsky PDE for UPO (periodic boundary)" << endl;
			cout << "   77 = 1D Kuramoto-Sivashinsky PDE with UPO (periodic boundary)" << endl;
			cout << "   84 = 1D Kuramoto-Sivashinsky PDE for RPO (periodic boundary)" << endl;
			cout << " : ";
			cin >> MapID;
			SetMap(Map, MapID);
			cout << "Spacial dimensions? (int) : ";
			cin >> (*Map)->Dim;
			cout << "Number of elements? (int) : ";
			cin >> (*Map)->RealSize;
			for (i=0;i<(*Map)->ParamNum;i++) {
				cout << (*Map)->ParamName[i] << "? (double) : ";
				cin >> (*Map)->Param[i];
			}
			cout << "Random seed? (unsigned long) : ";
			cin >> (*Map)->Seed;
			cout << "Number of Lyapunov vectors to calculate? (int) : ";
			cin >> (*Map)->NumLyap;
			cout << "Thermalization time? (int, +=per DOF, -=total) : ";
			cin >> (*Map)->ThermalizationTime;
			cout << "Transient time? (int, +=per DOF, -=total) : ";
			cin >> (*Map)->TransientTime;
			cout << "Recording period? (int) : ";
			cin >> (*Map)->RecordPeriod;
			cout << "Time interval between orthonormalizing vectors? (int) : ";
			cin >> (*Map)->GSInterval;
			cout << "Record time series? (0=no, 1=yes) : ";
			cin >> Record->isRecordingTimeSeries;
			cout << "Record Lyapunov exponents? (0=no, 1=yes) : ";
			cin >> Record->isRecordingLyap;
			cout << "Output mode? (int)" << endl;
			cout << "     0 = nothing" << endl;
			cout << "     1 = snapshot" << endl;
			cout << "     2 = gssnapshot" << endl;
			cout << "     4 = cvsnapshot" << endl;
			cout << "     8 = cumulants for time series" << endl;
			cout << "    16 = smallest angle, btwn expanding & contracting directions" << endl;
			cout << "    32 = all angles, btwn expanding & contracting directions" << endl;
			cout << "    64 = record time-series of exponents at every time step" << endl;
			cout << "   128 = angles at every time step" << endl;
			cout << "   256 = violation of DOS" << endl;
			cout << "   512 = angles between neighbors" << endl;
			cout << "  1024 = power spectra" << endl;
			cout << "  2048 = histogram of snapshot (only for # of local DOF = 1)" << endl;
			cout << "  4096 = histogram of inv. ptcp. ratio" << endl;
			cout << "  8192 = finite-amplitude perturbation" << endl;
			cout << " 16384 = histogram of instantaneous Lyap exponents" << endl;
			cout << " 32768 = power spectrum of dynamics" << endl;
			cout << " 65536 = output map-specific quantity" << endl;
			cout << "131072 = time series of density of dynamical variables" << endl;
			cout << "262144 = supplementary definitions for ptcp. ratio" << endl;
			cout << "524288 = angles between subspaces" << endl;
			cout << "1048576= histogram of map-specific quantity" << endl;
			cout << " : ";
			cin >> OutputMode;
			cout << "Backward Mode? (int)" << endl;
			cout << "    0 = no backward no data save" << endl;
			cout << "    1 = backward" << endl;
			cout << "    2 = no backward but data save" << endl;
			cout << "    3 = inverse dynamics" << endl;
			cout << "    4 = inverse dynamics + covariant vectors" << endl;
			cout << "    5 = backward from saved data" << endl;
			cout << " : ";
			cin >> (*Map)->BackwardMode;
			cout << "Check exponents for covariant vectors? (0=no, 1=yes) : ";
			cin >> (*Map)->isCheckingCVExp; 
			cout << "Record them? (0=no, 1=yes) : ";
			cin >> Record->isRecordingCVExp;
			cout << "Output file name? (string, without extention) : ";
			cin >> Record->OutFileName;
			cout << "Maximum occupied memory? (double, +=MB, -=GB, 0.0=default) : ";
			cin >> (*Map)->MaxMem;
			(*Map)->Ver = VERSION;
		} else {
			//Continue from saved data;
			cout << "File name? (string, without extention) : ";
			cin >> Record->OutFileName;
			ReadParameters(&MapID, Map, Record);
			if (!((*Map)->isFromBeginning))
				cout << "Already recorded over " << (*Map)->SavedRecordPeriod << " time steps." << endl;
			cout << "New recording period? (int) : ";
			cin >> (*Map)->RecordPeriod;
			cout << "Record time series? (0=no, 1=yes) : ";
			cin >> Record->isRecordingTimeSeries;
			cout << "Record Lyapunov exponents? (0=no, 1=yes) : ";
			cin >> Record->isRecordingLyap;
			cout << "Output mode? (int)" << endl;
			cout << "     0 = nothing" << endl;
			cout << "     1 = snapshot" << endl;
			cout << "     2 = gssnapshot" << endl;
			cout << "     4 = cvsnapshot" << endl;
			cout << "     8 = cumulants for time series" << endl;
			cout << "    16 = smallest angle, btwn expanding & contracting directions" << endl;
			cout << "    32 = all angles, btwn expanding & contracting directions" << endl;
			cout << "    64 = record time-series of exponents at every time step" << endl;
			cout << "   128 = angles at every time step" << endl;
			cout << "   256 = violation of DOS" << endl;
			cout << "   512 = angles between neighbors" << endl;
			cout << "  1024 = power spectra" << endl;
			cout << "  2048 = histogram of snapshot (only for # of local DOF = 1)" << endl;
			cout << "  4096 = histogram of inv. ptcp. ratio" << endl;
			cout << "  8192 = not used" << endl;
			cout << " 16384 = histogram of instantaneous Lyap exponents" << endl;
			cout << " 32768 = power spectrum of dynamics" << endl;
			cout << " 65536 = output map-specific quantity" << endl;
			cout << "131072 = time series of density of dynamical variables" << endl;
			cout << "262144 = supplementary definitions for ptcp. ratio" << endl;
			cout << "524288 = angles between subspaces" << endl;
			cout << "1048576= histogram of map-specific quantity" << endl;
			cout << " : ";
			cin >> OutputMode;
			cout << "Backward Mode? (int)" << endl;
			cout << "    0 = no backward no data save" << endl;
			cout << "    1 = backward" << endl;
			cout << "    2 = no backward but data save" << endl;
			cout << "    3 = inverse dynamics" << endl;
			cout << "    4 = inverse dynamics + covariant vectors" << endl;
			cout << "    5 = backward from saved data" << endl;
			cout << " : ";
			cin >> (*Map)->BackwardMode;
			cout << "Check exponents for covariant vectors? (0=no, 1=yes) : ";
			cin >> (*Map)->isCheckingCVExp; 
			cout << "Record them? (0=no, 1=yes) : ";
			cin >> Record->isRecordingCVExp;
		}
		//Convert outputmode;
		Record->isRecordingSnapshot = static_cast<bool>(OutputMode & 1);
		Record->isRecordingSnapshotGS = static_cast<bool>(OutputMode & 2);
		Record->isRecordingSnapshotCV = static_cast<bool>(OutputMode & 4);
		(*Map)->isRecordingTsCumulants = static_cast<bool>(OutputMode & 8);
		Record->isRecordingSmallestAngle = static_cast<bool>(OutputMode & 16);
		Record->isRecordingAllAngles = static_cast<bool>(OutputMode & 32);
		Record->isRecordingLyapTimeSeriesAlways = static_cast<bool>(OutputMode & 64);
		Record->isRecordingAngleTimeSeries = static_cast<bool>(OutputMode & 128);
		Record->isRecordingDOSViolation = static_cast<bool>(OutputMode & 256);
		Record->isRecordingNeighAngle = static_cast<bool>(OutputMode & 512);
		(*Map)->isRecordingPowerSpectra = static_cast<bool>(OutputMode & 1024);
		Record->isRecordingSnapHist = static_cast<bool>(OutputMode & 2048);
		Record->isRecordingPtcpHist = static_cast<bool>(OutputMode & 4096);
		Record->isRecordingLyapHist = static_cast<bool>(OutputMode & 16384);
		(*Map)->isRecordingPowSpecDyn = static_cast<bool>(OutputMode & 32768);
		Record->isRecordingOutputQuantity = static_cast<bool>(OutputMode & 65536);
		(*Map)->isRecordingOutputQuantity = Record->isRecordingOutputQuantity;
		Record->isRecordingDensityTimeSeries = static_cast<bool>(OutputMode & 131072);
		Record->isRecordingSupplPtcpRatio = static_cast<bool>(OutputMode & 262144);
		Record->isRecordingSubspaceAngle = static_cast<bool>(OutputMode & 524288);
		(*Map)->isRecordingQuantityHist = static_cast<bool>(OutputMode & 1048576);
		cout << "Time interval for recording time series? (int) : ";
		cin >> Record->RecordTimeSeriesInt;
		if (Record->isRecordingSnapshot || Record->isRecordingSnapshotGS || Record->isRecordingSnapshotCV) {
			cout << "Record the whole time-series (0) or only the first and last snapshots (1) ? (int) : ";
			cin >> Record->SnapshotFirstLastSw;
			if (Record->SnapshotFirstLastSw) {
				cout << "How many snapshots? (int) : ";
				cin >> Record->NumSnapshot;
			}
		}
		if (Record->isRecordingSnapshotGS) {
			cout << "How many GS vectors to take snapshots? (int) : ";
			cin >> Record->SnapshotGSNum;
			Record->InitSnapshotGS();
			for (i=0;i<Record->SnapshotGSNum;i++) {
				cout << "Index of GS vector #" << i << " (int) : ";
				cin >> Record->SnapshotGSInd[i];
			}
		}
		if (Record->isRecordingSnapshotCV) {
			cout << "How many Cov. vectors to take snapshots? (int) : ";
			cin >> Record->SnapshotCVNum;
			Record->InitSnapshotCV();
			for (i=0;i<Record->SnapshotCVNum;i++) {
				cout << "Index of Cov. vector #" << i << " (int) : ";
				cin >> Record->SnapshotCVInd[i];
			}
		}
		if (Record->isRecordingAllAngles) {
			cout << "Scan all directions? (0=no, 1=yes) : ";
			cin >> Record->AngleHistSeparate;
			if (!Record->AngleHistSeparate) {
				cout << "Maximum index of expanding directions? (int) : ";
				cin >> Record->MaxExpInd;
				cout << "Minimum index of contracting directions? (int) : ";
				cin >> Record->MinContInd;
			} else {
				Record->MaxExpInd = (*Map)->NumLyap - 1;
				Record->MinContInd = 0;
			}
		} else if (Record->isRecordingSmallestAngle) {
			Record->AngleHistSeparate = false;
			cout << "Maximum index of expanding directions? (int) : ";
			cin >> Record->MaxExpInd;
			cout << "Minimum index of contracting directions? (int) : ";
			cin >> Record->MinContInd;
		}
		if (Record->isRecordingSmallestAngle) {
			cout << "How many bins for smallest-angle distribution? (int) : ";
			cin >> Record->SmallestAngleBinNum;
		}
		if (Record->isRecordingAllAngles) {
			cout << "How many bins for all-angle distribution? (int) : ";
			cin >> Record->AllAngleBinNum;
		}
		if (Record->isRecordingLyapTimeSeriesAlways) {
			cout << "How many exponents to record? (int, 0=all) : ";
			cin >> Record->ExpTimeSeriesNum;
			if (Record->ExpTimeSeriesNum > 0) {
				Record->ExpTimeSeriesInd = new int[Record->ExpTimeSeriesNum];
				for (i=0;i<Record->ExpTimeSeriesNum;i++) {
					cout << "Index of exponent #" << i << " (int) : ";
					cin >> Record->ExpTimeSeriesInd[i];
				}
			}
		}
		if (Record->isRecordingAngleTimeSeries) {
			cout << "Specify pairs = 0, All neighbors = 1, Minimum angle = 2 : ";
			cin >> Record->AngleTimeSeriesNeighborSw;
			if (Record->AngleTimeSeriesNeighborSw == 0) {
				cout << "How many pairs? (int) : ";
				cin >> Record->AngleTimeSeriesPairNum;
				Record->InitAngleTimeSeries();
				for (i=0;i<Record->AngleTimeSeriesPairNum;i++) {
					cout << "Pair #" << i << ", first vector (int) : ";
					cin >> Record->AngleTimeSeriesInd[i*2];
					cout << "Pair #" << i << ", second vector (int) : ";
					cin >> Record->AngleTimeSeriesInd[i*2+1];
				}
			}
		}
		if (Record->isRecordingDOSViolation || Record->isRecordingNeighAngle) {
			cout << "Distance to neighbor to count? (int) : ";
			cin >> Record->NeighborDistance;
		}
		if (Record->isRecordingDOSViolation) {
			(*Map)->isCheckingCVExp = true;
			cout << "How many periods to use? (int) : ";
			cin >> Record->DOSPeriodNum;
			Record->InitDOS();
			for (i=0;i<Record->DOSPeriodNum;i++) {
				cout << "Period #" << i << " of finite-time LE for DOS (int) : ";
				cin >> Record->DOSPeriod[i];
			}
			cout << "Record all-pair DOS (1) + step order parameter (2) + finite-time LE fluctuations (4)? : ";
			cin >> TmpInt;
			Record->isRecordingAllDOSViolation = static_cast<bool>(TmpInt & 1);
			Record->isRecordingStepOrderParam = static_cast<bool>(TmpInt & 2);
			Record->isRecordingFTLEFluc = static_cast<bool>(TmpInt & 4);
		}
		if (Record->isRecordingNeighAngle) {
			cout << "How many bins for neighbors-angle distribution? (int) : ";
			cin >> Record->NeighAngleBinNum;
		}
		if (Record->isRecordingSnapHist) {
			cout << "How many bins for snapshot histograms? (int) : ";
			cin >> Record->SnapHistBinNum;
			cout << "Min of x for snapshot histogram? (double) : ";
			cin >> Record->SnapHistMinX;
			cout << "Max of x for snapshot histogram? (double) : ";
			cin >> Record->SnapHistMaxX;
			cout << "How many histograms to make? (int) : ";
			cin >> Record->SnapHistInterval;
		}
		if (Record->isRecordingPtcpHist) {
			cout << "How many bins for snapshot histograms? (int) : ";
			cin >> Record->PtcpHistBinNum;
		}
		if (Record->isRecordingLyapHist) {
			cout << "0 = GS exponent, 1 = CV exponent : ";
			cin >> Record->LyapHistGSCV;
			cout << "How many bins for Lyapunov histograms? (int) : ";
			cin >> Record->LyapHistBinNum;
			cout << "Min of Lyap exp. for histogram? (double) : ";
			cin >> Record->LyapHistMin;
			cout << "Max of Lyap exp. for histogram? (double) : ";
			cin >> Record->LyapHistMax;
		}
		if (Record->isRecordingOutputQuantity) {
			cout << "Record the whole time-series (0) or only the first and last parts (1) ? (int) : ";
			cin >> Record->OutputQuantityFirstLastSw;
			if (Record->OutputQuantityFirstLastSw) {
				cout << "How many time steps to record? (int) : ";
				cin >> Record->NumOutputQuantity;
			}
		}
		if (Record->isRecordingDensityTimeSeries) {
			Record->InitDensityTimeSeries(*Map);
			cout << "How many bins for each dimension? (int) : ";
			cin >> Record->DensityTimeSeriesBinNum;
			for (i=0;i<(*Map)->LocalDOF;i++) {
				cout << "Min in dimension " << i << " ? (double) : ";
				cin >> Record->DensityTimeSeriesMin[i];
				cout << "Max in dimension " << i << " ? (double) : ";
				cin >> Record->DensityTimeSeriesMax[i];
			}
		}
		if (Record->isRecordingSubspaceAngle) {
			cout << "How many subspaces? (int) : ";
			cin >> Record->NumSubspace;
			cout << "How many pairs of subspaces? (int) : ";
			cin >> Record->NumSubspacePair;
			Record->InitSubspaceAngle();
			for (i=0;i<Record->NumSubspace;i++) {
				cout << "First index of subspace #" << i << "? (int; negative = complementary subspace) : ";
				cin >> Record->SubspaceVecInd[i*2];
				cout << "Last index of subspace #" << i << "? (int; negative = complementary subspace) : ";
				cin >> Record->SubspaceVecInd[i*2+1];
			}
			for (i=0;i<Record->NumSubspacePair;i++) {
				cout << "First subspace of pair #" << i << "? (int) : ";
				cin >> Record->SubspacePairInd[i*2];
				cout << "Second subspace of pair #" << i << "? (int) : ";
				cin >> Record->SubspacePairInd[i*2+1];
			}
			cout << "How many bins for subspace-angle distributions? (int) : ";
			cin >> Record->SubspaceAngleBinNum;
			cout << "Record time-series as well? (bool) : ";
			cin >> Record->isRecordingSubspaceAngleTimeSeries;
		}
		if ((*Map)->isRecordingQuantityHist) {
			cout << "How many bins for histogram of map-specific quantity? (int) : ";
			cin >> (*Map)->QuantityHistBinNum;
			cout << "Min of x for histogram of map-specific quantity? (double) : ";
			cin >> (*Map)->QuantityHistMin;
			cout << "Max of x for histogram of map-specific quantity? (double) : ";
			cin >> (*Map)->QuantityHistMax;
			cout << "Recording interval for the histogram? 0 for the total histogram. (int) : ";
			cin >> (*Map)->QuantityHistRecInt;
		}
		//Clear buffer;
		cin.ignore(80, '\n');
		cout << endl;
	} else {
		//Set parameters from arguments;
		Ind = 1;
		isContinued = static_cast<bool>(atoi(argv[Ind++]));
		if (!isContinued) {
			//Start from beginning;
			MapID = atoi(argv[Ind++]);
			SetMap(Map, MapID);
			(*Map)->Dim = atoi(argv[Ind++]);
			(*Map)->RealSize = atoi(argv[Ind++]);
			for (i=0;i<(*Map)->ParamNum;i++) (*Map)->Param[i] = atof(argv[Ind++]);
			(*Map)->Seed = static_cast<unsigned long>(atol(argv[Ind++]));
			(*Map)->NumLyap = atoi(argv[Ind++]);
			(*Map)->ThermalizationTime = atoi(argv[Ind++]);
			(*Map)->TransientTime = atoi(argv[Ind++]);
			(*Map)->RecordPeriod = atoi(argv[Ind++]);
			(*Map)->GSInterval = atoi(argv[Ind++]);
			Record->isRecordingTimeSeries = static_cast<bool>(atoi(argv[Ind++]));
			Record->isRecordingLyap = static_cast<bool>(atoi(argv[Ind++]));
			OutputMode = atoi(argv[Ind++]);
			(*Map)->BackwardMode = atoi(argv[Ind++]);
			(*Map)->isCheckingCVExp = static_cast<bool>(atoi(argv[Ind++]));
			Record->isRecordingCVExp = static_cast<bool>(atoi(argv[Ind++]));
			Record->OutFileName = argv[Ind++];
			(*Map)->MaxMem = atof(argv[Ind++]);
			(*Map)->Ver = VERSION;
		} else {
			//Continue from saved data;
			Record->OutFileName = argv[Ind++];
			ReadParameters(&MapID, Map, Record);
			(*Map)->RecordPeriod = atoi(argv[Ind++]);
			Record->isRecordingTimeSeries = static_cast<bool>(atoi(argv[Ind++]));
			Record->isRecordingLyap = static_cast<bool>(atoi(argv[Ind++]));
			OutputMode = atoi(argv[Ind++]);
			(*Map)->BackwardMode = atoi(argv[Ind++]);
			(*Map)->isCheckingCVExp = static_cast<bool>(atoi(argv[Ind++]));
			Record->isRecordingCVExp = static_cast<bool>(atoi(argv[Ind++]));
		}			
		//Convert outputmode;
		Record->isRecordingSnapshot = static_cast<bool>(OutputMode & 1);
		Record->isRecordingSnapshotGS = static_cast<bool>(OutputMode & 2);
		Record->isRecordingSnapshotCV = static_cast<bool>(OutputMode & 4);
		(*Map)->isRecordingTsCumulants = static_cast<bool>(OutputMode & 8);
		Record->isRecordingSmallestAngle = static_cast<bool>(OutputMode & 16);
		Record->isRecordingAllAngles = static_cast<bool>(OutputMode & 32);
		Record->isRecordingLyapTimeSeriesAlways = static_cast<bool>(OutputMode & 64);
		Record->isRecordingAngleTimeSeries = static_cast<bool>(OutputMode & 128);
		Record->isRecordingDOSViolation = static_cast<bool>(OutputMode & 256);
		Record->isRecordingNeighAngle = static_cast<bool>(OutputMode & 512);
		(*Map)->isRecordingPowerSpectra = static_cast<bool>(OutputMode & 1024);
		Record->isRecordingSnapHist = static_cast<bool>(OutputMode & 2048);
		Record->isRecordingPtcpHist = static_cast<bool>(OutputMode & 4096);
		Record->isRecordingLyapHist = static_cast<bool>(OutputMode & 16384);
		(*Map)->isRecordingPowSpecDyn = static_cast<bool>(OutputMode & 32768);
		Record->isRecordingOutputQuantity = static_cast<bool>(OutputMode & 65536);
		(*Map)->isRecordingOutputQuantity = Record->isRecordingOutputQuantity;
		Record->isRecordingDensityTimeSeries = static_cast<bool>(OutputMode & 131072);
		Record->isRecordingSupplPtcpRatio = static_cast<bool>(OutputMode & 262144);
		Record->isRecordingSubspaceAngle = static_cast<bool>(OutputMode & 524288);
		(*Map)->isRecordingQuantityHist = static_cast<bool>(OutputMode & 1048576);
		Record->RecordTimeSeriesInt = atoi(argv[Ind++]);
		if (Record->isRecordingSnapshot || Record->isRecordingSnapshotGS || Record->isRecordingSnapshotCV) {
			Record->SnapshotFirstLastSw = static_cast<bool>(atoi(argv[Ind++]));
			if (Record->SnapshotFirstLastSw) {
				Record->NumSnapshot = atoi(argv[Ind++]);
			}
		}
		if (Record->isRecordingSnapshotGS) {
			Record->SnapshotGSNum = atoi(argv[Ind++]);
			Record->InitSnapshotGS();
			for (i=0;i<Record->SnapshotGSNum;i++) {
				Record->SnapshotGSInd[i] = atoi(argv[Ind++]);
			}
		}
		if (Record->isRecordingSnapshotCV) {
			Record->SnapshotCVNum = atoi(argv[Ind++]);
			Record->InitSnapshotCV();
			for (i=0;i<Record->SnapshotCVNum;i++) {
				Record->SnapshotCVInd[i] = atoi(argv[Ind++]);
			}
		}
		if (Record->isRecordingAllAngles) {
			Record->AngleHistSeparate = static_cast<bool>(atoi(argv[Ind++]));
			if (!Record->AngleHistSeparate) {
				Record->MaxExpInd = atoi(argv[Ind++]);
				Record->MinContInd = atoi(argv[Ind++]);
			} else {
				Record->MaxExpInd = (*Map)->NumLyap - 1;
				Record->MinContInd = 0;
			}
		} else if (Record->isRecordingSmallestAngle) {
			Record->AngleHistSeparate = false;
			Record->MaxExpInd = atoi(argv[Ind++]);
			Record->MinContInd = atoi(argv[Ind++]);
		}
		if (Record->isRecordingSmallestAngle) {
			Record->SmallestAngleBinNum = atoi(argv[Ind++]);
		}
		if (Record->isRecordingAllAngles) {
			Record->AllAngleBinNum = atoi(argv[Ind++]);
		}
		if (Record->isRecordingLyapTimeSeriesAlways) {
			Record->ExpTimeSeriesNum = atoi(argv[Ind++]);
			if (Record->ExpTimeSeriesNum > 0) {
				Record->ExpTimeSeriesInd = new int[Record->ExpTimeSeriesNum];
				for (i=0;i<Record->ExpTimeSeriesNum;i++) {
					Record->ExpTimeSeriesInd[i] = atoi(argv[Ind++]);
				}
			}
		}
		if (Record->isRecordingAngleTimeSeries) {
			Record->AngleTimeSeriesNeighborSw = atoi(argv[Ind++]);
			if (Record->AngleTimeSeriesNeighborSw == 0) {
				Record->AngleTimeSeriesPairNum = atoi(argv[Ind++]);
				Record->InitAngleTimeSeries();
				for (i=0;i<Record->AngleTimeSeriesPairNum;i++) {
					Record->AngleTimeSeriesInd[i*2] = atoi(argv[Ind++]);
					Record->AngleTimeSeriesInd[i*2+1] = atoi(argv[Ind++]);
				}
			}
		}
		if (Record->isRecordingDOSViolation || Record->isRecordingNeighAngle) {
			Record->NeighborDistance = atoi(argv[Ind++]);
		}
		if (Record->isRecordingDOSViolation) {
			(*Map)->isCheckingCVExp = true;
			Record->DOSPeriodNum = atoi(argv[Ind++]);
			Record->InitDOS();
			for (i=0;i<Record->DOSPeriodNum;i++) Record->DOSPeriod[i] = atoi(argv[Ind++]);
			TmpInt = atoi(argv[Ind++]);
			Record->isRecordingAllDOSViolation = static_cast<bool>(TmpInt & 1);
			Record->isRecordingStepOrderParam = static_cast<bool>(TmpInt & 2);
			Record->isRecordingFTLEFluc = static_cast<bool>(TmpInt & 4);
		}
		if (Record->isRecordingNeighAngle) {
			Record->NeighAngleBinNum = atoi(argv[Ind++]);
		}
		if (Record->isRecordingSnapHist) {
			Record->SnapHistBinNum = atoi(argv[Ind++]);
			Record->SnapHistMinX = atof(argv[Ind++]);
			Record->SnapHistMaxX = atof(argv[Ind++]);
			Record->SnapHistInterval = atoi(argv[Ind++]);
		}
		if (Record->isRecordingPtcpHist)
			Record->PtcpHistBinNum = atoi(argv[Ind++]);
		if (Record->isRecordingLyapHist) {
			Record->LyapHistGSCV = static_cast<bool>(atoi(argv[Ind++]));
			Record->LyapHistBinNum = atoi(argv[Ind++]);
			Record->LyapHistMin = atof(argv[Ind++]);
			Record->LyapHistMax = atof(argv[Ind++]);
		}
		if (Record->isRecordingOutputQuantity) {
			Record->OutputQuantityFirstLastSw = static_cast<bool>(atoi(argv[Ind++]));
			if (Record->OutputQuantityFirstLastSw)
				Record->NumOutputQuantity = atoi(argv[Ind++]);
		}
		if (Record->isRecordingDensityTimeSeries) {
			Record->InitDensityTimeSeries(*Map);
			Record->DensityTimeSeriesBinNum = atoi(argv[Ind++]);
			for (i=0;i<(*Map)->LocalDOF;i++) {
				Record->DensityTimeSeriesMin[i] = atof(argv[Ind++]);
				Record->DensityTimeSeriesMax[i] = atof(argv[Ind++]);
			}
		}
		if (Record->isRecordingSubspaceAngle) {
			Record->NumSubspace = atoi(argv[Ind++]);
			Record->NumSubspacePair = atoi(argv[Ind++]);
			Record->InitSubspaceAngle();
			for (i=0;i<Record->NumSubspace;i++) {
				Record->SubspaceVecInd[i*2] = atoi(argv[Ind++]);
				Record->SubspaceVecInd[i*2+1] = atoi(argv[Ind++]);
			}
			for (i=0;i<Record->NumSubspacePair;i++) {
				Record->SubspacePairInd[i*2] = atoi(argv[Ind++]);
				Record->SubspacePairInd[i*2+1] = atoi(argv[Ind++]);
			}
			Record->SubspaceAngleBinNum = atoi(argv[Ind++]);
			Record->isRecordingSubspaceAngleTimeSeries = static_cast<bool>(atoi(argv[Ind++]));
		}
		if ((*Map)->isRecordingQuantityHist) {
			(*Map)->QuantityHistBinNum = atoi(argv[Ind++]);
			(*Map)->QuantityHistMin = atof(argv[Ind++]);
			(*Map)->QuantityHistMax = atof(argv[Ind++]);
			(*Map)->QuantityHistRecInt = atoi(argv[Ind++]);
		}
	}
	cout << Record->OutFileName << endl;
	//Check required sw;
	if (Record->isRecordingSubspaceAngle && !ARMASW) ErrorExit("Armadillo is required for subspace-angle");
	//Convert values;
	(*Map)->isContinued = isContinued;
	if (!((*Map)->BackwardMode == 1 || (*Map)->BackwardMode == 4 || (*Map)->BackwardMode == 5)) {
		(*Map)->isCheckingCVExp = false;
		Record->isRecordingCVExp = false;
	}
	if (!(*Map)->isContinued) {
		if ((*Map)->MaxMem < 0) (*Map)->MaxMem = fabs((*Map)->MaxMem) * 1024.0;
			else if ((*Map)->MaxMem == 0.0) (*Map)->MaxMem = static_cast<double>(DEFAULTMEMORY);
		(*Map)->MaxMem *= 1024.0 * 1024.0;
	} else {
		(*Map)->RecordPeriod += (*Map)->SavedRecordPeriod;
	}
	//Initialize map & recording;
	(*Map)->OutFileName = Record->OutFileName;
	(*Map)->Init();
	Record->Init(*Map);
}


//Main routine without backward evolution
void MainRoutineWithoutBackward(map *Map, record *Record) {
	//Thermalization;
	//Pretend to have NumLyap=0;
	int NumLyapBak = Map->NumLyap;
	Map->NumLyap = 0;
	while (Map->t < Map->ThermalizationTime) {
		//Time evolution;
		Map->f();
		//Display progress bar;
		Record->ProgBar();
	}
	//Initial transient;
	Map->NumLyap = NumLyapBak;
	Map->NextGSTime = Map->ThermalizationTime + Map->GSInterval;
	while (Map->t < Map->TransientEndTime) {
		//Time evolution;
		Map->f();
		//Gram-Schmidt orthonormalization;
		if (Map->t == Map->NextGSTime) {
			Map->GramSchmidt();
			Map->NextGSTime += Map->GSInterval;
		}
		//Display progress bar;
		Record->ProgBar();
	}
	//Time evolution;
	while (Map->t < Map->RecordEndTime) {
		//Write variables;
		if (Map->BackwardMode == 2 && Map->t == Map->LoopTimes[Map->LoopInd]) Map->WriteVar();
		//Time evolution;
		Map->f();
		Record->SetRecordTimeSeriesSw(Map);
		//Time series record;
		if (Record->isRecordingTimeSeries || Map->isRecordingTsCumulants) Map->CalcMeanField();
		if (Record->isRecordingTimeSeries) Record->RecordTimeSeries(Map);
		if (Map->isRecordingTsCumulants) Map->CalcTsCumulants();
		if (Record->isRecordingSnapHist) Record->AccSnapHist(Map);
		if (Map->isRecordingPowSpecDyn) Map->AccPowSpecDyn();
		//Gram-Schmidt orthonormalization;
		Map->GramSchmidt();
		//Record it;
		Map->CalcLyap();
		Map->CalcPtcpRatio(Map->dX, Map->GSPtcpRatio);
		if (Record->isRecordingLyap && (Record->isRecordingLyapTimeSeriesAlways || !(Map->BackwardMode == 2 && Map->t != (Map->LoopTimes)[Map->LoopInd])))
			Record->RecordLyap(Map);
		if (Record->isRecordingLyapHist && !Record->LyapHistGSCV) Record->AccLyapHist(Map);
		if (Record->isRecordingSnapshot && (!Record->SnapshotFirstLastSw || (Map->t <= Map->TransientEndTime+Record->NumSnapshot || Map->t > Map->RecordEndTime-Record->NumSnapshot))) Record->RecordSnapshot(Map);
		if (Record->isRecordingSnapshotGS && (!Record->SnapshotFirstLastSw || (Map->t <= Map->TransientEndTime+Record->NumSnapshot || Map->t > Map->RecordEndTime-Record->NumSnapshot))) Record->RecordSnapshotGS(Map);
		if (Record->isRecordingDensityTimeSeries) Record->RecordDensityTimeSeries(Map);
		if (Record->isRecordingOutputQuantity && Map->ExistOutputQuantity && (!Record->OutputQuantityFirstLastSw || (Map->t <= Map->TransientEndTime+Record->NumOutputQuantity || Map->t > Map->RecordEndTime-Record->NumOutputQuantity))) {Map->CalcOutputQuantity(); Record->RecordOutputQuantity(Map);}
		if (Map->isRecordingQuantityHist && !(Map->isQuantityHistBack)) Map->AccQuantityHist();
		//Display progress bar;
		Record->ProgBar();
	}
	//Write final values of variables;
	if (Map->BackwardMode == 2) Map->WriteVar();
}

//Main routine without loop division
void MainRoutineWithoutDivision(map *Map, record *Record) {
	//Thermalization;
	//Pretend to have NumLyap=0;
	int NumLyapBak = Map->NumLyap;
	Map->NumLyap = 0;
	while (Map->t < Map->ThermalizationTime) {
		//Time evolution;
		Map->f();
		//Display progress bar;
		Record->ProgBar();
	}
	//Initial transient;
	Map->NumLyap = NumLyapBak;
	Map->NextGSTime = Map->ThermalizationTime + Map->GSInterval;
	while (Map->t < Map->TransientEndTime) {
		//Time evolution;
		Map->f();
		//Gram-Schmidt orthonormalization;
		if (Map->t == Map->NextGSTime) {
			Map->GramSchmidt();
			Map->NextGSTime += Map->GSInterval;
		}
		//Display progress bar;
		Record->ProgBar();
	}
	//Time evolution;
	while (Map->t < Map->RecordEndTime) {
		//Time evolution;
		Map->f();
		Record->SetRecordTimeSeriesSw(Map);
		//Time series record;
		if (Record->isRecordingTimeSeries || Map->isRecordingTsCumulants) Map->CalcMeanField();
		if (Record->isRecordingTimeSeries) Record->RecordTimeSeries(Map);
		if (Map->isRecordingTsCumulants) Map->CalcTsCumulants();
		if (Record->isRecordingSnapHist) Record->AccSnapHist(Map);
		if (Map->isRecordingPowSpecDyn) Map->AccPowSpecDyn();
		//Gram-Schmidt orthonormalization;
		Map->GramSchmidt();
		//Record it;
		Map->CalcLyap();
		Map->CalcPtcpRatio(Map->dX, Map->GSPtcpRatio);
		if (Record->isRecordingLyap) Record->RecordLyap(Map);
		if (Record->isRecordingLyapHist && !Record->LyapHistGSCV) Record->AccLyapHist(Map);
		if (Record->isRecordingSnapshot && (!Record->SnapshotFirstLastSw || (Map->t <= Map->TransientEndTime+Record->NumSnapshot || Map->t > Map->RecordEndTime-Record->NumSnapshot))) Record->RecordSnapshot(Map);
		if (Record->isRecordingSnapshotGS && (!Record->SnapshotFirstLastSw || (Map->t <= Map->TransientEndTime+Record->NumSnapshot || Map->t > Map->RecordEndTime-Record->NumSnapshot))) Record->RecordSnapshotGS(Map);
		if (Record->isRecordingDensityTimeSeries) Record->RecordDensityTimeSeries(Map);
		if (Record->isRecordingOutputQuantity && Map->ExistOutputQuantity && (!Record->OutputQuantityFirstLastSw || (Map->t <= Map->TransientEndTime+Record->NumOutputQuantity || Map->t > Map->RecordEndTime-Record->NumOutputQuantity))) {Map->CalcOutputQuantity(); Record->RecordOutputQuantity(Map);}
		if (Map->isRecordingQuantityHist && !(Map->isQuantityHistBack)) Map->AccQuantityHist();
		//Store data;
		Map->StoreXdX();
		Map->StoreMatR();
		//Display progress bar;
		Record->ProgBar();
	}
	//Final transient;
	Map->isNowBackward = true;
	Map->NextGSTime = Map->t + Map->GSInterval;
	while (Map->t < Map->StoreEndTime) {
		//Time evolution;
		Map->f();
		//Gram-Schmidt orthonormalization;
		if (Map->t == Map->NextGSTime) {
			Map->GramSchmidt();
			Map->NextGSTime += Map->GSInterval;
			//StoreData;
			Map->StoreMatR();
		}
		//Display progress bar;
		Record->ProgBar();
	}
	//Switch to backward evolution;
	Map->SetInitialMatC();
	//Backward transient;
	while (Map->t > Map->RecordEndTime) {
		Map->BackwardEvolution();
		Map->t -= Map->GSInterval;
	}
	//Backward evolution;
	while (Map->t > Map->TransientEndTime) {
		//Calculate covariant vectors;
		Record->SetRecordTimeSeriesSw(Map);
		Map->RestoreXdX();
		Map->CovariantVectors();
		Map->CalcPtcpRatio(Map->CV, Map->CVPtcpRatio);
		if (Record->isRecordingSnapshotCV && (!Record->SnapshotFirstLastSw || (Map->t <= Map->TransientEndTime+Record->NumSnapshot || Map->t > Map->RecordEndTime-Record->NumSnapshot))) Record->RecordSnapshotCV(Map);
		if (Record->isRecordingSmallestAngle || Record->isRecordingAllAngles) Record->AccAngleHist(Map);
		if (Record->isRecordingNeighAngle) Record->AccNeighAngleHist(Map);
		if (Record->isRecordingSubspaceAngle) Record->AccSubspaceAngleHist(Map);
		if (Record->isRecordingAngleTimeSeries) Record->RecordAngleTimeSeries(Map);
		if (Map->isRecordingPowerSpectra) Map->AccPowerSpectra();
		if (Record->isRecordingOutputQuantity && Map->ExistOutputQuantityBack && (!Record->OutputQuantityFirstLastSw || (Map->t <= Map->TransientEndTime+Record->NumOutputQuantity || Map->t > Map->RecordEndTime-Record->NumOutputQuantity))) {Map->CalcOutputQuantityBack(); Record->RecordOutputQuantityBack(Map);}
		if (Record->isRecordingSupplPtcpRatio) Record->CalcSupplPtcpRatio(Map);
		if (Map->isRecordingQuantityHist && Map->isQuantityHistBack) Map->AccQuantityHist();
		//Check exponents for covariant vectors (Caution! It changes X and CV);
		Map->CheckCVExp();
		Record->SetRecordTimeSeriesSw(Map, -1);
		if (Record->isRecordingCVExp) Record->RecordCVExp(Map);
		if (Record->isRecordingDOSViolation) Record->CountDOSViolation(Map);
		if (Record->isRecordingPtcpHist) Record->AccPtcpHist(Map);
		if (Map->isCheckingCVExp && Record->isRecordingLyapHist && Record->LyapHistGSCV) Record->AccLyapHist(Map);
		//Backward evolution;
		Map->BackwardEvolution();
		Map->t--;
		//Display progress bar;
		Record->ProgBar();
	}
}

//Main routine with loop division;
void MainRoutineWithDivision(map *Map, record *Record) {
	//Thermalization;
	//Pretend to have NumLyap=0;
	int NumLyapBak = Map->NumLyap;
	Map->NumLyap = 0;
	while (Map->t < Map->ThermalizationTime) {
		//Time evolution;
		Map->f();
		//Display progress bar;
		Record->ProgBar();
	}
	//Initial transient;
	Map->NumLyap = NumLyapBak;
	Map->NextGSTime = Map->ThermalizationTime + Map->GSInterval;
	while (Map->t < Map->TransientEndTime) {
		//Time evolution;
		Map->f();
		//Gram-Schmidt orthonormalization;
		if (Map->t == Map->NextGSTime) {
			Map->GramSchmidt();
			Map->NextGSTime += Map->GSInterval;
		}
		//Display progress bar;
		Record->ProgBar();
	}
	//Time evolution, just store variables at junctions;
	while (Map->t < Map->RecordEndTime) {
		//Write variables;
		if (Map->t == Map->LoopTimes[Map->LoopInd]) Map->WriteVar();
		//Time evolution;
		Map->f();
		Record->SetRecordTimeSeriesSw(Map);
		//Time series record;
		if (Record->isRecordingTimeSeries || Map->isRecordingTsCumulants) Map->CalcMeanField();
		if (Record->isRecordingTimeSeries) Record->RecordTimeSeries(Map);
		if (Map->isRecordingTsCumulants) Map->CalcTsCumulants();
		if (Record->isRecordingSnapHist) Record->AccSnapHist(Map);
		if (Map->isRecordingPowSpecDyn) Map->AccPowSpecDyn();
		//Gram-Schmidt orthonormalization;
		Map->GramSchmidt();
		//Record it;
		Map->CalcLyap();
		Map->CalcPtcpRatio(Map->dX, Map->GSPtcpRatio);
		if (Record->isRecordingLyap && (Record->isRecordingLyapTimeSeriesAlways || Map->t == (Map->LoopTimes)[Map->LoopInd]))
			Record->RecordLyap(Map);
		if (Record->isRecordingLyapHist && !Record->LyapHistGSCV) Record->AccLyapHist(Map);
		if (Record->isRecordingSnapshot && (!Record->SnapshotFirstLastSw || (Map->t <= Map->TransientEndTime+Record->NumSnapshot || Map->t > Map->RecordEndTime-Record->NumSnapshot))) Record->RecordSnapshot(Map);
		if (Record->isRecordingSnapshotGS && (!Record->SnapshotFirstLastSw || (Map->t <= Map->TransientEndTime+Record->NumSnapshot || Map->t > Map->RecordEndTime-Record->NumSnapshot))) Record->RecordSnapshotGS(Map);
		if (Record->isRecordingDensityTimeSeries) Record->RecordDensityTimeSeries(Map);
		if (Record->isRecordingOutputQuantity && Map->ExistOutputQuantity && (!Record->OutputQuantityFirstLastSw || (Map->t <= Map->TransientEndTime+Record->NumOutputQuantity || Map->t > Map->RecordEndTime-Record->NumOutputQuantity))) {Map->CalcOutputQuantity(); Record->RecordOutputQuantity(Map);}
		if (Map->isRecordingQuantityHist && !(Map->isQuantityHistBack)) Map->AccQuantityHist();
		//Display progress bar;
		Record->ProgBar();
	}
	//Final transient, just store variables at junctions;
	Map->isNowBackward = true;
	Map->NextGSTime = Map->t + Map->GSInterval;
	while (Map->t < Map->StoreEndTime) {
		//Write variables;
		if (Map->t == Map->LoopTimes[Map->LoopInd]) Map->WriteVar();
		//Time evolution;
		Map->f();
		//Gram-Schmidt orthonormalization;
		if (Map->t == Map->NextGSTime) {
			Map->GramSchmidt();
			Map->NextGSTime += Map->GSInterval;
		}
		//Display progress bar;
		Record->ProgBar();
	}
	//Switch to backward evolution;
	Map->SetInitialMatC();
	//Backward transient, repeat forward-backward loops;
	while (Map->LoopInd > Map->NumRecLoop) {
		//Read variables;
		Map->ReadVar();
		Map->NextGSTime = Map->t + Map->GSInterval;
		//Move forward;
		while (Map->t < (Map->LoopTimes)[Map->LoopInd+1]) {
			//Time evolution;
			Map->f();
			//Gram-Schmidt orthonormalization;
			if (Map->t == Map->NextGSTime) {
				Map->GramSchmidt();
				Map->NextGSTime += Map->GSInterval;
				//StoreData;
				Map->StoreMatR();
			}
			//Display progress bar;
			Record->ProgBar();
		}
		//Move backward;
		while (Map->t > (Map->LoopTimes)[Map->LoopInd]) {
			Map->BackwardEvolution();
			Map->t -= Map->GSInterval;
		}
	}
	//Save data here;
	Record->WriteAfterBackwardTransient(Map);
	//Backward evolution, repeat forward-backward loops;
	while (Map->LoopInd > 0) {
		//Read variables;
		Map->ReadVar();
		//Move forward;
		while (Map->t < (Map->LoopTimes)[Map->LoopInd+1]) {
			//Time evolution;
			Map->f();
			//Gram-Schmidt orthonormalization;
			Map->GramSchmidt();
			//Store data;
			Map->StoreXdX();
			Map->StoreMatR();
			//Display progress bar;
			Record->ProgBar();
		}
		//Move backward;
		while (Map->t > (Map->LoopTimes)[Map->LoopInd]) {
			//Calculate covariant vectors;
			Record->SetRecordTimeSeriesSw(Map);
			Map->RestoreXdX();
			Map->CovariantVectors();
			Map->CalcPtcpRatio(Map->CV, Map->CVPtcpRatio);
			if (Record->isRecordingSnapshotCV && (!Record->SnapshotFirstLastSw || (Map->t <= Map->TransientEndTime+Record->NumSnapshot || Map->t > Map->RecordEndTime-Record->NumSnapshot))) Record->RecordSnapshotCV(Map);
			if (Record->isRecordingSmallestAngle || Record->isRecordingAllAngles) Record->AccAngleHist(Map);
			if (Record->isRecordingNeighAngle) Record->AccNeighAngleHist(Map);
			if (Record->isRecordingSubspaceAngle) Record->AccSubspaceAngleHist(Map);
			if (Record->isRecordingAngleTimeSeries) Record->RecordAngleTimeSeries(Map);
			if (Map->isRecordingPowerSpectra) Map->AccPowerSpectra();
			if (Record->isRecordingOutputQuantity && Map->ExistOutputQuantityBack && (!Record->OutputQuantityFirstLastSw || (Map->t <= Map->TransientEndTime+Record->NumOutputQuantity || Map->t > Map->RecordEndTime-Record->NumOutputQuantity))) {Map->CalcOutputQuantityBack(); Record->RecordOutputQuantityBack(Map);}
			if (Record->isRecordingSupplPtcpRatio) Record->CalcSupplPtcpRatio(Map);
			if (Map->isRecordingQuantityHist && Map->isQuantityHistBack) Map->AccQuantityHist();
			//Check exponents for covariant vectors (Caution! It changes X and CV);
			Map->CheckCVExp();
			Record->SetRecordTimeSeriesSw(Map, -1);
			if (Record->isRecordingCVExp && (Record->isRecordingLyapTimeSeriesAlways || Map->t - 1 == (Map->LoopTimes)[Map->LoopInd]))
				Record->RecordCVExp(Map);
			if (Record->isRecordingDOSViolation) Record->CountDOSViolation(Map);
			if (Record->isRecordingPtcpHist) Record->AccPtcpHist(Map);
			if (Map->isCheckingCVExp && Record->isRecordingLyapHist && Record->LyapHistGSCV) Record->AccLyapHist(Map);
			//Backward evolution;
			Map->BackwardEvolution();
			Map->t--;
			//Display progress bar;
			Record->ProgBar();
		}
	}
}

//Main routine of inverse dynamics without loop division
void MainRoutineInverseWithoutDivision(map *Map, record *Record) {
	//Pretend to have NumLyap=0;
	int NumLyapBak = Map->NumLyap;
	Map->NumLyap = 0;
	//Initial transient;
	while (Map->t < Map->TransientEndTime) {
		//Time evolution;
		Map->f();
	}
	//Time evolution;
	Map->StoreX();
	while (Map->t < Map->RecordEndTime) {
		//Time evolution;
		Map->f();
		Record->SetRecordTimeSeriesSw(Map);
		//Time series record;
		if (Record->isRecordingTimeSeries || Map->isRecordingTsCumulants) Map->CalcMeanField();
		if (Record->isRecordingTimeSeries) Record->RecordTimeSeries(Map);
		if (Map->isRecordingTsCumulants) Map->CalcTsCumulants();
		if (Record->isRecordingSnapHist) Record->AccSnapHist(Map);
		if (Map->isRecordingPowSpecDyn) Map->AccPowSpecDyn();
		//Store data;
		Map->StoreX();
	}
	//Final transient;
	Map->isNowBackward = true;
	while (Map->t < Map->StoreEndTime) {
		//Time evolution;
		Map->f();
		//Store data;
		Map->StoreX();
	}
	//Switch to backward evolution;
	Map->NumLyap = NumLyapBak;
	Map->NextGSTime = Map->StoreEndTime - Map->GSInterval;
	Map->RestoreX();
	//Backward transient;
	while (Map->t > Map->RecordEndTime) {
		//Backward evolution;
		Map->BackDF(Map->dX);
		//Gram-Schmidt orthonormalization;
		if (Map->t == Map->NextGSTime) {
			Map->GramSchmidt();
			Map->NextGSTime -= Map->GSInterval;
		}
		//Display progress bar;
		Record->ProgBar();
	}
	//Backward evolution;
	while (Map->t > Map->TransientEndTime) {
		//Backward evolution;
		Map->BackDF(Map->dX);
		Record->SetRecordTimeSeriesSw(Map);
		//Gram-Schmidt orthonormalization;
		Map->GramSchmidt();
		//Record it;
		Map->CalcLyap();
		Map->CalcPtcpRatio(Map->dX, Map->GSPtcpRatio);
		if (Record->isRecordingLyap) Record->RecordLyap(Map);
		if (Record->isRecordingLyapHist && !Record->LyapHistGSCV) Record->AccLyapHist(Map);
		if (Record->isRecordingSnapshot && (!Record->SnapshotFirstLastSw || (Map->t < Map->TransientEndTime+Record->NumSnapshot || Map->t >= Map->RecordEndTime-Record->NumSnapshot))) Record->RecordSnapshot(Map);
		if (Record->isRecordingSnapshotGS && (!Record->SnapshotFirstLastSw || (Map->t < Map->TransientEndTime+Record->NumSnapshot || Map->t >= Map->RecordEndTime-Record->NumSnapshot))) Record->RecordSnapshotGS(Map);
		if (Record->isRecordingDensityTimeSeries) Record->RecordDensityTimeSeries(Map);
		if (Record->isRecordingOutputQuantity && Map->ExistOutputQuantity && (!Record->OutputQuantityFirstLastSw || (Map->t <= Map->TransientEndTime+Record->NumOutputQuantity || Map->t > Map->RecordEndTime-Record->NumOutputQuantity))) {Map->CalcOutputQuantity(); Record->RecordOutputQuantity(Map);}
		if (Map->isRecordingQuantityHist && !(Map->isQuantityHistBack)) Map->AccQuantityHist();
		//Display progress bar;
		Record->ProgBar();
	}
}

//Main routine of inverse dynamics with loop division;
void MainRoutineInverseWithDivision(map *Map, record *Record) {
	//Pretend to have NumLyap=0;
	int NumLyapBak = Map->NumLyap;
	Map->NumLyap = 0;
	//Initial transient;
	while (Map->t < Map->TransientEndTime) {
		//Time evolution;
		Map->f();
	}
	//Time evolution, just store variables at junctions;
	while (Map->t < Map->RecordEndTime) {
		//Write variables;
		if (Map->t == Map->LoopTimes[Map->LoopInd]) Map->WriteVarX();
		//Time evolution;
		Map->f();
		Record->SetRecordTimeSeriesSw(Map);
		//Time series record;
		if (Record->isRecordingTimeSeries || Map->isRecordingTsCumulants) Map->CalcMeanField();
		if (Record->isRecordingTimeSeries) Record->RecordTimeSeries(Map);
		if (Map->isRecordingTsCumulants) Map->CalcTsCumulants();
		if (Record->isRecordingSnapHist) Record->AccSnapHist(Map);
		if (Map->isRecordingPowSpecDyn) Map->AccPowSpecDyn();
	}
	//Final transient, just store variables at junctions;
	Map->isNowBackward = true;
	while (Map->t < Map->StoreEndTime) {
		//Write variables;
		if (Map->t == Map->LoopTimes[Map->LoopInd]) Map->WriteVarX();
		//Time evolution;
		Map->f();
	}
	//Switch to backward evolution;
	//Backward transient, repeat forward-backward loops;
	while (Map->LoopInd > Map->NumRecLoop) {
		//Read variables;
		Map->ReadVarX();
		//Move forward;
		Map->NumLyap = 0;
		Map->StoreX();
		while (Map->t < (Map->LoopTimes)[Map->LoopInd+1]) {
			//Time evolution;
			Map->f();
			//Store data;
			Map->StoreX();
		}
		//Move backward;
		Map->NumLyap = NumLyapBak;
		Map->NextGSTime = Map->t - Map->GSInterval;
		Map->RestoreX();
		while (Map->t > (Map->LoopTimes)[Map->LoopInd]) {
			//Backward evolution;
			Map->BackDF(Map->dX);
			//Gram-Schmidt orthonormalization;
			if (Map->t == Map->NextGSTime) {
				Map->GramSchmidt();
				Map->NextGSTime -= Map->GSInterval;
			}
			//Display progress bar;
			Record->ProgBar();
		}
	}
	//Backward evolution, repeat forward-backward loops;
	while (Map->LoopInd > 0) {
		//Read variables;
		Map->ReadVarX();
		//Move forward;
		Map->NumLyap = 0;
		Map->StoreX();
		while (Map->t < (Map->LoopTimes)[Map->LoopInd+1]) {
			//Time evolution;
			Map->f();
			//Store data;
			Map->StoreX();
		}
		//Move backward;
		Map->NumLyap = NumLyapBak;
		Map->RestoreX();
		while (Map->t > (Map->LoopTimes)[Map->LoopInd]) {
			//Backward evolution;
			Map->BackDF(Map->dX);
			Record->SetRecordTimeSeriesSw(Map);
			//Gram-Schmidt orthonormalization;
			Map->GramSchmidt();
			//Record it;
			Map->CalcLyap();
			Map->CalcPtcpRatio(Map->dX, Map->GSPtcpRatio);
			if (Record->isRecordingLyap && (Record->isRecordingLyapTimeSeriesAlways || Map->t == (Map->LoopTimes)[Map->LoopInd]))
				Record->RecordLyap(Map);
			if (Record->isRecordingLyapHist && !Record->LyapHistGSCV) Record->AccLyapHist(Map);
			if (Record->isRecordingSnapshot && (!Record->SnapshotFirstLastSw || (Map->t < Map->TransientEndTime+Record->NumSnapshot || Map->t >= Map->RecordEndTime-Record->NumSnapshot))) Record->RecordSnapshot(Map);
			if (Record->isRecordingSnapshotGS && (!Record->SnapshotFirstLastSw || (Map->t < Map->TransientEndTime+Record->NumSnapshot || Map->t >= Map->RecordEndTime-Record->NumSnapshot))) Record->RecordSnapshotGS(Map);
			if (Record->isRecordingDensityTimeSeries) Record->RecordDensityTimeSeries(Map);
			if (Record->isRecordingOutputQuantity && Map->ExistOutputQuantity && (!Record->OutputQuantityFirstLastSw || (Map->t <= Map->TransientEndTime+Record->NumOutputQuantity || Map->t > Map->RecordEndTime-Record->NumOutputQuantity))) {Map->CalcOutputQuantity(); Record->RecordOutputQuantity(Map);}
			if (Map->isRecordingQuantityHist && !(Map->isQuantityHistBack)) Map->AccQuantityHist();
			//Display progress bar;
			Record->ProgBar();
		}
	}
}

//Main routine of inverse dynamics & forward CV without loop division
void MainRoutineInverseCVWithoutDivision(map *Map, record *Record) {
	int i;
	//Pretend to have NumLyap=0;
	int NumLyapBak = Map->NumLyap;
	Map->NumLyap = 0;
	//Thermalization;
	while (Map->t < Map->ThermalizationTime) {
		//Time evolution;
		Map->f();
	}
	//Initial transient;
	Map->StoreX();
	while (Map->t < Map->TransientEndTime) {
		//Time evolution;
		Map->f();
		//Store data;
		Map->StoreX();
	}
	//Time evolution;
	while (Map->t < Map->RecordEndTime) {
		//Time evolution;
		Map->f();
		Record->SetRecordTimeSeriesSw(Map);
		//Time series record;
		if (Record->isRecordingTimeSeries || Map->isRecordingTsCumulants) Map->CalcMeanField();
		if (Record->isRecordingTimeSeries) Record->RecordTimeSeries(Map);
		if (Map->isRecordingTsCumulants) Map->CalcTsCumulants();
		if (Record->isRecordingSnapHist) Record->AccSnapHist(Map);
		if (Map->isRecordingPowSpecDyn) Map->AccPowSpecDyn();
		//Store data;
		Map->StoreX();
	}
	//Final transient;
	Map->isNowBackward = true;
	while (Map->t < Map->StoreEndTime) {
		//Time evolution;
		Map->f();
		//Store data;
		Map->StoreX();
	}
	//Switch to backward evolution;
	Map->NumLyap = NumLyapBak;
	Map->RestoreX();
	Map->NextGSTime = Map->StoreEndTime - Map->GSInterval;
	//Backward transient;
	while (Map->t > Map->RecordEndTime) {
		//Backward evolution;
		Map->BackDF(Map->dX);
		//Gram-Schmidt orthonormalization;
		if (Map->t == Map->NextGSTime) {
			Map->GramSchmidt();
			Map->NextGSTime -= Map->GSInterval;
		}
		//Display progress bar;
		Record->ProgBar();
	}
	//Backward evolution;
	while (Map->t > Map->TransientEndTime) {
		//Backward evolution;
		Map->BackDF(Map->dX);
		Record->SetRecordTimeSeriesSw(Map);
		//Gram-Schmidt orthonormalization;
		Map->GramSchmidt();
		//Record it;
		Map->CalcLyap();
		Map->CalcPtcpRatio(Map->dX, Map->GSPtcpRatio);
		if (Record->isRecordingLyap) Record->RecordLyap(Map);
		if (Record->isRecordingLyapHist && !Record->LyapHistGSCV) Record->AccLyapHist(Map);
		if (Record->isRecordingSnapshot && (!Record->SnapshotFirstLastSw || (Map->t < Map->TransientEndTime+Record->NumSnapshot || Map->t >= Map->RecordEndTime-Record->NumSnapshot))) Record->RecordSnapshot(Map);
		if (Record->isRecordingSnapshotGS && (!Record->SnapshotFirstLastSw || (Map->t < Map->TransientEndTime+Record->NumSnapshot || Map->t >= Map->RecordEndTime-Record->NumSnapshot))) Record->RecordSnapshotGS(Map);
		if (Record->isRecordingDensityTimeSeries) Record->RecordDensityTimeSeries(Map);
		if (Record->isRecordingOutputQuantity && Map->ExistOutputQuantity && (!Record->OutputQuantityFirstLastSw || (Map->t <= Map->TransientEndTime+Record->NumOutputQuantity || Map->t > Map->RecordEndTime-Record->NumOutputQuantity))) {Map->CalcOutputQuantity(); Record->RecordOutputQuantity(Map);}
		if (Map->isRecordingQuantityHist && !(Map->isQuantityHistBack)) Map->AccQuantityHist();
		//Store data;
		Map->StoreDX();
		Map->StoreMatR();
		//Display progress bar;
		Record->ProgBar();
	}
	//Backward final transient;
	Map->NextGSTime = Map->TransientEndTime - Map->GSInterval;
	while (Map->t > Map->ThermalizationTime) {
		//Backward evolution;
		Map->BackDF(Map->dX);
		//Gram-Schmidt orthonormalization;
		if (Map->t == Map->NextGSTime) {
			Map->GramSchmidt();
			Map->NextGSTime -= Map->GSInterval;
			//StoreData;
			Map->StoreMatR();
		}
		//Display progress bar;
		Record->ProgBar();
	}
	//Switch to 2nd forward evolution;
	Map->SetInitialMatC();
	//2nd forward transient;
	while (Map->t < Map->TransientEndTime) {
		Map->BackwardEvolution();
		Map->t += Map->GSInterval;
		for (i=0;i<Map->GSInterval;i++) Map->ReRestoreX();
	}
	//2nd forward evolution;
	while (Map->t < Map->RecordEndTime) {
		//Calculate covariant vectors;
		Record->SetRecordTimeSeriesSw(Map);
		Map->RestoreDX();
		Map->CovariantVectors();
		Map->CalcPtcpRatio(Map->CV, Map->CVPtcpRatio);
		if (Record->isRecordingSnapshotCV && (!Record->SnapshotFirstLastSw || (Map->t < Map->TransientEndTime+Record->NumSnapshot || Map->t >= Map->RecordEndTime-Record->NumSnapshot))) Record->RecordSnapshotCV(Map);
		if (Record->isRecordingSmallestAngle || Record->isRecordingAllAngles) Record->AccAngleHist(Map);
		if (Record->isRecordingNeighAngle) Record->AccNeighAngleHist(Map);
		if (Record->isRecordingSubspaceAngle) Record->AccSubspaceAngleHist(Map);
		if (Record->isRecordingAngleTimeSeries) Record->RecordAngleTimeSeries(Map);
		if (Map->isRecordingPowerSpectra) Map->AccPowerSpectra();
		if (Record->isRecordingOutputQuantity && Map->ExistOutputQuantityBack && (!Record->OutputQuantityFirstLastSw || (Map->t <= Map->TransientEndTime+Record->NumOutputQuantity || Map->t > Map->RecordEndTime-Record->NumOutputQuantity))) {Map->CalcOutputQuantityBack(); Record->RecordOutputQuantityBack(Map);}
		if (Record->isRecordingSupplPtcpRatio) Record->CalcSupplPtcpRatio(Map);
		if (Map->isRecordingQuantityHist && Map->isQuantityHistBack) Map->AccQuantityHist();
		//Check exponents for covariant vectors (Caution! It changes X and CV);
		Map->CheckBackCVExp();
		Record->SetRecordTimeSeriesSw(Map, +1);
		if (Record->isRecordingCVExp) Record->RecordBackCVExp(Map);
		if (Record->isRecordingDOSViolation) Record->CountDOSViolation(Map);
		if (Record->isRecordingPtcpHist) Record->AccPtcpHist(Map);
		if (Map->isCheckingCVExp && Record->isRecordingLyapHist && Record->LyapHistGSCV) Record->AccLyapHist(Map);
		//Backward evolution;
		Map->BackwardEvolution();
		Map->ReRestoreX();
		Map->t++;
		//Display progress bar;
		Record->ProgBar();
	}
}

//Main routine of inverse dynamics & forward CV with loop division;
void MainRoutineInverseCVWithDivision(map *Map, record *Record) {
	int i;
	//Pretend to have NumLyap=0;
	int NumLyapBak = Map->NumLyap;
	Map->NumLyap = 0;
	//Thermalization;
	while (Map->t < Map->ThermalizationTime) {
		//Time evolution;
		Map->f();
	}
	//Initial transient;
	while (Map->t < Map->TransientEndTime) {
		//Write variables;
		if (Map->t == Map->LoopTimes[Map->LoopInd] || Map->t == Map->LoopTimes[Map->NumTransLoop]-1) Map->WriteVar();
		//Time evolution;
		Map->f();
	}
	//Time evolution, just store variables at junctions;
	while (Map->t < Map->RecordEndTime) {
		//Write variables;
		if (Map->LoopInd < Map->NumRecLoop && Map->t == Map->LoopTimes[Map->LoopInd]-1) Map->WriteVar();
		//Time evolution;
		Map->f();
		Record->SetRecordTimeSeriesSw(Map);
		//Time series record;
		if (Record->isRecordingTimeSeries || Map->isRecordingTsCumulants) Map->CalcMeanField();
		if (Record->isRecordingTimeSeries) Record->RecordTimeSeries(Map);
		if (Map->isRecordingTsCumulants) Map->CalcTsCumulants();
		if (Record->isRecordingSnapHist) Record->AccSnapHist(Map);
		if (Map->isRecordingPowSpecDyn) Map->AccPowSpecDyn();
	}
	//Final transient, just store variables at junctions;
	Map->isNowBackward = true;
	while (Map->t < Map->StoreEndTime) {
		//Write variables;
		if (Map->t == Map->LoopTimes[Map->LoopInd])	Map->WriteVar();
		//Time evolution;
		Map->f();
	}
	//Switch to backward evolution;
	//Backward transient, repeat forward-backward loops;
	while (Map->LoopInd > Map->NumRecLoop) {
		//Read variables;
		Map->ReadVarNoDX();
		//Move forward;
		Map->NumLyap = 0;
		Map->StoreX();
		while (Map->t < (Map->LoopTimes)[Map->LoopInd+1]) {
			//Time evolution;
			Map->f();
			//Store data;
			Map->StoreX();
		}
		//Move backward;
		Map->NumLyap = NumLyapBak;
		Map->RestoreX();
		Map->NextGSTime = Map->t - Map->GSInterval;
		while (Map->t > (Map->LoopTimes)[Map->LoopInd]) {
			//Backward evolution;
			Map->BackDF(Map->dX);
			//Gram-Schmidt orthonormalization;
			if (Map->t == Map->NextGSTime) {
				Map->GramSchmidt();
				Map->NextGSTime -= Map->GSInterval;
			}
			//Display progress bar;
			Record->ProgBar();
		}
	}
	//Backward evolution, repeat forward-backward loops;
	while (Map->LoopInd > Map->NumTransLoop) {
		//Read variables;
		Map->ReadVarNoDX();
		Map->t--;	//Stored data are at time LoopTimes-1;
		//Write variables;
		Map->AppendVarDX();
		//Move forward;
		Map->NumLyap = 0;
		Map->StoreX();
		while (Map->t < (Map->LoopTimes)[Map->LoopInd+1]) {
			//Time evolution;
			Map->f();
			//Store data;
			Map->StoreX();
		}
		//Move backward;
		Map->NumLyap = NumLyapBak;
		Map->RestoreX();
		while (Map->t > (Map->LoopTimes)[Map->LoopInd]) {
			//Backward evolution;
			Map->BackDF(Map->dX);
			Record->SetRecordTimeSeriesSw(Map);
			//Gram-Schmidt orthonormalization;
			Map->GramSchmidt();
			//Record it;
			Map->CalcLyap();
			Map->CalcPtcpRatio(Map->dX, Map->GSPtcpRatio);
			if (Record->isRecordingLyap && (Record->isRecordingLyapTimeSeriesAlways || Map->t == (Map->LoopTimes)[Map->LoopInd]))
				Record->RecordLyap(Map);
			if (Record->isRecordingLyapHist && !Record->LyapHistGSCV) Record->AccLyapHist(Map);
			if (Record->isRecordingSnapshot && (!Record->SnapshotFirstLastSw || (Map->t < Map->TransientEndTime+Record->NumSnapshot || Map->t >= Map->RecordEndTime-Record->NumSnapshot))) Record->RecordSnapshot(Map);
			if (Record->isRecordingSnapshotGS && (!Record->SnapshotFirstLastSw || (Map->t < Map->TransientEndTime+Record->NumSnapshot || Map->t >= Map->RecordEndTime-Record->NumSnapshot))) Record->RecordSnapshotGS(Map);
			if (Record->isRecordingDensityTimeSeries) Record->RecordDensityTimeSeries(Map);
			if (Record->isRecordingOutputQuantity && Map->ExistOutputQuantity && (!Record->OutputQuantityFirstLastSw || (Map->t <= Map->TransientEndTime+Record->NumOutputQuantity || Map->t > Map->RecordEndTime-Record->NumOutputQuantity))) {Map->CalcOutputQuantity(); Record->RecordOutputQuantity(Map);}
			if (Map->isRecordingQuantityHist && !(Map->isQuantityHistBack)) Map->AccQuantityHist();
			//Display progress bar;
			Record->ProgBar();
		}
		//Discard excessively stored data;
		Map->RestoreX();
	}
	//Backward final transient, repeat forward-backward loops;
	while (Map->LoopInd > 0) {
		//Read variables;
		Map->ReadVarNoDX();
		//Write variables;
		Map->AppendVarDX();
		//Move forward;
		Map->NumLyap = 0;
		Map->StoreX();
		while (Map->t < (Map->LoopTimes)[Map->LoopInd+1]) {
			//Time evolution;
			Map->f();
			//Store data;
			Map->StoreX();
		}
		//Move backward;
		Map->NumLyap = NumLyapBak;
		Map->RestoreX();
		Map->NextGSTime = Map->t - Map->GSInterval;
		while (Map->t > (Map->LoopTimes)[Map->LoopInd]) {
			//Backward evolution;
			Map->BackDF(Map->dX);
			//Gram-Schmidt orthonormalization;
			if (Map->t == Map->NextGSTime) {
				Map->GramSchmidt();
				Map->NextGSTime -= Map->GSInterval;
			}
			//Display progress bar;
			Record->ProgBar();
		}
	}
	//Switch to 2nd forward evolution;
	Map->SetInitialMatC();
	//2nd forward transient, repeat forward-backward-forward loops;
	while (Map->LoopInd < Map->NumTransLoop) {
		//Read variables;
		Map->LoopInd++;
		Map->ReadVar();
		//Move forward;
		Map->NumLyap = 0;
		Map->StoreX();
		while (Map->t < (Map->LoopTimes)[Map->LoopInd+1]) {
			//Time evolution;
			Map->f();
			//Store data;
			Map->StoreX();
		}
		//Move backward;
		Map->NumLyap = NumLyapBak;
		Map->RestoreX();
		Map->NextGSTime = Map->t - Map->GSInterval;
		while (Map->t > (Map->LoopTimes)[Map->LoopInd]) {
			//Backward evolution;
			Map->BackDF(Map->dX);
			//Gram-Schmidt orthonormalization;
			if (Map->t == Map->NextGSTime) {
				Map->GramSchmidt();
				Map->NextGSTime -= Map->GSInterval;
				//StoreData;
				Map->StoreMatR();
			}
			//Display progress bar;
			Record->ProgBar();
		}
		//Move again forward;
		while (Map->t < (Map->LoopTimes)[Map->LoopInd+1]) {
			Map->BackwardEvolution();
			Map->t += Map->GSInterval;
		}
		//Increment LoopInd;
		Map->LoopInd++;
	}
	//2nd forward evolution, repeat forward-backward-forward loops;
	while (Map->LoopInd < Map->NumRecLoop) {
		//Read variables;
		Map->LoopInd++;
		Map->ReadVar();
		Map->t--;	//Stored data are at time LoopTimes-1;
		//Move forward;
		Map->NumLyap = 0;
		Map->StoreX();
		while (Map->t < (Map->LoopTimes)[Map->LoopInd+1]) {
			//Time evolution;
			Map->f();
			//Store data;
			Map->StoreX();
		}
		//Move backward;
		Map->NumLyap = NumLyapBak;
		Map->RestoreX();
		while (Map->t > (Map->LoopTimes)[Map->LoopInd]) {
			//Backward evolution;
			Map->BackDF(Map->dX);
			//Gram-Schmidt orthonormalization;
			Map->GramSchmidt();
			//Store data;
			Map->StoreDX();
			Map->StoreMatR();
			//Display progress bar;
			Record->ProgBar();
		}
		//Move again forward;
		while (Map->t < (Map->LoopTimes)[Map->LoopInd+1]) {
			//Calculate covariant vectors;
			Record->SetRecordTimeSeriesSw(Map);
			Map->RestoreDX();
			Map->CovariantVectors();
			Map->CalcPtcpRatio(Map->CV, Map->CVPtcpRatio);
			if (Record->isRecordingSnapshotCV && (!Record->SnapshotFirstLastSw || (Map->t < Map->TransientEndTime+Record->NumSnapshot || Map->t >= Map->RecordEndTime-Record->NumSnapshot))) Record->RecordSnapshotCV(Map);
			if (Record->isRecordingSmallestAngle || Record->isRecordingAllAngles) Record->AccAngleHist(Map);
			if (Record->isRecordingNeighAngle) Record->AccNeighAngleHist(Map);
			if (Record->isRecordingSubspaceAngle) Record->AccSubspaceAngleHist(Map);
			if (Record->isRecordingAngleTimeSeries) Record->RecordAngleTimeSeries(Map);
			if (Map->isRecordingPowerSpectra) Map->AccPowerSpectra();
			if (Record->isRecordingOutputQuantity && Map->ExistOutputQuantityBack && (!Record->OutputQuantityFirstLastSw || (Map->t <= Map->TransientEndTime+Record->NumOutputQuantity || Map->t > Map->RecordEndTime-Record->NumOutputQuantity))) {Map->CalcOutputQuantityBack(); Record->RecordOutputQuantityBack(Map);}
			if (Record->isRecordingSupplPtcpRatio) Record->CalcSupplPtcpRatio(Map);
			if (Map->isRecordingQuantityHist && Map->isQuantityHistBack) Map->AccQuantityHist();
			//Check exponents for covariant vectors (Caution! It changes X and CV);
			Map->CheckBackCVExp();
			Record->SetRecordTimeSeriesSw(Map, +1);
			if (Record->isRecordingCVExp && (Record->isRecordingLyapTimeSeriesAlways || Map->t + 1 == (Map->LoopTimes)[Map->LoopInd+1]))
				Record->RecordBackCVExp(Map);
			if (Record->isRecordingDOSViolation) Record->CountDOSViolation(Map);
			if (Record->isRecordingPtcpHist) Record->AccPtcpHist(Map);
			if (Map->isCheckingCVExp && Record->isRecordingLyapHist && Record->LyapHistGSCV) Record->AccLyapHist(Map);
			//Backward evolution;
			Map->BackwardEvolution();
			Map->ReRestoreX();
			Map->t++;
			//Display progress bar;
			Record->ProgBar();
		}
		//Increment LoopInd & reset pStoredX;
		Map->LoopInd++;
		Map->ResetPStored();
	}
}

//Main routine with loop division from saved data;
void MainRoutineWithDivisionAfterTransient(map *Map, record *Record) {
	//Setup
	Map->isNowBackward = true;
	//Load data of MatC;
	Record->ReadAfterBackwardTransient(Map);
	//Backward evolution, repeat forward-backward loops;
	while (Map->LoopInd > 0) {
		//Read variables;
		Map->ReadVar();
		//Move forward;
		while (Map->t < (Map->LoopTimes)[Map->LoopInd+1]) {
			//Time evolution;
			Map->f();
			//Gram-Schmidt orthonormalization;
			Map->GramSchmidt();
			//Store data;
			Map->StoreXdX();
			Map->StoreMatR();
			//Display progress bar;
			Record->ProgBar();
		}
		//Move backward;
		while (Map->t > (Map->LoopTimes)[Map->LoopInd]) {
			//Calculate covariant vectors;
			Record->SetRecordTimeSeriesSw(Map);
			Map->RestoreXdX();
			Map->CovariantVectors();
			Map->CalcPtcpRatio(Map->CV, Map->CVPtcpRatio);
			if (Record->isRecordingSnapshotCV && (!Record->SnapshotFirstLastSw || (Map->t <= Map->TransientEndTime+Record->NumSnapshot || Map->t > Map->RecordEndTime-Record->NumSnapshot))) Record->RecordSnapshotCV(Map);
			if (Record->isRecordingSmallestAngle || Record->isRecordingAllAngles) Record->AccAngleHist(Map);
			if (Record->isRecordingNeighAngle) Record->AccNeighAngleHist(Map);
			if (Record->isRecordingSubspaceAngle) Record->AccSubspaceAngleHist(Map);
			if (Record->isRecordingAngleTimeSeries) Record->RecordAngleTimeSeries(Map);
			if (Map->isRecordingPowerSpectra) Map->AccPowerSpectra();
			if (Record->isRecordingOutputQuantity && Map->ExistOutputQuantityBack && (!Record->OutputQuantityFirstLastSw || (Map->t <= Map->TransientEndTime+Record->NumOutputQuantity || Map->t > Map->RecordEndTime-Record->NumOutputQuantity))) {Map->CalcOutputQuantityBack(); Record->RecordOutputQuantityBack(Map);}
			if (Record->isRecordingSupplPtcpRatio) Record->CalcSupplPtcpRatio(Map);
			if (Map->isRecordingQuantityHist && Map->isQuantityHistBack) Map->AccQuantityHist();
			//Check exponents for covariant vectors (Caution! It changes X and CV);
			Map->CheckCVExp();
			Record->SetRecordTimeSeriesSw(Map, -1);
			if (Record->isRecordingCVExp && (Record->isRecordingLyapTimeSeriesAlways || Map->t - 1 == (Map->LoopTimes)[Map->LoopInd]))
				Record->RecordCVExp(Map);
			if (Record->isRecordingDOSViolation) Record->CountDOSViolation(Map);
			if (Record->isRecordingPtcpHist) Record->AccPtcpHist(Map);
			if (Map->isCheckingCVExp && Record->isRecordingLyapHist && Record->LyapHistGSCV) Record->AccLyapHist(Map);
			//Backward evolution;
			Map->BackwardEvolution();
			Map->t--;
			//Display progress bar;
			Record->ProgBar();
		}
	}
}


//main function;
int main(int argc, char *argv[]) {
	//Define variables;
	map *Map;
	record *Record = new record;
	//Initialization;
	omp::Init();
	//Input parameters;
	InputParameters(argc, argv, &Map, Record);
#if OMPSW
	cout << "using " << omp::NumThreads << " processors" << endl;
#endif
	//Set initial conditions;
	Map->SetInitialX();
	Map->SetInitialVar();
	//Main routine;
	switch (Map->BackwardMode) {
		case 0:
		case 2:
			Map->SignDynamics = 1.0;
			MainRoutineWithoutBackward(Map, Record);
			break;
		case 1:
			Map->SignDynamics = 1.0;
			if (Map->isDividingLoop) MainRoutineWithDivision(Map, Record);
				else MainRoutineWithoutDivision(Map, Record);
			break;
		case 3:
			Map->SignDynamics = -1.0;
			if (Map->isDividingLoop) MainRoutineInverseWithDivision(Map, Record);
				else MainRoutineInverseWithoutDivision(Map, Record);
			break;
		case 4:
			Map->SignDynamics = -1.0;
			if (Map->isDividingLoop) MainRoutineInverseCVWithDivision(Map, Record);
				else MainRoutineInverseCVWithoutDivision(Map, Record);
			break;
		case 5:
			Map->SignDynamics = 1.0;
			MainRoutineWithDivisionAfterTransient(Map, Record);
			break;
		default:
			ErrorExit("Undefined mode.");
			break;
	}
	//Record results;
	if (Map->BackwardMode == 2) Record->WriteResults(Map);
	Map->CalcResults();
	Record->RecordResults(Map);
	cout << endl;
	//Free memory;
	delete Map;
	delete Record;
	//End program;
	return 0;
}
