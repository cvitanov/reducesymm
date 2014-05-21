//
//  main.cpp
//  LorentzGas
//
//  Created by Hana&Tina on 5/19/14.
//  Copyright (c) 2014 Georgia Institute of Technology. All rights reserved.
//

#include <iostream>
#include <string>
#include <vector>
#include <list>
#include <cmath>
#include <nag.h>
#include <nage04.h>
#include <nag_stdlib.h>

using std::cout;
using std::endl;
using std::list;
using std::string;
using std::vector;

template<typename T>
std::ostream& operator << (std::ostream& out, const vector<T>& vec)
{
    for(auto& it : vec)
        out << it << " ";
    return out;
}

typedef list<vector<int> > iveclist;
typedef list<vector<double> > dveclist;
void genNecklace(int n, int k, iveclist& out)
{
    
    int i = 0;
    vector<int> tmp(n);
    for(; i < n; ++i)
        tmp[i] = 0;
    out.push_back(tmp);
    i = n - 1;
    do
    {
        tmp[i] = tmp[i] + 1;
        int j = 0;
        while(j < n - i)
        {
            tmp[j + i + 1] = tmp[j];
            ++j;
        }
        if (n % (i + 1) == 0)
            out.push_back(tmp);
        i = n - 1;
        while( tmp[i] == k - 1)
            i = i - 1;
    } while(i >= 0);
}

class LorentzGasElCells
{
private:
    static const int mNdim = 2;
    int mNsyms;
    iveclist mSymbols; // initial list of symbols
    vector<int>* mPtrCurSymbols; // to be used internally with pathLength
    iveclist mLabels; // final list of symbols
    dveclist mThetas; // final list of thetas

    // geometric params
    const double mRadius = 1;
    const double mWidth = 0.3;
    vector<double> mCenterListX;
    vector<double> mCenterListY;
    // internal functions 
    bool pruneRule(const vector<int>& curSymbols);
    bool testLink();
    double pathLength(const double* x);
public:
    void init(int n);
    void mainLoop();
    LorentzGasElCells();
    ~LorentzGasElCells() {};
    // for calling the gsl routine
    static void minSrchFunc(Integer n, const double *xc, double *fc, Nag_Comm* comm)
    {
        // n is not used here because we internnaly recorded
        *fc = static_cast<LorentzGasElCells*>(comm->p)->pathLength(xc);
    };
};

LorentzGasElCells::LorentzGasElCells(): mNsyms(0)
{
    const double M_SQRT3 = 1.73205080756887729352744634151;
    mCenterListX.resize(12);
    mCenterListY.resize(12);
    const double distEven = 2 * mRadius + mWidth;
    const double distOdd = distEven * M_SQRT3;
    // the x array
    const double arrayX[12] = {distEven, distOdd * M_SQRT3 / 2, distEven / 2, 0, -distEven / 2, -distOdd * M_SQRT3 / 2, -distEven, -distOdd * M_SQRT3 / 2, -distEven / 2, 0, distEven / 2, distOdd * M_SQRT3 / 2};
    mCenterListX.assign(arrayX, arrayX + 12);
    // the y array
    const double arrayY[12] = {0, distOdd / 2, distEven * M_SQRT3 / 2, distOdd, distEven * M_SQRT3 / 2, distOdd / 2, 0, -distOdd / 2, -distEven * M_SQRT3 / 2, -distOdd, -distEven * M_SQRT3 / 2, -distOdd / 2};
    mCenterListY.assign(arrayY, arrayY + 12);
    // cout << mCenterListX << endl;
    // cout << mCenterListY << endl;
}

bool LorentzGasElCells::pruneRule(const vector<int>& curSymbols)
{
    // do sth to the mSymbols;
    return false;
}

void LorentzGasElCells::init(int n)
{
    mNsyms = n;
    // generate length n topolotical cycle out of 12 symbols
    genNecklace(n, 12, mSymbols);
    mSymbols.remove_if([this](const vector<int>& curSymbols){return pruneRule(curSymbols);});
    /*
    for (auto it : mSymbols)
    {
        cout << it << endl;
    }
    */
    // the path depend on Ndim * n variables
}

// this is going to be our path function
// depends on internal symbol list
double LorentzGasElCells::pathLength(const double * thetas)
{
    double sLength = 0;
    int idx = 0;
    vector<double> thvec(mNsyms + 1, thetas[0]);
    thvec.assign(thetas, thetas + mNsyms);
    for (; idx < mNsyms; ++idx)
    {
        // do sth, according to symbols
        int tmpidx = (*mPtrCurSymbols)[idx]; 
        // with tmpidx we can find the seperation of the two disks for the next bounce
        double diskX = mCenterListX[tmpidx];
        double diskY = mCenterListY[tmpidx];
        double deltaX = mRadius * (cos(thvec[idx + 1]) - cos(thvec[idx]));
        double deltaY = mRadius * (sin(thvec[idx + 1]) - sin(thvec[idx]));
        double tmpLength = sqrt((deltaX + diskX) * (deltaX + diskX) + (deltaY + diskY) * (deltaY + diskY));
        sLength = sLength + tmpLength;

    }
    sLength = (thetas[0] - 2) * (thetas[0] - 2) + (thetas[1] - 1) * (thetas[1] - 1);
    return sLength;
}

void monit(const Nag_Search_State *st, Nag_Comm *comm)
{
    cout << st->iter << endl;
}

void LorentzGasElCells::mainLoop()
{
    NagError fail;
    Nag_Comm comm;
    Nag_E04_Opt options;
    double objf, *x;
    vector<double> xvec(mNsyms, 0.0);
    x = &xvec[0];

    INIT_FAIL(fail);
    // automatically initialized;
    nag_opt_init(&options);
    options.print_fun = monit;
    nag_opt_read("e04ccc", "config.txt", &options, Nag_TRUE, "stdout", &fail);
    if (fail.code != NE_NOERROR)
    {
        cout << fail.message << endl;
        return;
    }
    comm.p = static_cast<void*>(this);
    for (auto& pSymbol : mSymbols)
    {
        cout << "current symbols: " << pSymbol << endl;
        // set the current symbol to evaluate
        nag_opt_simplex(mNsyms, minSrchFunc, x, &objf, &options, NAGCOMM_NULL, &fail);
        cout << "done" << endl;
        if (fail.code == NE_NOERROR)
        {

            vector<double> tmpvec(mNsyms);
            int i = 0;
            for ( ; i < mNsyms; ++i)
            {
                tmpvec[i] = x[i];
            }
            mLabels.push_back(pSymbol);
            mThetas.push_back(tmpvec);
            cout << tmpvec << endl;
        }
        // check if it is a physically valid solution by testing intersection
        // store the solution and symbol sequence as well
    }
}


int main(int argc, const char * argv[])
{
    LorentzGasElCells billiardSystem;
    billiardSystem.init(2);
    billiardSystem.mainLoop();
    cout << "done\n";
    return 0;
}

