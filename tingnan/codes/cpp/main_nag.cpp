//
//  main.cpp
//  LorentzGas
//
//  Created by Hana&Tina on 5/19/14.
//  Copyright (c) 2014 Georgia Institute of Technology. All rights reserved.
//

#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <list>
#include <cmath>
#include <nag.h>
#include <nage04.h>
#include <nag_stdlib.h>

using std::cout;
using std::cerr;
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
        while(j < n - i - 1)
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
    template<typename T>
    class Vector2D
    {
    public:
        T x;
        T y;
        Vector2D() : x(T(0)), y(T(0)) {};
        Vector2D(T xx) : x(xx), y(xx) {};
        Vector2D(T xx, T yy) : x(xx), y(yy) {};
        T magnitude() const { return sqrt(x * x + y * y); };
        void normalize() { T mag = magnitude(); if (mag) *this *= 1 / mag; };
        T dot(const Vector2D<T> &v) const { return x * v.x + y * v.y; };
        Vector2D operator + (const Vector2D<T> &v) const { return Vector2D<T>(x + v.x, y + v.y); };
        Vector2D operator * (const T &val) const { return Vector2D<T>(x * val, y * val); };
        Vector2D operator / (const T &val) const { T invVal = T(1) / val; return Vector2D<T>(x * invVal, y * invVal); };
        Vector2D operator / (const Vector2D<T> &v) const { return Vector2D<T>(x / v.x, y / v.y); };
        Vector2D operator * (const Vector2D<T> &v) const { return Vector2D<T>(x * v.x, y * v.y); };
        Vector2D operator - (const Vector2D<T> &v) const { return Vector2D<T>(x - v.x, y - v.y); };
        Vector2D operator - () const { return Vector2D<T>(-x, -y); }
        Vector2D& operator *= (const T &val) { x *= val; y *= val; return *this; };
    };
    typedef Vector2D<double> vec2;
private:
    std::ofstream myFile;
    std::ofstream myTest;
    static const int mNdim = 2;
    int mNsyms;
    iveclist mSymbols; // initial list of symbols
    vector<int>* mPtrCurSymbols; // to be used internally with pathLength
    iveclist mLabels; // final list of symbols
    dveclist mThetas; // final list of thetas

    // geometric params
    const double mRadius = 1;
    const double mWidth = 0.05;
    vector<vec2> mCenterList;
    vector<double> mCenterListX;
    vector<double> mCenterListY;
    // internal functions 
    bool pruneRule(const vector<int>& curSymbols);
    bool testLink(const vector<int>& curSymbols, const vector<double>& thetas);
    double pathLength(const double* x);
public:
    void mainLoop(int n);
    LorentzGasElCells();
    ~LorentzGasElCells() {};
    // for calling the gsl routine
    static void minSrchFunc(Integer n, const double *xc, double *fc, Nag_Comm* comm)
    {
        // n is not used here because we internnaly recorded
        *fc = reinterpret_cast<LorentzGasElCells*>(comm->p)->pathLength(xc);
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
    vector<int> sbvec(mNsyms + 1, curSymbols[0]);
    sbvec.assign(curSymbols.begin(), curSymbols.end());
    int i = 0;
    for( ; i < mNsyms; ++i)
    {
        if (sbvec[i] == sbvec[i + 1])
            return true;
        if (sbvec[i] % 2 == 0 && abs(sbvec[i + 1] - sbvec[i]) < 2)
            return true;
        if (sbvec[i] % 2 == 1 && abs(sbvec[i + 1] - sbvec[i]) < 3)
            return true;
    }
    return false;
}

// this is going to be our path function
// depends on internal symbol list
double LorentzGasElCells::pathLength(const double * thetas)
{
    double sLength = 0;
    vector<double> thvec(mNsyms + 1, thetas[0]);
    thvec.assign(thetas, thetas + mNsyms);
    int idx = 0;
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
    return sLength;
}

/*
static void minSrchFuncTest(Integer n, const double *xc, double *fc, Nag_Comm* comm)
{
    // n is not used here because we internnaly recorded
    *fc = exp(xc[0]) * (xc[0] * 4.0 * (xc[0] + xc[1]) +
                        xc[1] * 2.0 * (xc[1] + 1.0) + 1.0);
};
*/

bool LorentzGasElCells::testLink(const vector<int>& curSymbols, const vector<double>& thetas)
{
    // loop assignment
    vector<double> thvec(mNsyms + 1, thetas[0]);
    thvec.assign(thetas.begin(), thetas.end());
    vector<int> sbvec(mNsyms + 1, curSymbols[0]);
    sbvec.assign(curSymbols.begin(), curSymbols.end());
    int idx = 0;
    for(; idx < mNsyms; ++idx)
    {
        int sbl = sbvec[idx]; 
        vector<vec2> dkList(4);
        dkList[0] = vec2();
        dkList[1] = vec2(mCenterListX[sbl], mCenterListY[sbl]);
        // with tmpidx we can find the seperation of the two disks for the next bounce
        vec2 pt1(mRadius * cos(thvec[idx]), mRadius * sin(thvec[idx]));
        vec2 pt2(mRadius * cos(thvec[idx + 1]), mRadius * sin(thvec[idx + 1]));
        vec2 segHat = dkList[1] + pt2 - pt1;
        double segLength = segHat.magnitude();
        segHat.normalize();
        vector<vec2> crList(4);
        crList[0] = -pt1;
        crList[1] = dkList[1] - pt1;
        auto mod = [](int a, int b) { int m = a % b; return m < 0 ? m + b : m;};
        int i2, i3;
        if(sbl % 2 == 1)
        {
            i2 = mod(sbl - 1, 12);
            i3 = mod(sbl + 1, 12);
        }
        else
        {
            i2 = mod(sbl - 2, 12);
            i3 = mod(sbl + 2, 12);
        }
        dkList[2] = vec2(mCenterListX[i2], mCenterListY[i2]);
        dkList[3] = vec2(mCenterListX[i3], mCenterListY[i3]);
        crList[2] = dkList[2] - pt1;
        crList[3] = dkList[3] - pt1;
        int lidx = 0;
        for(; lidx < 4; ++lidx)
        {
            double projLength = segHat.dot(crList[lidx]);
            if (projLength <= 0 || projLength >= segLength);
            else
            {
                vec2 projVec = segHat * projLength;
                vec2 distVec = pt1 + projVec - dkList[lidx];
                double distVert = distVec.magnitude();
                if (distVert < mRadius)
                {
                    /*
                    cout << "debug: neigh disk symbosl: ";
                    cout << idx << " " << i2 << " " << i3 << endl;
                    cout << "debug: problem disk: ";
                    cout << lidx << " " << distVert << " " << projLength << " " << segLength << endl;
                    cout << "debug: other: ";
                    cout << curSymbols << " " << thetas << endl;
                    */
                    return false;
                }
            }
        }
    }

    return true;
}

void LorentzGasElCells::mainLoop(int n)
{
    mNsyms = n;
    // generate length n topolotical cycle out of 12 symbols
    genNecklace(n, 12, mSymbols);
    std::string fname = "l" + std::to_string(n) + ".txt";
    myFile.open(fname.c_str());
    myTest.open("converge.txt");
    cout << "init success" << endl;
    mSymbols.remove_if([this](const vector<int>& curSymbols){return pruneRule(curSymbols);});
    //
    NagError fail;
    Nag_Comm comm;
    Nag_E04_Opt options;
    double objf, *x;
    vector<double> xvec(mNsyms, 0.0);
    x = &xvec[0];
    INIT_FAIL(fail);
    // automatically initialized;
    nag_opt_init(&options);
    // options.print_fun = monit;
    options.print_level = Nag_NoPrint;
    options.max_iter = 100000;
    // options.optim_tol = 1e-10;
    // nag_opt_read("e04ccc", "config.txt", &options, Nag_FALSE, "stdout", &fail);
    if (fail.code != NE_NOERROR)
    {
        cout << fail.message << endl;
        return;
    }
    comm.p = reinterpret_cast<void*>(this);
    for (auto& pSymbol : mSymbols)
    {
        mPtrCurSymbols = &pSymbol;
        xvec.assign(mNsyms, 0.0);
        // set the current symbol to evaluate
        nag_opt_simplex(mNsyms, minSrchFunc, x, &objf, &options, &comm, &fail);
        bool noIntersect = testLink(pSymbol, xvec);
        if (fail.code != NE_NOERROR)
        {
            cerr << fail.message << endl;
            vector<double> tmpvec(mNsyms);
            int i = 0;
            for ( ; i < mNsyms; ++i)
            {
                tmpvec[i] = x[i];
            }
            myTest << pSymbol << " ";
            myTest << std::fixed << std::setprecision(14) << tmpvec << endl;
        }
        if (fail.code == NE_NOERROR && noIntersect)
        {

            vector<double> tmpvec(mNsyms);
            int i = 0;
            for ( ; i < mNsyms; ++i)
            {
                tmpvec[i] = x[i];
            }
            // mLabels.push_back(pSymbol);
            // mThetas.push_back(tmpvec);
            cout << "current symbols: " << pSymbol << " ";
            cout << tmpvec << endl;
            myFile << pSymbol << " ";
            myFile << std::fixed << std::setprecision(14) << tmpvec << endl;
        }
        // check if it is a physically valid solution by testing intersection
        // store the solution and symbol sequence as well
    }
    myFile.close();
    myTest.close();
    mSymbols.resize(0);
}


int main(int argc, const char * argv[])
{
    /*
    LorentzGasElCells billiardSystem;
    for(int i = 2; i < 9; ++i)
    {
        billiardSystem.mainLoop(i); 
    }
    */
    for (int i = 2; i < 3; ++i)
    {
        string fname =  "w=0.30/l" + std::to_string(i) + ".txt";
        std::ifstream inputFile(fname);
        while (inputFile.good())
        {
            string strvalue;
            getline(inputFile, strvalue, '\n');
            std::istringstream line_stream(strvalue);
            for (int j = 0; j < i; ++j)
            {
                string symbolval;
                getline(line_stream, symbolval, ' ');
                cout << symbolval << "\t";
            }
            cout << endl;
        }
    }
    return 0;
}

