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
#include "gsl/gsl_multimin.h"

using std::cout;
using std::endl;
using std::list;
using std::string;
using std::vector;

std::ostream& operator << (std::ostream& out, const vector<int>& vec)
{
	for(auto it : vec)
		out << it << " ";
	return out;
}

typedef list<vector<int> > veclist;
void genNecklace(int n, int k, veclist& out)
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
	const gsl_multimin_fminimizer_type * mType;
	gsl_multimin_fminimizer * mMethod;
	gsl_vector * mSs;
	gsl_vector * mXvec;
	veclist mSymbols;
	void pruneRule();
    double minSrchFunc(const gsl_vector * x, void * params);
public:
	void init(int n);
	void mainLoop();
    LorentzGasElCells();
    ~LorentzGasElCells();
	
};

LorentzGasElCells::LorentzGasElCells(): mType(gsl_multimin_fminimizer_nmsimplex2)
{
}

void LorentzGasElCells::init(int n)
{
	// generate length n topolotical cycle out of 12 symbols
	genNecklace(n, 12, mSymbols);
}


double LorentzGasElCells::minSrchFunc(const gsl_vector * x, void * params)
{
	return 0;
}



int main(int argc, const char * argv[])
{

    // insert code here...
    //;
	veclist mylist;
	genNecklace(4, 12, mylist);
	for (auto it : mylist)
	{
		cout << it << endl;
	}
	cout << "done\n";
    return 0;
}

