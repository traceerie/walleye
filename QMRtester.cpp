
//============================================================================
// Name        : ff.cpp
// Author      : v
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================
using namespace std;

#include <iostream>

#include "QMR1.h"


#include <iostream>

#include <vector>

#include <cmath>

#include <limits>



int main() {


    QMR1 a;

	vector<vector<double>> A = {

        {4, 1},

        {1, 3}

    };

    vector<double> b = {1, 2};

    vector<double> x;



    if (a.qmr(A,b,x)) {

        cout << "QMR converged.\nSolution x:\n";

        for (double xi : x)

            cout << xi << " ";

        cout << endl;

    } else {

        cout << "QMR failed to converge." << endl;

    }


}
