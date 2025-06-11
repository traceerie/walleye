/*
 * QMR.h
 *
 *  Created on: Jun 9, 2025
 *      Author: Trace
 */
using namespace std;

#ifndef QMR1_H_
#define QMR1_H_
#include <iostream>

#include <vector>

#include <cmath>

#include <limits>

class QMR1 {
public:
     QMR1();
    virtual  ~QMR1();
    vector<double> matVecMul(const vector<vector<double>>& A, const vector<double>& x);
	bool qmr(const vector<vector<double>>& A, const vector<double>& b, vector<double>& x);
	bool hi();
	 double dot(const vector<double>& a, const vector<double>& b);
	 double norm(const vector<double>& v);
};

#endif /* QMR1_H_ */
