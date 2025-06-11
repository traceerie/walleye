
/*
 * QMR.cpp
 *
 *  Created on: Jun 9, 2025
 *      Author: Trace
 */
using namespace std;

#include "QMR1.h"

#include <iostream>

#include <vector>

#include <cmath>

#include <limits>

const double EPS = 1e-3;

const int MAX_ITERS = 1000;


QMR1::QMR1() {
	// TODO Auto-generated constructor stub

}

QMR1::~QMR1() {
	// TODO Auto-generated destructor stub
}

// Matrix-vector product

vector<double> QMR1::matVecMul(const vector<vector<double>>& A, const vector<double>& x) {

    int n = A.size();

    vector<double> result(n, 0.0);

    for (int i = 0; i < n; ++i)

        for (int j = 0; j < x.size(); ++j)

            result[i] += A[i][j] * x[j];

    return result;

}



// Dot product

 double QMR1::dot(const vector<double>& a, const vector<double>& b) {

    double sum = 0.0;

    for (int i = 0; i < a.size(); ++i)

        sum += a[i] * b[i];

    return sum;

}
//in eclipse C++ how to create a library


// Norm

 double QMR1::norm(const vector<double>& v) {

    return sqrt(dot(v, v));

}

bool QMR1::hi(){
	cout<<"hi";
}

// QMR Algorithm

 bool QMR1::qmr(const vector<vector<double>>& A, const vector<double>& b, vector<double>& x) {

    int n = b.size();

    x.assign(n, 0.0);



    vector<double> r = b;

    vector<double> v_tld = r;

    vector<double> y = r;



    vector<double> w = r;

    vector<double> z = r;



    double rho = norm(r);

    double rho_prev = 1.0, alpha = 1.0, beta = 0.0, c = 1.0, c_prev = 1.0;

    double eta = -1.0, theta = 0.0, theta_prev = 0.0;

    vector<double> p(n, 0.0), q(n, 0.0), d(n, 0.0), s(n, 0.0);



    for (int iter = 0; iter < MAX_ITERS; ++iter) {

        if (rho < EPS) break;



        vector<double> v = v_tld;

        for (int i = 0; i < n; ++i)

            v[i] /= rho;



        vector<double> u = matVecMul(A, v);

        alpha = dot(u, v);



        if (fabs(alpha) < EPS) {

            cerr << "Breakdown: alpha too small." << endl;

            return false;

        }



        vector<double> u_hat(n);

        for (int i = 0; i < n; ++i)

            u_hat[i] = u[i] - alpha * v[i];



        rho_prev = rho;

        rho = norm(u_hat);

        if (rho < EPS) break;



        for (int i = 0; i < n; ++i)

            u_hat[i] /= rho;



        theta = rho / (alpha * c_prev);

        c = 1.0 / sqrt(1.0 + theta * theta);

        eta = c * c * alpha;

        if (iter == 0) {

            for (int i = 0; i < n; ++i)

                d[i] = eta * v[i];

        } else {

            for (int i = 0; i < n; ++i)

                d[i] = eta * v[i] + (theta_prev * theta_prev * c * c) * d[i];

        }



        for (int i = 0; i < n; ++i)

            x[i] += d[i];



        v_tld = u_hat;

        theta_prev = theta;

        c_prev = c;



        // Check convergence

        vector<double> r_new = b;

        vector<double> Ax = matVecMul(A, x);

        for (int i = 0; i < n; ++i)

            r_new[i] -= Ax[i];

        double res_norm = norm(r_new);



        cout << "Iter " << iter + 1 << " | Residual norm = " << res_norm << endl;

        if (res_norm < EPS)

            return true;

    }



    return false;

}


//int main() {

//	vector<vector<double>> A = {
//
//        {4, 1},
//
//        {1, 3}
//
//    };
//
//    vector<double> b = {1, 2};
//
//    vector<double> x;
//
//
//
//    if (qmr(A, b, x)) {
//
//        cout << "QMR converged.\nSolution x:\n";
//
//        for (double xi : x)
//
//            cout << xi << " ";
//
//        cout << endl;
//
//    } else {
//
//        cout << "QMR failed to converge." << endl;
//
//    }



//    return 0;
//
//}
//
//
// Sample usage



