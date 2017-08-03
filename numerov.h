//
// Created by martas on 11.01.17.
//numerov dla dwoch jonow
//

#ifndef NUMEROV1D_NUMEROV_H
#define NUMEROV1D_NUMEROV_H

#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <sstream>
#include"Parameters.h"

using namespace std;

double v(double *z, int i, double d, const Parameters &params);

double r_max(double E, const Parameters &params, double d);

int n(double E, const Parameters &params, double d);

void num_left(double E, int N, double *z, double *fl, double d, const Parameters &params);

void num_right(double E, double *z, double *fr, int *wskN, double d, const Parameters &params);

double help_function(double f_right_r, double f_right_dr, double f_left_r, double f_left_dr);

double HF(double E, double d, const Parameters &params);

double bisection(double E_1, double E_2, double d, const Parameters &params);

void calculate(double d, const Parameters &params);

double v(double *z, int i, double d, const Parameters &params) {
    double trap = params.getAlpha() * z[i] * z[i];
    double interaction =
            -1.0 / (((z[i] + d) * (z[i] + d) + params.getB() * params.getB()) *
                    ((z[i] + d) * (z[i] + d) + params.getB() * params.getB()))
            - 1.0 / (((z[i] - d) * (z[i] - d) + params.getB() * params.getB()) *
                     ((z[i] - d) * (z[i] - d) + params.getB() * params.getB()));
    return trap + interaction;
}

double r_max(double E, const Parameters &params, double d) {
    if (E <= -0.1) return 2 + d;
    else if (E > -0.1 && E < 0.1) return 20 + d;
    else if (E >= 0.1 && E < 0.3) return 30 + d;
    else if (E >= 0.3 && E < 0.5) return 2.2 * sqrt(E / params.getAlpha()) + d;
    else if (E >= 0.5 && E < 4) return (2 - 0.1 * E) * sqrt(E / params.getAlpha()) + d;
    else if (E >= 4 && E < 10) return (2 - 0.05 * E) * sqrt(E / params.getAlpha()) + d;
    else if (E >= 10) return (2 - 0.01 * E) * sqrt(E / params.getAlpha()) - 40 + d;


}

int n(double E, const Parameters &params, double d) {
    int nn = (2 * r_max(E, params, d)) / params.dx(E);
    if (nn % 2 == 0) nn = nn + 1;
    return nn;
}

void num_left(double E, int N, double *z, double *fl, double d, const Parameters &params) {
    z[0] = -r_max(E, params, d);
    fl[0] = -params.getFx_max();

    z[1] = -r_max(E, params, d) + params.dx(E);
    fl[1] = -0.00001;

    double dr2 = params.dx(E) * params.dx(E);

    for (int i = 2; i < N + 1; i++) {
        z[i] = -r_max(E, params, d) + i * params.dx(E);
        double q_i = -v(z, i, d, params) + E;  //Q(i+1)
        double q_i1 = (-v(z, i - 1, d, params) + E);    // Q(i)
        double q_i2 = (-v(z, i - 2, d, params) + E);    //Q(i-1)


        fl[i] = ((2 - 5 * dr2 * q_i1 / 6) * fl[i - 1] - (1 + dr2 * q_i2 / 12) * fl[i - 2]) / (1 + dr2 * q_i / 12);
    }

}

void num_right(double E, double *z, double *fr, int *wskN, double d, const Parameters &params) {

    int nn = n(E, params, d);
    z[nn - 1] = r_max(E, params, d);
    fr[nn - 1] = params.getFx_max();

    z[nn - 2] = r_max(E, params, d) - params.dx(E);


    fr[nn - 2] = 0.00001;//0.5*pow(k, -0.5)*exp(-abs(k)*(z[n-2]));

    double dr2 = params.dx(E) * params.dx(E);

    for (int i = nn - 3; i > 0; i--) {
        z[i] = r_max(E, params, d) - (nn - i) * params.dx(E);

        double q_i = -v(z, i, d, params) + E;  //Q(i+1)
        double q_i1 = (-v(z, i + 1, d, params) + E);    // Q(i)
        double q_i2 = (-v(z, i + 2, d, params) + E);    //Q(i-1)

        fr[i] = ((2 - 5 * dr2 * q_i1 / 6) * fr[i + 1] - (1 + dr2 * q_i2 / 12) * fr[i + 2]) / (1 + dr2 * q_i / 12);

        if (fr[i] < fr[i + 1]) {
            *wskN = i + 1;
            break;
        }
    }

}

double help_function(double f_right_r, double f_right_dr, double f_left_r, double f_left_dr) {
    double help;
    help = f_right_r / f_right_dr - f_left_r / f_left_dr;
    return help;
}

double HF(double E, double d, const Parameters &params) {

    double *r = NULL;
    r = new double[n(E, params, d)];


    int N = 0;
    double psi_right_r, psi_right_dr, psi_left_r, psi_left_dr;

    double *psi_right;
    psi_right = new double[n(E, params, d)];

    num_right(E, r, psi_right, &N, d, params);

    if (N == 0) N = 1;


    psi_right_r = psi_right[N - 1];
    psi_right_dr = psi_right[N];

    double *psi_left = NULL;
    psi_left = new double[N + 1];

    num_left(E, N, r, psi_left, d, params);

    psi_left_r = psi_left[N - 1];
    psi_left_dr = psi_left[N];


    double hf = help_function(psi_right_r, psi_right_dr, psi_left_r, psi_left_dr);

    delete psi_left;
    delete psi_right;
    delete r;


    return hf;


}

double bisection(double E_1, double E_2, double d, const Parameters &params) {
    double s;
    if (HF(E_1, d, params) * HF(E_2, d, params) < 0) {
        s = E_1;
        while ((E_2 - E_1) > params.getEpsilon()) {
            s = (E_1 + E_2) / 2;
            if (HF(E_1, d, params) * HF(s, d, params) < 0) {
                E_2 = s;
            } else if (HF(E_2, d, params) * HF(s, d, params) < 0) {
                E_1 = s;
            } else break;
        }
    } else {
        s = 0;;
    }
    return s;
}

void calculate(double d, const Parameters &params) {


    int N = 0;
    double psi_right_r, psi_right_dr, psi_left_r, psi_left_dr;


    ostringstream energy_levels;
    energy_levels << "//home//martas//Dropbox//PracaMagisterska//Numerov1D//Wyniki//Ed=" << d << "b=" << params.getB()
                  << "alpha=" << params.getAlpha() << ".dat";
    ofstream f2(energy_levels.str().c_str());

    ostringstream plot_eigenfunctions;
    plot_eigenfunctions << "//home//martas//Dropbox//PracaMagisterska//Numerov1D//Wyniki//plotEigenfd=" << d << "alpha="
                        << params.getAlpha() << ".txt";
    ofstream f3(plot_eigenfunctions.str().c_str());


    f3 << "set terminal png size 1000, 600";
    f3 << "\n";
    f3 << "unset key";
    f3 << "\n";

    int s = 0;
    for (double E = params.getE_min(); E < params.getE_max(); E = E + params.deltaE(E)) {

        double K = HF(E, d, params);


        if (K / HF(E + params.deltaE(E), d, params) < 0) {
            s = s + 1;
            if (abs(HF(E - params.deltaE(E), d, params)) / abs(K) > 1) {
                double E_1 = E;
                double E_2 = E + params.deltaE(E);

                double eig = bisection(E_1, E_2, d, params);
                double *r = NULL;
                r = new double[n(eig, params, d)];

                std::cout << "Eig(d = " << d << "): " << setprecision(8) << eig << ", r_max: " << r_max(E, params, d) <<
                          endl;


                f2 << eig;
                f2 << "\n";

                double *psi_right = NULL;
                psi_right = new double[n(eig, params, d)];

                num_right(eig, r, psi_right, &N, d, params);
                psi_right_r = psi_right[N - 1];
                psi_right_dr = psi_right[N];

                double *psi_left = NULL;
                psi_left = new double[N + 1];

                num_left(eig, N, r, psi_left, d, params);
                psi_left_r = psi_left[N - 1];
                psi_left_dr = psi_left[N];

                if (d == 1 || d == 5 || d == 10 || d == 20 || d == 25 || d == 15) {
                    ostringstream function_values;
                    function_values << "//home//martas//Dropbox//PracaMagisterska//Numerov1D//Wyniki//psid=" << d << "E"
                                    << eig << ".dat";
                    ofstream f4(function_values.str().c_str());


                    f3 << "set output \"psi_d=";
                    f3 << d;
                    f3 << "_";
                    f3 << eig;
                    f3 << ".png\" ";
                    f3 << "\n";
                    f3 << "set title ";
                    f3 << "\"l = ";
                    f3 << d;
                    f3 << ", E = ";
                    f3 << eig;
                    f3 << "E*\" ";
                    f3 << "\n";
                    f3 << "plot \"/home//martas//Dropbox//PracaMagisterska//Numerov1D//Wyniki//psid=";
                    f3 << d;
                    f3 << "E";
                    f3 << eig;
                    f3 << ".dat\" with lines";
                    f3 << "\n";
                    f3 << "unset title";
                    f3 << "\n";
                    f3 << "unset output";
                    f3 << "\n";


                    int add = 10;
                    if (r_max(eig, params, d) > 50) add = 10;
                    for (int i = 0; i < N; i = i + add) {

                        f4 << r[i];

                        f4 << "\t";
                        f4 << psi_left[i] / psi_left_r;
                        f4 << "\n";

                    }
                    for (int i = N; i < n(eig, params, d); i = i + add) {
                        f4 << r[i];
                        f4 << "\t";
                        f4 << psi_right[i] / psi_right_r;
                        f4 << "\n";
                    }

                    f4.close();
                }
                delete psi_right;
                delete psi_left;
                delete r;
            }
        }
    }

    f3 << "unset terminal";
    f3.close();
    f2.close();


}


#endif //NUMEROV1D_NUMEROV_H
