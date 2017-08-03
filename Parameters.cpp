//
// Created by martas on 10.01.17.
//

#include "Parameters.h"
int Parameters::n_d() const {
    return (d_max - d_min) / delta_d;
}

double Parameters::dx(double E) const{
    if (E < 0)
        return 0.00001;
    if (E >= 0)
        return 0.0002;
}

Parameters::Parameters() {
    b = 0.2;
    fx_max = 0.000001;
    alpha = 0.002;
    epsilon = 0.0001;
    d_min = 0.05;
    d_max = 3;
    delta_d = 0.05;
    E_min = -0.5;
    E_max = 0.2;
}

double Parameters::deltaE(double E) const{
    if (E < -20) return 10;
    else if (E >= -20) return 0.005;}