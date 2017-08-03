//
// Created by martas on 10.01.17.
//

#ifndef NUMEROV1D_PARAMETERS_H
#define NUMEROV1D_PARAMETERS_H


class Parameters {
    double b;
    double fx_max;
    double alpha;
    double epsilon;
    double d_min;
    double d_max;
    double delta_d;
    double E_min;
    double E_max;

public:

    Parameters();

    int n_d() const;

    double dx(double E) const;

    double deltaE(double E) const;


    double getB() const {
        return b;
    }

    double getFx_max() const {
        return fx_max;
    }

    double getAlpha() const {
        return alpha;
    }

    double getEpsilon() const {
        return epsilon;
    }

    double getD_min() const {
        return d_min;
    }

    double getD_max() const {
        return d_max;
    }

    double getDelta_d() const {
        return delta_d;
    }

    double getE_min() const {
        return E_min;
    }

    double getE_max() const {
        return E_max;
    }
};


#endif //NUMEROV1D_PARAMETERS_H
