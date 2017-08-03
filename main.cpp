#include <iostream>
#include <sstream>
#include "Parameters.h"
#include "numerov.h"

using namespace std;

int main() {

    Parameters params;

    ostringstream ds;
    ds<< "//home//martas//Dropbox//PracaMagisterska//Numerov1D//Wyniki//dListb=" << params.getB() << "alpha=" << params.getAlpha() << ".dat";
    ofstream f2(ds.str().c_str());
    for (int j = 0; j < params.n_d(); ++j) {
        ds << params.getD_min() + j * params.getDelta_d()<<'\n';
    }
    f2.close();

    int i;
    #pragma omp parallel for
    for (i = 0; i < params.n_d(); i++) {
        double d = params.getD_min() + i * params.getDelta_d();
        calculate(d, 3, params);
    }
}

