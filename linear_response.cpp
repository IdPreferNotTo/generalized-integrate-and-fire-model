//
// Created by lukassid on 19.02.20.
//

#include "linear_response.h"
#include "neuron_models.h"

double calculateQifLinResp(std::vector<double> &noises, double stepSize, double mu) {
    double linResp = 0;
    for (size_t x = 0; x < noises.size(); x++) {
        linResp += (1 / mu) * (1. / 2.) * (1. - cos(2 * sqrt(mu) * x * stepSize)) * noises[x] * stepSize;
    }
    return linResp;
}

double calculateLifLinResp(std::vector<double> &noises, double stepSize) {
    double linResp = 0;
    for (size_t x = 0; x < noises.size(); x++) {
        linResp += exp((x * stepSize)) * noises[x] * stepSize;
    }
    return linResp;
}

double
calculateGifLinResp(std::vector<double> &noises, double stepSize, double mu, double beta, double tauW, double delta,
                    double tDet, double aDet, double wDet) {
    double vT = 1;
    double nu = 1 + 1 / tauW;
    double omega = sqrt((beta + 1) / tauW - nu * nu / 4);
    double denom = mu - vT - beta * wDet - aDet + delta;

    double linResp = 0;
    for (size_t x = 0; x < noises.size(); x++) {
        double t = x * stepSize;
        linResp += noises[x] * stepSize * exp(nu * (t - tDet) / 2) *
                   (cos(omega * (t - tDet)) - ((1. - tauW) / (2 * tauW * omega)) * sin(omega * (t - tDet)));
    }
    return linResp / denom;
}
