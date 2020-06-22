//
// Created by lukassid on 19.02.20.
//

#include "linear_response.h"
#include "neuron_models.h"

double calculate_qif_lin_resp(std::vector<double> &noises, double stepSize, double mu) {
    double linResp = 0;
    for (size_t x = 0; x < noises.size(); x++) {
        linResp += (1 / mu) * (1. / 2.) * (1. - cos(2 * sqrt(mu) * x * stepSize)) * noises[x] * stepSize;
    }
    return linResp;
}

double calculate_theta_lin_resp(std::vector<double> &noises, double stepSize, double t_det) {
    double linResp = 0;
    int precision = 100;
    int totalSteps = noises.size();
    for (size_t x = 0; x < totalSteps; x += precision) {
        linResp += (1 - cos(2 * pi * (float(x) / totalSteps))) * noises[x] * precision * stepSize;
    }
    return linResp;
}

double calculate_lif_lin_resp(std::vector<double> &noises, double stepSize) {
    double linResp = 0;
    for (size_t x = 0; x < noises.size(); x++) {
        linResp += exp((x * stepSize)) * noises[x] * stepSize;
    }
    return linResp;
}

double
calculate_gif_lin_resp(std::vector<double> &noises, double stepSize, double mu, double beta, double tau_w, double Delta,
                       double t_det, double a_det, double w_det) {
    double vT = 1;
    double nu = 1 + 1 / tau_w;
    double omega = sqrt((beta + 1) / tau_w - nu * nu / 4);
    double denom = mu - vT - beta * w_det - a_det + Delta;

    double linResp = 0;
    for (size_t x = 0; x < noises.size(); x++) {
        double t = x * stepSize;
        linResp += noises[x] * stepSize * exp(nu * (t - t_det) / 2) *
                   (cos(omega * (t - t_det)) - ((1. - tau_w) / (2 * tau_w * omega)) * sin(omega * (t - t_det)));
    }
    return linResp / denom;
}
