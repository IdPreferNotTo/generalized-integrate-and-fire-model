//
// Created by lukassid on 19.02.20.
//

#ifndef SCH13_LINEAR_RESPONSE_H
#define SCH13_LINEAR_RESPONSE_H

#include <cmath>
#include <vector>

double
calculate_qif_lin_resp(std::vector<double> &noises, double stepSize, double mu);

double
calculate_theta_lin_resp(std::vector<double> &noises, double stepSize, double t_det);

double
calculate_lif_lin_resp(std::vector<double> &noises, double stepSize);

double
calculate_gif_lin_resp(std::vector<double> &noises, double stepSize, double mu, double beta, double tau_w, double Delta,
                       double t_det, double a_det, double w_det);

#endif //SCH13_LINEAR_RESPONSE_H
