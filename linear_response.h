//
// Created by lukassid on 19.02.20.
//

#ifndef SCH13_LINEAR_RESPONSE_H
#define SCH13_LINEAR_RESPONSE_H

#include <cmath>
#include <vector>

double
calculateQifLinResp(std::vector<double> &noises, double stepSize, double mu);

double
calculateThetaLinResp(std::vector<double> &noises, double stepSize, double tDet);

double
calculateLifLinResp(std::vector<double> &noises, double stepSize);

double
calculateGifLinResp(std::vector<double> &noises, double stepSize, double mu, double beta, double tauW, double delta,
                    double tDet, double aDet, double wDet);

#endif //SCH13_LINEAR_RESPONSE_H
