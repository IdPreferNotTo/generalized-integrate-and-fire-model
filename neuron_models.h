//
// Created by lukassid on 15.07.19.
//

#ifndef SCH13_NEURON_MODELS_H
#define SCH13_NEURON_MODELS_H

#include <cmath>

const double pi = 3.14159265358979323846;

bool
QIF(double *v, double *t, double gamma, double mu, double vT, double D, double xi, double stepSize);


bool
cnoise_QIF(double *v, double *n, double *t, double gamma, double mu, double vT, double tau_n, double D_n, double D_w,
           double xi_n, double xi_w, double stepSize);


bool
adap_QIF(double *v, double *a, double *t, double vT, double gamma, double mu, double tau_a, double delta,
         double D, double xi, double stepSize);


bool
cnoise_adap_QIF(double *v, double *a, double *n, double *t, double vT, double gamma, double mu,
                double tau_a, double delta, double tau_n, double D_n, double D_w, double xi_n, double xi_w,
                double stepSize);


bool
theta(double *theta, double *t, double D, double xi, double stepSize);


bool
cnoise_theta(double *theta, double *n, double *t, double tau_n, double D, double xi, double stepSize);


bool
cnoise_adap_theta(double *theta, double *a, double *n, double *t, double tau_a,
                  double delta, double tau_n, double D, double xi, double stepSize);


bool
LIF(double *v, double *t, double gamma, double mu, double D, double xi, double stepSize);


bool
adap_LIF(double *v, double *a, double *t, double gamma, double mu, double tau_a, double delta, double D, double xi,
         double stepSize);


bool
cnoise_LIF(double *v, double *n, double *t, double gamma, double mu, double tau_n,
           double D_n, double D_w, double xi_n, double xi_w, double stepSize);

bool
cnoise_adap_LIF(double *v, double *a, double *n, double *t, double gamma, double mu, double tau_a, double delta,
                double tau_n, double D_n, double D_w, double xi_n, double xi_w, double stepSize);

bool
GIF(double *v, double *w, double *t, double gamma, double mu, double beta, double tau_w, double w_reset,
    double D, double xi, double stepSize);


bool
adap_GIF(double *v, double *w, double *a, double *t, double gamma, double mu, double beta, double tau_w, double w_reset,
         double tau_a,
         double delta, double D, double xi, double stepSize);


bool
cnoise_GIF(double *v, double *w, double *n, double *t, double gamma, double mu, double beta, double tau_w,
           double w_reset, double tau_n, double D_n, double D_w, double xi_n, double xi_w, double stepSize);


bool
cnoise_adap_GIF(double *v, double *w, double *a, double *n, double *t, double gamma, double mu, double beta,
                double tau_w, double w_reset, double tau_a, double delta, double tau_n, double D_n, double D_w,
                double xi_n, double xi_w, double stepSize);

#endif //SCH13_NEURON_MODELS_H
