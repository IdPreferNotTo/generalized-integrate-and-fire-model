//
// Created by lukassid on 15.07.19.
//

#ifndef SCH13_NEURON_MODELS_H
#define SCH13_NEURON_MODELS_H

#include <cmath>

const double pi = 3.14159265358979323846;

bool
qif(double *v, double *t, double gamma, double mu, double vT, double d, double xi, double stepSize);


bool
cnoiseQif(double *v, double *n, double *t, double gamma, double mu, double vT, double tauN, double dN, double dW,
          double xiN, double xiW, double stepSize);


bool
adapQif(double *v, double *a, double *t, double vT, double gamma, double mu, double tauA, double delta,
        double d, double xi, double stepSize);


bool
cnoiseAdapQif(double *v, double *a, double *n, double *t, double vT, double gamma, double mu,
              double tauA, double delta, double tauN, double dN, double dW, double xiN, double xiW,
              double stepSize);

bool
lif(double *v, double *t, double gamma, double mu, double d, double xi, double stepSize);


bool
adapLif(double *v, double *a, double *t, double gamma, double mu, double tauA, double delta, double d, double xi,
        double stepSize);


bool
cnoiseLif(double *v, double *n, double *t, double gamma, double mu, double tauN,
          double dN, double dW, double xiN, double xiW, double stepSize);

bool
cnoiseAdapLif(double *v, double *a, double *n, double *t, double gamma, double mu, double tauA, double delta,
              double tauN, double dN, double dW, double xiN, double xiW, double stepSize);

bool
gif(double *v, double *w, double *t, double gamma, double mu, double beta, double tauW, double wReset,
    double d, double xi, double stepSize);


bool
adapGif(double *v, double *w, double *a, double *t, double gamma, double mu, double beta, double tauW, double wReset,
        double tauA,
        double delta, double d, double xi, double stepSize);


bool
cnoiseGif(double *v, double *w, double *n, double *t, double gamma, double mu, double beta, double tauW,
          double wReset, double tauN, double dN, double dW, double xiN, double xiW, double stepSize);


bool
cnoiseAdapGif(double *v, double *w, double *a, double *n, double *t, double gamma, double mu, double beta,
              double tauW, double wReset, double tauA, double delta, double tauN, double dN, double dW,
              double xiN, double xiW, double stepSize);

#endif //SCH13_NEURON_MODELS_H
