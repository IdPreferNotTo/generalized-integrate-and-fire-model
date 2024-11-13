//
// Created by lukassid on 15.07.19.
//

#include "neuron_models.h"

bool
qif(double *v, double *t, double gamma, double mu, double vT, double d, double xi, double stepSize) {
    // The qif model the threshold and reset are problematic as they are at infinity.
    // However since the model is solve-able one can:
    // 1. caluclate a threshold vT so that the spike time is shortened less than some error t_error
    // 2. Introduce a threshold and calculate how long it would have taken the
    // det. qif model to reach the actual threshold.
    // For now I go with 1. vT > tan(pi/2 - t_error*sqrt(beta)) with t_error = 10^-5 and beta = 1.
    *t += stepSize;
    if (*v >= vT) {
        *v = -vT;
        return true;
    }
    *v += stepSize * (mu + gamma*pow(*v, 2.)) + sqrt(stepSize * 2 * d) * xi;
    return false;
}


bool
cnoiseQif(double *v, double *n, double *t, double gamma, double mu, double vT, double tauN, double dN, double dW, double xiN, double xiW, double stepSize) {
    // The the qif model the threshold and reset are problematic as they are at infinity.
    // However since the model is solve-able one can:
    // 1. caluclate a threshold vT so that the spike time is shortened less than some error t_error
    // 2. Introduce a threshold and calculate how long it would have taken the
    // det. qif model to reach the actual threshold.
    // For now I go with 1. vT > tan(pi/2 - t_error*sqrt(beta)) with t_error = 10^-5 and beta = 1.
    *t += stepSize;
    if (*v >= vT) {
        *v = -vT;
        return true;
    }
    *v += stepSize * (mu + gamma*pow(*v, 2.) - *n) + sqrt(stepSize * 2 * dW) * xiW;
    *n += stepSize * (-*n / tauN) + sqrt(stepSize * 2 * dN) * (xiN / tauN);
    return false;
}


bool
adapQif(double *v, double *a, double *t, double vT, double gamma, double mu, double tauA, double delta,
        double d, double xi, double stepSize) {
    *t += stepSize;
    if (*v >= vT) {
        *v = -vT;
        *a += delta / tauA;
        return true;
    }
    *v += stepSize * (mu + gamma*pow(*v, 2.) - *a) + sqrt(stepSize * 2 * d) * xi;
    *a += stepSize * (-*a / tauA);
    return false;
}


bool
cnoiseAdapQif(double *v, double *a, double *n, double *t, double vT, double gamma, double mu,
              double tauA, double delta, double tauN, double dN, double dW, double xiN, double xiW, double stepSize) {
    *t += stepSize;
    if (*v >= vT) {
        *v = -vT;
        *a += delta / tauA;
        return true;
    }
    *v += stepSize * (mu + gamma*pow(*v, 2.) - *a - *n) + sqrt(stepSize * 2 * dW) * xiW;
    *a += stepSize * (-*a / tauA);
    *n += stepSize * (-*n / tauN) + sqrt(stepSize * 2 * dN) * xiN / tauN;
    return false;
}


bool
lif(double *v, double *t, double gamma, double mu, double d, double xi, double stepSize) {

    /* Differential equation: lif Model
   * \dot{v} = -γ*v + μ - a + √2D*ξ
   */

    *t += stepSize;
    if (*v >= 1.) {
        *v = 0;
        return true;
    }
    *v += stepSize * (-gamma*(*v) + mu) + sqrt(stepSize * 2 * d) * xi;
    return false;
}


bool
adapLif(double *v, double *a, double *t, double gamma, double mu, double tauA, double delta, double d, double xi, double stepSize) {

    /* Differential equation: Exponential lif Model
   * \dot{v} = -γ*v + μ - a + √2D*ξ
   * \dot{a} = -a/τ_a + Δ*\sum δ(t - t_i)
   */
    *t += stepSize;
    if (*v >= 1.) {
        *v = 0;
        *a += delta / tauA;
        return true;
    }
    *v += stepSize * (-gamma*(*v) + mu - *a) + sqrt(stepSize * 2 * d) * xi;
    *a += stepSize * (-*a / tauA);
    return false;
}


bool
cnoiseLif(double *v, double *n, double *t, double gamma, double mu, double tauN,
          double dN, double dW, double xiN, double xiW, double stepSize) {

    /* Differential equation: lif Model with colored noise
    * \dot{v} = -γ*v + μ + η
    * \dot{η} = -η/τ_η + √2D * ξ/τ_η
    */

    *t += stepSize;
    if (*v >= 1.) {
        *v = 0;
        return true;
    }
    *v += stepSize * (-gamma*(*v) + mu - *n) + sqrt(stepSize * 2 * dW) * xiW;
    *n += stepSize * (-*n / tauN) + sqrt(stepSize * 2 * dN) * (xiN / tauN);
    return false;
}


bool
cnoiseAdapLif(double *v, double *a, double *n, double *t, double gamma, double mu, double tauA, double delta, double tauN, double dN, double dW, double xiN, double xiW, double stepSize) {

    /* Differential equation: Exponential lif Model
    * \dot{v} = -γ*v + μ - a + η
    * \dot{η} = -η + √2D*ξ
    * \dot{a} = -a/τ_a + Δ*\sum δ(t - t_i)
    */
    *t += stepSize;
    if (*v >= 1.) {
        *v = 0;
        *a += delta / tauA;
        return true;
    }
    *v += stepSize * (-gamma*(*v) + mu - *a - *n) + sqrt(stepSize * 2 * dW) * xiW;
    *n += stepSize * (-*n / tauN) + sqrt(stepSize * 2 * dN) * xiN / tauN;
    *a += stepSize * (-*a / tauA);
    return false;
}


bool
gif(double *v, double *w, double *t, double gamma, double mu, double beta, double tauW, double wReset,
    double d, double xi, double stepSize) {

    /* Differential equation: Exponential lif Model
    * \dot{v} = -γ*v -β*w + μ + √2D*ξ
    * \dot{w} = (v-w)/τ_w
    */

    *t += stepSize;
    if (*v >= 1.) {
        *v = 0;
        *w = wReset;
        return true;
    }
    double vTmp = *v;
    double wTmp = *w;

    *v += stepSize * (-gamma*(vTmp) - beta * (wTmp) + mu) + sqrt(stepSize * 2 * d) * xi;
    *w += stepSize * (vTmp - wTmp) / tauW;
    return false;
}


bool
adapGif(double *v, double *w, double *a, double *t, double gamma, double mu, double beta, double tauW, double wReset, double tauA,
        double delta, double d, double xi, double stepSize) {


    /* Differential equation: Exponential lif Model
    * \dot{v} = -γ*v -β*w + μ - a + √2D*ξ
    * \dot{w} = (v-w)/τ_w
    * \dot{a} = -a/τ_a + Δ*\sum δ(t - t_i)
    */

    *t += stepSize;
    if (*v >= 1.) {
        *v = 0;
        *w = wReset;
        *a += delta / tauA;
        return true;
    }
    double vTmp = *v;
    double wTmp = *w;

    *v += stepSize * (-gamma*(vTmp) - beta * (wTmp) + mu - *a) + sqrt(stepSize * 2 * d) * xi;
    *w += stepSize * (vTmp - wTmp) / tauW;
    *a += stepSize * (-*a / tauA);
    return false;
}


bool
cnoiseGif(double *v, double *w, double *n, double *t, double gamma, double mu, double beta, double tauW, double wReset, double tauN, double dN, double dW, double xiN, double xiW, double stepSize) {

    /* Differential equation: Exponential lif Model
   * \dot{v} = -γ*v -β*w + μ + √2D*ξ
   * \dot{w} = (v-w)/τ_w
   */

    *t += stepSize;
    if (*v >= 1.) {
        *v = 0;
        *w = wReset;
        return true;
    }
    double vTmp = *v;
    double wTmp = *w;

    *v += stepSize * (-gamma*(vTmp) - beta * (wTmp) + mu - *n) + sqrt(stepSize * 2 * dW) * xiW;
    *w += stepSize * (vTmp - wTmp) / tauW;
    *n += stepSize * (-*n / tauN) + sqrt(stepSize * 2 * dN) * xiN / tauN;
    return false;
}


bool
cnoiseAdapGif(double *v, double *w, double *a, double *n, double *t, double gamma, double mu, double beta, double tauW, double wReset, double tauA, double delta, double tauN, double dN, double dW, double xiN, double xiW, double stepSize) {


    /* Differential equation: Exponential lif Model
    * \dot{v} = -γ*v -β*w + μ + η
    * \dot{w} = (v-w)/τ_w
    * \dot{η} = -η + √2D*ξ
    * \dot{a} = -a/τ_a + Δ*\sum δ(t - t_i)
    */

    *t += stepSize;
    if (*v >= 1.) {
        *v = 0;
        *w = wReset;
        *a += delta / tauA;
        return true;
    }
    double vTmp = *v;
    double wTmp = *w;

    *v += stepSize * (-gamma*(vTmp) - beta * (wTmp) + mu - (*a) - *n) + sqrt(stepSize * 2 * dW) * xiW;
    *w += stepSize * (vTmp - wTmp) / tauW;
    *n += stepSize * (-*n / tauN) + sqrt(stepSize * 2 * dN) * xiN / tauN;
    *a += stepSize * (-*a / tauA);
    return false;
}