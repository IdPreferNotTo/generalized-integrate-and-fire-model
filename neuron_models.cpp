//
// Created by lukassid on 15.07.19.
//

#include "neuron_models.h"

bool
QIF(double *v, double *t, double gamma, double mu, double vT, double D, double xi, double stepSize) {
    // The the QIF model the threshold and reset are problematic as they are at infinity.
    // However since the model is solve-able one can:
    // 1. caluclate a threshold vT so that the spike time is shortened less than some error t_error
    // 2. Introduce a threshold and calculate how long it would have taken the
    // det. QIF model to reach the actual threshold.
    // For now I go with 1. vT > tan(pi/2 - t_error*sqrt(beta)) with t_error = 10^-5 and beta = 1.
    *t += stepSize;
    if (*v >= vT) {
        *v = -vT;
        return true;
    }
    *v += stepSize * (mu + gamma*pow(*v, 2.)) + sqrt(stepSize * 2 * D) * xi;
    return false;
}


bool
cnoise_QIF(double *v, double *n, double *t, double gamma, double mu, double vT, double tau_n, double D_n, double D_w, double xi_n, double xi_w, double stepSize) {
    // The the QIF model the threshold and reset are problematic as they are at infinity.
    // However since the model is solve-able one can:
    // 1. caluclate a threshold vT so that the spike time is shortened less than some error t_error
    // 2. Introduce a threshold and calculate how long it would have taken the
    // det. QIF model to reach the actual threshold.
    // For now I go with 1. vT > tan(pi/2 - t_error*sqrt(beta)) with t_error = 10^-5 and beta = 1.
    *t += stepSize;
    if (*v >= vT) {
        *v = -vT;
        return true;
    }
    *v += stepSize * (mu + gamma*pow(*v, 2.) - *n) + sqrt(stepSize * 2 * D_w) * xi_w;
    *n += stepSize * (-*n / tau_n) + sqrt(stepSize * 2 * D_n) * (xi_n / tau_n);
    return false;
}


bool
adap_QIF(double *v, double *a, double *t, double vT, double gamma, double mu, double tau_a, double delta,
         double D, double xi, double stepSize) {
    *t += stepSize;
    if (*v >= vT) {
        *v = -vT;
        *a += delta/tau_a;
        return true;
    }
    *v += stepSize * (mu + gamma*pow(*v, 2.) - *a) + sqrt(stepSize * 2 * D) * xi;
    *a += stepSize * (-*a / tau_a);
    return false;
}


bool
cnoise_adap_QIF(double *v, double *a, double *n, double *t, double vT, double gamma, double mu,
                double tau_a, double delta, double tau_n, double D_n, double D_w, double xi_n, double xi_w, double stepSize) {
    *t += stepSize;
    if (*v >= vT) {
        *v = -vT;
        *a += delta/tau_a;
        return true;
    }
    *v += stepSize * (mu + gamma*pow(*v, 2.) - *a - *n) + sqrt(stepSize * 2 * D_w) * xi_w;
    *a += stepSize * (-*a / tau_a);
    *n += stepSize * (-*n / tau_n) + sqrt(stepSize * 2 * D_n) * xi_n / tau_n;
    return false;
}


bool
LIF(double *v, double *t, double gamma, double mu, double D, double xi, double stepSize) {

    /* Differential equation: LIF Model
   * \dot{v} = -γ*v + μ - a + √2D*ξ
   */

    *t += stepSize;
    if (*v >= 1.) {
        *v = 0;
        return true;
    }
    *v += stepSize * (-gamma*(*v) + mu) + sqrt(stepSize * 2 * D) * xi;
    return false;
}


bool
adap_LIF(double *v, double *a, double *t, double gamma,  double mu, double tau_a, double delta, double D, double xi, double stepSize) {

    /* Differential equation: Exponential LIF Model
   * \dot{v} = -γ*v + μ - a + √2D*ξ
   * \dot{a} = -a/τ_a + Δ*\sum δ(t - t_i)
   */
    *t += stepSize;
    if (*v >= 1.) {
        *v = 0;
        *a += delta/tau_a;
        return true;
    }
    *v += stepSize * (-gamma*(*v) + mu - *a) + sqrt(stepSize * 2 * D) * xi;
    *a += stepSize * (-*a / tau_a);
    return false;
}


bool
cnoise_LIF(double *v, double *n, double *t, double gamma, double mu, double tau_n,
           double D_n, double D_w, double xi_n, double xi_w, double stepSize) {

    /* Differential equation: LIF Model with colored noise
    * \dot{v} = -γ*v + μ + η
    * \dot{η} = -η/τ_η + √2D * ξ/τ_η
    */

    *t += stepSize;
    if (*v >= 1.) {
        *v = 0;
        return true;
    }
    *v += stepSize * (-gamma*(*v) + mu - *n) + sqrt(stepSize * 2 * D_w) * xi_w;
    *n += stepSize * (-*n / tau_n) + sqrt(stepSize * 2 * D_n) * (xi_n / tau_n);
    return false;
}


bool
cnoise_adap_LIF(double *v, double *a, double *n, double *t, double gamma, double mu, double tau_a, double delta, double tau_n, double D_n, double D_w, double xi_n, double xi_w, double stepSize) {

    /* Differential equation: Exponential LIF Model
    * \dot{v} = -γ*v + μ - a + η
    * \dot{η} = -η + √2D*ξ
    * \dot{a} = -a/τ_a + Δ*\sum δ(t - t_i)
    */
    *t += stepSize;
    if (*v >= 1.) {
        *v = 0;
        *a += delta/tau_a;
        return true;
    }
    *v += stepSize * (-gamma*(*v) + mu - *a - *n) + sqrt(stepSize * 2 * D_w) * xi_w;
    *n += stepSize * (-*n / tau_n) + sqrt(stepSize * 2 * D_n) * xi_n / tau_n;
    *a += stepSize * (-*a / tau_a);
    return false;
}


bool
GIF(double *v, double *w, double *t, double gamma, double mu, double beta, double tau_w, double w_reset,
    double D, double xi, double stepSize) {

    /* Differential equation: Exponential LIF Model
    * \dot{v} = -γ*v -β*w + μ + √2D*ξ
    * \dot{w} = (v-w)/τ_w
    */

    *t += stepSize;
    if (*v >= 1.) {
        *v = 0;
        *w = w_reset;
        return true;
    }
    double v_tmp = *v;
    double w_tmp = *w;

    *v += stepSize * (-gamma*(v_tmp) - beta * (w_tmp) + mu) + sqrt(stepSize * 2 * D) * xi;
    *w += stepSize * (v_tmp - w_tmp) / tau_w;
    return false;
}


bool
adap_GIF(double *v, double *w, double *a, double *t, double gamma, double mu, double beta, double tau_w, double w_reset, double tau_a,
         double delta, double D, double xi, double stepSize) {


    /* Differential equation: Exponential LIF Model
    * \dot{v} = -γ*v -β*w + μ - a + √2D*ξ
    * \dot{w} = (v-w)/τ_w
    * \dot{a} = -a/τ_a + Δ*\sum δ(t - t_i)
    */

    *t += stepSize;
    if (*v >= 1.) {
        *v = 0;
        *w = w_reset;
        *a += delta/tau_a;
        return true;
    }
    double v_tmp = *v;
    double w_tmp = *w;

    *v += stepSize * (-gamma*(v_tmp) - beta * (w_tmp) + mu - *a) + sqrt(stepSize * 2 * D) * xi;
    *w += stepSize * (v_tmp - w_tmp) / tau_w;
    *a += stepSize * (-*a / tau_a);
    return false;
}


bool
cnoise_GIF(double *v, double *w, double *n, double *t, double gamma, double mu, double beta, double tau_w, double w_reset, double tau_n, double D_n, double D_w, double xi_n, double xi_w, double stepSize) {

    /* Differential equation: Exponential LIF Model
   * \dot{v} = -γ*v -β*w + μ + √2D*ξ
   * \dot{w} = (v-w)/τ_w
   */

    *t += stepSize;
    if (*v >= 1.) {
        *v = 0;
        *w = w_reset;
        return true;
    }
    double v_tmp = *v;
    double w_tmp = *w;

    *v += stepSize * (-gamma*(v_tmp) - beta * (w_tmp) + mu - *n) + sqrt(stepSize * 2 * D_w) * xi_w;
    *w += stepSize * (v_tmp - w_tmp) / tau_w;
    *n += stepSize * (-*n / tau_n) + sqrt(stepSize * 2 * D_n) * xi_n / tau_n;
    return false;
}


bool
cnoise_adap_GIF(double *v, double *w, double *a, double *n, double *t, double gamma,  double mu, double beta, double tau_w, double w_reset, double tau_a, double delta, double tau_n, double D_n, double D_w, double xi_n, double xi_w, double stepSize) {


    /* Differential equation: Exponential LIF Model
    * \dot{v} = -γ*v -β*w + μ + η
    * \dot{w} = (v-w)/τ_w
    * \dot{η} = -η + √2D*ξ
    * \dot{a} = -a/τ_a + Δ*\sum δ(t - t_i)
    */

    *t += stepSize;
    if (*v >= 1.) {
        *v = 0;
        *w = w_reset;
        *a += delta/tau_a;
        return true;
    }
    double v_tmp = *v;
    double w_tmp = *w;

    *v += stepSize * (-gamma*(v_tmp) - beta * (w_tmp) + mu - (*a) - *n) + sqrt(stepSize * 2 * D_w) * xi_w;
    *w += stepSize * (v_tmp - w_tmp) / tau_w;
    *n += stepSize * (-*n / tau_n) + sqrt(stepSize * 2 * D_n) * xi_n / tau_n;
    *a += stepSize * (-*a / tau_a);
    return false;
}