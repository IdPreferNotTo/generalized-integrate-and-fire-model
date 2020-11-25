#include <iostream>
#include <iomanip>
#include <random>
#include <unistd.h>
#include <pwd.h>
#include <cmath>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include "neuron_models.h"
#include "linear_response.h"

using namespace std;
namespace pt = boost::property_tree;

tuple<double, double, double>
get_deterministic_values(const string *model, double *v, double *w, double *a, double *t, double gamma, double input, double beta_gif, double tau_gif, double gif_reset,
                         double tau_a, double jump_a, double dt) {
    double T = 0;
    double t_tmp = *t;
    double w_tmp = *w;
    double dif_t = 10;
    double precision_delta_t = 0.00001;
    int spikes = 0;
    bool fire = false;

    struct passwd *pw = getpwuid(getuid());
    const char *homedir = pw->pw_dir;

    string path = "/Data/" + *model + "/data/";
    char parameters[100];
    if((*model=="LIF") || (*model=="QIF")){
        std::sprintf(parameters, "mu%.2f_taua%.1f_Delta%.1f.txt", input,
                     tau_a, jump_a);
    }
    if(*model=="GIF"){
        std::sprintf(parameters, "mu%.2f_beta%.1f_tauw%.1f_taua%.1f_Delta%.1f.txt", input, beta_gif, tau_gif,
                     tau_a, jump_a);
    }

    string dataFile;
    dataFile = string(homedir) + path + parameters;
    ofstream file_a;
    file_a.open(dataFile);
    int N = 0;
    while (dif_t > precision_delta_t) {
        N +=1;
        if (spikes < 1000) {
            if (tau_a > 0.001) {
                w_tmp = *w;
                if (*model == "LIF") {
                    fire = adap_LIF(v, a, t, gamma, input, tau_a, jump_a, 0, 0, dt);
                }
                if (*model == "GIF") {
                    fire = adap_GIF(v, w, a, t, gamma, input, beta_gif, tau_gif, gif_reset, tau_a, jump_a, 0, 0, dt);
                }
                if (*model == "QIF") {
                    fire = adap_QIF(v, a, t, 10000, gamma,  input, tau_a, jump_a, 0, 0, dt);
                }
                if (N%100==0){
                    file_a << *t << ' ' << *a << ' ' << "\n";
                }
            } else {
                w_tmp = *w;
                //fire = theta(v, t, 0, stepSize, 0.);
                if (*model == "LIF") {
                    fire = LIF(v, t, gamma,  input, 0, 0, dt);
                }
                if (*model == "GIF") {
                    fire = GIF(v, w, t, gamma, input, beta_gif, tau_gif, gif_reset, 0, 0, dt);
                }
                if (*model == "QIF") {
                    fire = QIF(v, t, gamma, input, 10000, 0, 0, dt);
                }
            }
            if (fire) {
                dif_t = abs(T - (*t - t_tmp));
                T = *t - t_tmp;
                t_tmp = *t;
                spikes++;
            }
        } else {
            spikes = 0;
            dt /= 2;
        }
    }
    file_a.close();
    return make_tuple(T, w_tmp, *a);
}

void
save_noise(size_t bufferSize, vector<vector<double>> &x, vector<int> &stepCounts, int spikeCount, int nonLinRespCount,
           double noise) {

    if (spikeCount >= nonLinRespCount) {
        if (spikeCount % bufferSize < nonLinRespCount % bufferSize) {
            spikeCount %= bufferSize;
            nonLinRespCount %= bufferSize;
            for (int i = nonLinRespCount; i < bufferSize; i++) {
                x[i][stepCounts[i]] = noise;
                stepCounts[i] += 1;
            }
            for (int i = 0; i < spikeCount + 1; i++) {
                x[i][stepCounts[i]] = noise;
                stepCounts[i] += 1;
            }
        } else {
            spikeCount %= bufferSize;
            nonLinRespCount %= bufferSize;
            for (int i = nonLinRespCount; i < spikeCount + 1; i++) {
                x[i][stepCounts[i]] = noise;
                stepCounts[i] += 1;
            }
        }
    }
}


int main(int argc, char *argv[]) {
    //---------- read parameters form json file with boost -------------------------------------------------------------
    struct passwd *pw = getpwuid(getuid());
    const char *homedir = pw->pw_dir;

    string paraFile = string(homedir) + "/Parameter/" + string(argv[1]) + ".json";
    pt::ptree desc;
    pt::json_parser::read_json(paraFile, desc);

    // General parameters
    const auto model = desc.get<string>("neuron model.model");
    const auto input = desc.get<double>("neuron model.mu");
    const auto gamma = desc.get<double>("neuron model.gamma");
    // GIF parameters
    const auto beta_gif = desc.get<double>("neuron model.GIF.beta");
    const auto tau_gif = desc.get<double>("neuron model.GIF.time constant");
    const auto reset_gif = desc.get<double>("neuron model.GIF.reset");

    // Adaptation parameters
    const auto tau_a = desc.get<double>("adaptation.time constant");
    const auto jump_a = desc.get<double>("adaptation.strength");

    // Noise parameters
    const auto D_wn = desc.get<double>("noise.intensity");

    // OUP parameters
    const auto D_ou = desc.get<double>("OU.intensity");
    const auto tau_ou = desc.get<double>("OU.time constant");

    // Numerical parameters
    auto dt = desc.get<double>("num parameter.step size");
    const auto maxSpikeCount = desc.get<int>("num parameter.max spikes");
    const auto N = desc.get<int>("num parameter.run");

    // ------------------------ Parameters for output file -------------------------------------------------------------
    string path = "/Data/" + model + "/data/";
    char parameters[100];
    if((model=="LIF") || (model=="QIF")){
        std::sprintf(parameters, "mu%.2f_taua%.1f_Delta%.1f_taun%.3f_Dn%.2e_Dw%.2e_%d.txt", input,
                 tau_a, jump_a, tau_ou, D_ou, D_wn, N);
    }
    if(model=="GIF"){
        std::sprintf(parameters, "mu%.2f_beta%.1f_tauw%.1f_taua%.1f_Delta%.1f_taun%.3f_Dn%.2e_Dw%.2e.txt", input, beta_gif, tau_gif,
                     tau_a, jump_a, tau_ou, D_ou, D_wn);
    }
    string dataFile;
    dataFile = string(homedir) + path + parameters;
    ofstream file;
    file.open(dataFile);
    if (!file.is_open()) {
        cout << "Could not open file at: " << dataFile << endl;
        return 1;
    }
    file << "# t delta_a eta lin.resp" << endl;

    // -------------------------- Noise parameters & random number generator -------------------------------------------
    // This is how random number should be generated. Best practice!
    double rng_ou = 0;
    double rng_wn = 0;
    const double mean = 0.0;
    const double stddev = 1.0;
    std::random_device rd;
    std::mt19937 generator(rd());
    //Better seed from random_device instead of clock in case one runs many simulations in a short periode of time
    std::normal_distribution<double> dist(mean, stddev);
    vector<vector<double>> dataBuffer(maxSpikeCount, vector<double>(4));
    // dataBuffer = [[t1, δa_1, χ_1, η_1], [t2, δa_2, χ_2, η_2], ...]


    // -------------------------- Get deterministic values -------------------------------------------------------------
    double v = 0; // v
    double w = 0;
    double noise = 0.; // η
    double a = jump_a; // a
    double t = 0; // t
    double tDet, wDet, aDet;
    tie(tDet, wDet, aDet) = get_deterministic_values(&model, &v, &w, &a, &t, gamma, input, beta_gif, tau_gif, reset_gif, tau_a, jump_a, dt);
    file << "# tDet: " << tDet << '\n' << "# wDet: " << wDet << '\n' << "# aDet: " << aDet << endl;
    file << "# dt: " << dt << endl;

    // ------------------------- Simulate neuron dynamics --------------------------------------------------------------
    // Create a buffer that saves noise values between two subsequent spikes. This becomes necessary if one wants to
    // calculate a non-linear response for which the instantaneous noise is not sufficient due to the occuring
    // double integral over \eta(t_i + t2)*\eta(t_i + t1).
    size_t bufferSize = 50;
    int tDetSteps = int(tDet / dt);
    vector<vector<double>> noiseBuffer(bufferSize, vector<double>(tDetSteps));
    vector<int> noiseStepCounter(100, 0);
    int noiseBufferCount = 0;

    t = 0; // Reset time. Other values i.e. v, a, (w) remain unchanged and thus are v=0, a=a_det, w=w_det
    double tFire = 0;
    int spikeCount = 0;
    bool neuronFired = false;

    const bool isAdaptiv = (tau_a > 0.001);
    const bool isColored = (tau_ou > 0.001);
    while (spikeCount < maxSpikeCount) {
        rng_ou = dist(generator);
        rng_wn = dist(generator);
        /*
         * depending on tau_a and tau_n different model models have to be used as tau -> 0 implies a division
         * by 0 which can not be handled.
         */
        if (isColored) {
            if (isAdaptiv) {
                if (model == "LIF") {
                    neuronFired = cnoise_adap_LIF(&v, &a, &noise, &t, gamma, input, tau_a, jump_a, tau_ou, D_ou, D_wn, rng_ou, rng_wn, dt);
                    //neuronFired = cnoise_extra_noisy_adap_LIF(&v, &a, &noise, &t, rng_ou, rng_wn, dt, gamma, input, tau_a, jump_a, tau_ou, D_ou, D_wn);
                }
                if (model == "GIF") {
                    neuronFired = cnoise_adap_GIF(&v, &w, &a, &noise, &t, gamma, input, beta_gif, tau_gif, reset_gif, tau_a, jump_a, tau_ou, D_ou, D_wn, rng_ou, rng_wn, dt);
                }
                if (model == "QIF") {
                    neuronFired = cnoise_adap_QIF(&v, &a, &noise, &t, 10000, gamma, input, tau_a, jump_a, tau_ou, D_ou, D_wn, rng_ou, rng_wn, dt);
                }
            } else {
                if (model == "LIF") {
                    neuronFired = cnoise_LIF(&v, &noise, &t, gamma, input, tau_ou, D_ou, D_wn, rng_ou, rng_wn, dt);
                }
                if (model == "GIF") {
                    neuronFired = cnoise_GIF(&v, &w, &noise, &t, gamma, input, beta_gif, tau_gif, reset_gif, tau_ou, D_ou, D_wn, rng_ou, rng_wn, dt);
                }
                if (model == "QIF") {
                    neuronFired = cnoise_QIF(&v, &noise, &t, 10000, gamma, input, tau_ou, D_ou, D_wn, rng_ou, rng_wn, dt);
                }
            }
        } else {
            if (isAdaptiv) {
                if (model == "LIF") {
                    neuronFired = adap_LIF(&v, &a, &t, gamma, input, tau_a, jump_a, D_wn, rng_wn, dt);
                }
                if (model == "GIF") {
                    neuronFired = adap_GIF(&v, &w, &a, &t, gamma, input, beta_gif, tau_gif, reset_gif, tau_a, jump_a, D_wn, rng_wn, dt);
                }
                if (model == "QIF") {
                    neuronFired = adap_QIF(&v, &a, &t, 10000, gamma, input, tau_a, jump_a, D_wn, rng_wn, dt);
                }
                noise = rng_wn * sqrt(dt * 2 * D_wn);
            } else {
                if (model == "LIF") {
                    neuronFired = LIF(&v, &t, gamma, input, D_wn, rng_wn, dt);
                }
                if (model == "GIF") {
                    neuronFired = GIF(&v, &w, &t, gamma, input, beta_gif, tau_gif, reset_gif, D_wn, rng_wn, dt);
                }
                if (model == "QIF") {
                    neuronFired = QIF(&v, &t, gamma, input, 10000, D_wn, rng_wn, dt);
                }
                noise = rng_wn * sqrt(dt * 2 * D_wn);
            }
        }

        if (neuronFired) {
            if (spikeCount % 100 == 0) {
                cout << spikeCount << endl;
            }
            dataBuffer[spikeCount][0] = t - tFire;
            dataBuffer[spikeCount][1] = a - aDet;
            dataBuffer[spikeCount][2] = noise;
            spikeCount++;
            tFire = t;
        }

        save_noise(bufferSize, noiseBuffer, noiseStepCounter, spikeCount, noiseBufferCount, noise);
        if (noiseStepCounter[noiseBufferCount % bufferSize] >= tDetSteps) {
            double linResp = 0;
            if (model == "LIF") {
                linResp = calculate_lif_lin_resp(noiseBuffer[noiseBufferCount % bufferSize], dt);
            }
            if (model == "GIF") {
                linResp = calculate_gif_lin_resp(noiseBuffer[noiseBufferCount % bufferSize], dt, input, beta_gif, tau_gif, jump_a, tDet, aDet, wDet);
            }
            if (model == "QIF") {
                linResp = calculate_qif_lin_resp(noiseBuffer[noiseBufferCount % bufferSize], dt, input);
            }
            dataBuffer[noiseBufferCount][3] = linResp;
            noiseStepCounter[noiseBufferCount % bufferSize] = 0;
            noiseBufferCount++;
        }
    }

    for (const auto &data: dataBuffer) {
        for (const auto &x: data) {
            file << std::setprecision(8) << x << ' ';
        }
        file << '\n';
    }

    file.close();
    return 0;
}
