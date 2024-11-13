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
getDeterministicValues(const string *model, double *v, double *w, double *a, double *t, double gamma, double input, double betaGif, double tauGif, double gifReset,
                       double tauA, double deltaA, double dt) {
    double tInit = 0;
    double tTmp = *t;
    double wTmp = *w;
    double difT = 10;
    double tStepSize = 0.00001;
    int spikeNumber = 0;
    bool hasFired = false;

    struct passwd *pw = getpwuid(getuid());
    const char *homedir = pw->pw_dir;

    string path = "/Data/" + *model + "/data/";
    char parameters[100];
    if((*model=="LIF") || (*model=="QIF")){
        std::sprintf(parameters, "mu%.2f_taua%.1f_Delta%.1f.txt", input,
                     tauA, deltaA);
    }
    if(*model=="GIF"){
        std::sprintf(parameters, "mu%.2f_beta%.1f_tauw%.1f_taua%.1f_Delta%.1f.txt", input, betaGif, tauGif,
                     tauA, deltaA);
    }

    string dataFile;
    dataFile = string(homedir) + path + parameters;
    ofstream fileA;
    fileA.open(dataFile);
    int n = 0;
    while (difT > tStepSize) {
        n +=1;
        if (spikeNumber < 1000) {
            if (tauA > 0.001) {
                wTmp = *w;
                if (*model == "LIF") {
                    hasFired = adapLif(v, a, t, gamma, input, tauA, deltaA, 0, 0, dt);
                }
                if (*model == "GIF") {
                    hasFired = adapGif(v, w, a, t, gamma, input, betaGif, tauGif, gifReset, tauA, deltaA, 0, 0, dt);
                }
                if (*model == "QIF") {
                    hasFired = adapQif(v, a, t, 10000, gamma, input, tauA, deltaA, 0, 0, dt);
                }
                if (n % 100 == 0){
                    fileA << *t << ' ' << *a << ' ' << "\n";
                }
            } else {
                wTmp = *w;
                //hasFired = theta(v, t, 0, stepSize, 0.);
                if (*model == "LIF") {
                    hasFired = lif(v, t, gamma, input, 0, 0, dt);
                }
                if (*model == "GIF") {
                    hasFired = gif(v, w, t, gamma, input, betaGif, tauGif, gifReset, 0, 0, dt);
                }
                if (*model == "QIF") {
                    hasFired = qif(v, t, gamma, input, 10000, 0, 0, dt);
                }
            }
            if (hasFired) {
                difT = abs(tInit - (*t - tTmp));
                tInit = *t - tTmp;
                tTmp = *t;
                spikeNumber++;
            }
        } else {
            spikeNumber = 0;
            dt /= 2;
        }
    }
    fileA.close();
    return make_tuple(tInit, wTmp, *a);
}

void
saveNoise(size_t bufferSize, vector<vector<double>> &x, vector<int> &stepCounts, int spikeCount, int nonLinRespCount,
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
    // const char *homedir = pw->pw_dir;

    string paraFile = "../parameter/" + string(argv[1]) + ".json";
    pt::ptree desc;
    pt::json_parser::read_json(paraFile, desc);

    // General parameters
    const auto modelName = desc.get<string>("neuron model.model");
    const auto input = desc.get<double>("neuron model.mu");
    const auto gamma = desc.get<double>("neuron model.gamma");
    // gif parameters
    const auto betaGif = desc.get<double>("neuron model.GIF.beta");
    const auto tauGif = desc.get<double>("neuron model.GIF.time constant");
    const auto resetGif = desc.get<double>("neuron model.GIF.reset");

    // Adaptation parameters
    const auto tauA = desc.get<double>("adaptation.time constant");
    const auto deltaA = desc.get<double>("adaptation.strength");

    // Noise parameters
    const auto dWn = desc.get<double>("noise.intensity");

    // OUP parameters
    const auto dOu = desc.get<double>("OU.intensity");
    const auto tauOu = desc.get<double>("OU.time constant");

    // Numerical parameters
    auto dt = desc.get<double>("num parameter.step size");
    const auto maxSpikeCount = desc.get<int>("num parameter.max spikes");
    const auto nrRun = desc.get<int>("num parameter.run");

    // ------------------------ Parameters for output file -------------------------------------------------------------
    std::string path = "../out/";
    char parameters[100];
    if((modelName == "LIF") || (modelName == "QIF")){
        std::sprintf(parameters, "mu%.2f_taua%.1f_Delta%.1f_taun%.3f_Dn%.2e_Dw%.2e_%d.txt", input,
                     tauA, deltaA, tauOu, dOu, dWn, nrRun);
    }
    if(modelName == "GIF"){
        std::sprintf(parameters, "mu%.2f_beta%.1f_tauw%.1f_taua%.1f_Delta%.1f_taun%.3f_Dn%.2e_Dw%.2e.txt", input, betaGif, tauGif,
                     tauA, deltaA, tauOu, dOu, dWn);
    }
    string dataFile;
    dataFile  = path + parameters;
    ofstream file;
    file.open(dataFile);
    if (!file.is_open()) {
        cout << "Could not open file at: " << dataFile << endl;
        return 1;
    }
    file << "# t deltaA eta lin.resp" << endl;

    // -------------------------- Noise parameters & random number generator -------------------------------------------
    // This is how random number should be generated. Best practice!
    double rngOu;
    double rngWn;
    const double mean = 0.0;
    const double stdDev = 1.0;
    std::random_device rd;
    std::mt19937 generator(rd());
    //Better seed from random_device instead of clock in case one runs many simulations in a short periode of time
    std::normal_distribution<double> dist(mean, stdDev);
    vector<vector<double>> dataBuffer(maxSpikeCount, vector<double>(4));
    // dataBuffer = [[t1, δa_1, χ_1, η_1], [t2, δa_2, χ_2, η_2], ...]


    // -------------------------- Get deterministic values -------------------------------------------------------------
    double v = 0;
    double w = 0;
    double noise = 0.;
    double a = deltaA;
    double t = 0;
    double tDet, wDet, aDet;
    tie(tDet, wDet, aDet) = getDeterministicValues(&modelName, &v, &w, &a, &t, gamma, input, betaGif, tauGif, resetGif,
                                                   tauA, deltaA, dt);
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

    t = 0; // Reset time. Other values i.e. v, a, (w) remain unchanged and thus are v=0, a=aDet, w=wDet
    double tFire = 0;
    int spikeCount = 0;
    bool neuronFired = false;

    const bool isAdaptive = (tauA > 0.001);
    const bool isColored = (tauOu > 0.001);
    while (spikeCount < maxSpikeCount) {
        rngOu = dist(generator);
        rngWn = dist(generator);
        /*
         * depending on tauA and tau_n different modelName models have to be used as tau -> 0 implies a division
         * by 0 which can not be handled.
         */
        if (isColored) {
            if (isAdaptive) {
                if (modelName == "LIF") {
                    neuronFired = cnoiseAdapLif(&v, &a, &noise, &t, gamma, input, tauA, deltaA, tauOu, dOu, dWn, rngOu,
                                                rngWn, dt);
                    //neuronFired = cnoise_extra_noisy_adap_LIF(&v, &a, &noise, &t, rngOu, rngWn, dt, gamma, input, tauA, deltaA, tauOu, dOu, dWn);
                }
                if (modelName == "GIF") {
                    neuronFired = cnoiseAdapGif(&v, &w, &a, &noise, &t, gamma, input, betaGif, tauGif, resetGif,
                                                tauA, deltaA, tauOu, dOu, dWn, rngOu, rngWn, dt);
                }
                if (modelName == "QIF") {
                    neuronFired = cnoiseAdapQif(&v, &a, &noise, &t, 10000, gamma, input, tauA, deltaA, tauOu, dOu, dWn,
                                                rngOu, rngWn, dt);
                }
            } else {
                if (modelName == "LIF") {
                    neuronFired = cnoiseLif(&v, &noise, &t, gamma, input, tauOu, dOu, dWn, rngOu, rngWn, dt);
                }
                if (modelName == "GIF") {
                    neuronFired = cnoiseGif(&v, &w, &noise, &t, gamma, input, betaGif, tauGif, resetGif, tauOu, dOu,
                                            dWn, rngOu, rngWn, dt);
                }
                if (modelName == "QIF") {
                    neuronFired = cnoiseQif(&v, &noise, &t, gamma, input, 10000, tauOu, dOu, dWn, rngOu, rngWn, dt);
                }
            }
        } else {
            if (isAdaptive) {
                if (modelName == "LIF") {
                    neuronFired = adapLif(&v, &a, &t, gamma, input, tauA, deltaA, dWn, rngWn, dt);
                }
                if (modelName == "GIF") {
                    neuronFired = adapGif(&v, &w, &a, &t, gamma, input, betaGif, tauGif, resetGif, tauA, deltaA, dWn,
                                          rngWn, dt);
                }
                if (modelName == "QIF") {
                    neuronFired = adapQif(&v, &a, &t, 10000, gamma, input, tauA, deltaA, dWn, rngWn, dt);
                }
                noise = rngWn * sqrt(dt * 2 * dWn);
            } else {
                if (modelName == "LIF") {
                    neuronFired = lif(&v, &t, gamma, input, dWn, rngWn, dt);
                }
                if (modelName == "GIF") {
                    neuronFired = gif(&v, &w, &t, gamma, input, betaGif, tauGif, resetGif, dWn, rngWn, dt);
                }
                if (modelName == "QIF") {
                    neuronFired = qif(&v, &t, gamma, input, 10000, dWn, rngWn, dt);
                }
                noise = rngWn * sqrt(dt * 2 * dWn);
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

        saveNoise(bufferSize, noiseBuffer, noiseStepCounter, spikeCount, noiseBufferCount, noise);
        if (noiseStepCounter[noiseBufferCount % bufferSize] >= tDetSteps) {
            double linResp = 0;
            if (modelName == "lif") {
                linResp = calculateLifLinResp(noiseBuffer[noiseBufferCount % bufferSize], dt);
            }
            if (modelName == "gif") {
                linResp = calculateGifLinResp(noiseBuffer[noiseBufferCount % bufferSize], dt, input, betaGif, tauGif,
                                              deltaA, tDet, aDet, wDet);
            }
            if (modelName == "qif") {
                linResp = calculateQifLinResp(noiseBuffer[noiseBufferCount % bufferSize], dt, input);
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
