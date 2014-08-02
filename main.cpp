#include <iostream>
#include <stdexcept>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include "fouriere.h"

template<class res> void generate(unsigned const N, double const magnitude, double const frequency,
        double const phase, double const dt, double const lvl, res* samples) {
    double const dphi = 2.0*std::atan(1.0)*4.0*frequency*dt;
    res* const end = samples + N;
    for (double phi = phase; samples < end; ++samples, phi += dphi) {
        *samples = lvl + magnitude*std::sin(phi);
    }
} //generate()

void messageDev(std::ostream& os, unsigned const N, unsigned const m,
        unsigned const k, double const dev, std::string const& hint) {
    os << "Maximum deviation " << dev << " of " << hint << " is for N = " << N
            << ", k = " << k;
    if (m < N) {
        os <<", m = " << m;
    }
    os << "\n";
            ;
} //messageDev()

double rnd(double const from, double const to) {
    return from + std::rand() / static_cast<double>(RAND_MAX) * (to - from);
} //rnd()

float checkDFT100();
void calculateDFT100(unsigned k, float const* const samples, float* const re, float* const im);

float calculateFrequency(float const* const samples) {
    static float const ALFA = 2.0f;
    static float const DF = 6.0f;
    float re, im, dft[5][2];
    float magnitude = 0.0f, lev = 0.0f;
    for (unsigned k = 0; k < 100; ++k) {
        magnitude += samples[k]*samples[k];
        lev += samples[k];
    }
    for (unsigned k = 0; k < 5; ++k) {
        calculateDFT100(k + 3, samples, &re, &im);
        dft[k][0] = re;
        dft[k][1] = im;
    }
    dft[0][1] = dft[0][0]*dft[0][0] + dft[0][1]*dft[0][1];
    dft[1][1] = dft[1][0]*dft[1][0] + dft[1][1]*dft[1][1];
    dft[2][1] = dft[2][0]*dft[2][0] + dft[2][1]*dft[2][1];
    dft[3][1] = dft[3][0]*dft[3][0] + dft[3][1]*dft[3][1];
    unsigned const DFT = dft[2][1] > dft[3][1] ? 2 : 3;
    float const F0 = ( DFT == 2 ? 50.0 : 60.0);
    dft[DFT][1] /= 50.0f*magnitude - 0.5f*lev*lev;
    if (dft[DFT][1] > 1.0f) {
        dft[DFT][1] = 1.0f;
    }
    if ((dft[DFT - 2][0] - dft[DFT - 1][0]) * (dft[DFT - 1][0] - dft[DFT][0]) < 0) {
        return F0 - DF*std::sqrt(1.0 - dft[DFT][1]);
    }
    return F0 + DF*std::sqrt(1.0 - dft[DFT][1]);
} //calculateFrequency()

int main() {
    std::srand(std::time(0));
    unsigned const NN = 1000;
    double samples[NN], re[NN], im[NN];
    struct {
        unsigned N, k, m;
        double dev, tobe;
    } F = {0}, D = F, current = F;
    try {
        for (unsigned N = 1; N < NN/7; ++N) {
            std::cerr << N;
            if (!(N % 30)) {
                std::cerr << "/" << NN << "\n";
            }
            std::cerr << " ";
            DFT dft(N);
            //dft.printInfo(std::cerr);
            current.N = N;
            current.dev = dft.dev(&current.k, &current.m, &current.tobe);
            if (std::fabs(current.dev) > std::fabs(F.dev)) {
                F = current;
            }
            generate(N, rnd(0.1, 0.5), rnd(10.0, 100.0)
                    , rnd(0.0, 8.0*std::atan(1.0)), rnd(0.0005, 0.05)
                    , rnd(-1.0, 1.0), samples);
            current.dev = dft.apply(samples, re, im, &current.k);
            current.m = N;
            if (std::fabs(current.dev) > std::fabs(D.dev)) {
                D = current;
            }
        } //for
        std::cout << "\n";
        messageDev(std::cout, F.N, F.m, F.k, F.dev, "DFT matrix entries");
        messageDev(std::cout, D.N, D.m, D.k, D.dev, "DFT result");
        DFT dft100(100);
        dft100.code(std::cout, "float");
        std::cout << "checkDFT100(): " << checkDFT100() << "\n";
        std::cout << "frequency\tdeviation\t40\t50\t60\t70\n";
        std::cout << std::setprecision(6);
        float max_dev = 0.0f, fcalc, f_max_dev;
        for (double frequency = 57.5; frequency < 63.0; frequency += 0.01) {
            float magnitude = rnd(0.5, 1.5);
            float level = rnd(-1.0, 1.0);
            float phase = rnd(0.0, 8.0*std::atan(1.0));
            float samples[100];
            generate<float>(100, magnitude, frequency
                    , phase, 0.001
                    , level, samples);
            fcalc = 10.0f;
            fcalc = calculateFrequency(samples) - frequency;
            float re30, re40, re50, re60, re70, re80, im;
            calculateDFT100(3, samples, &re30, &im);
            re30 = std::sqrt(re30*re30 + im*im)/100.0f/magnitude;
            calculateDFT100(4, samples, &re40, &im);
            re40 = std::sqrt(re40*re40 + im*im)/100.0f/magnitude;
            calculateDFT100(5, samples, &re50, &im);
            re50 = std::sqrt(re50*re50 + im*im)/100.0f/magnitude;
            calculateDFT100(6, samples, &re60, &im);
            re60 = std::sqrt(re60*re60 + im*im)/100.0f/magnitude;
            calculateDFT100(7, samples, &re70, &im);
            re70 = std::sqrt(re70*re70 + im*im)/100.0f/magnitude;
            calculateDFT100(8, samples, &re80, &im);
            re80 = std::sqrt(re80*re80 + im*im)/100.0f/magnitude;
            std::cout << frequency << "\t" << fcalc << "\t" << re40 << "\t" << re50 << "\t" << re60 << "\t" << re70 << "\t" << re80 <<"\n";
            if (std::fabs(max_dev) < std::fabs(fcalc)) {
                max_dev = fcalc;
                f_max_dev = frequency;
            }
        }
        std::cout << "Maximum deviation is " << max_dev << "Hz for frequency "
                << f_max_dev << "Hz\n";
    } catch (std::exception& e) {
        std::cerr << "Caught std::exception: " << e.what() << "\n";
    } catch (...) {
        std::cerr << "Caught unknown exception\n";
    }
    return 0;
}
