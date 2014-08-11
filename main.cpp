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
        //std::cout << "checkDFT100(): " << checkDFT100() << "\n";
    } catch (std::exception& e) {
        std::cerr << "Caught std::exception: " << e.what() << "\n";
    } catch (...) {
        std::cerr << "Caught unknown exception\n";
    }
    return 0;
}
