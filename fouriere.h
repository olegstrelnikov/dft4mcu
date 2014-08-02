/**
 * DFC header
 *
 */

#ifndef FOURIERE_H_
#define FOURIERE_H_

#include <iostream>
#include <vector>
#include <map>
#include <string>

/**
 * DFT matrix class
 */
class DFT {
public:
    enum PART {RE, IM};
private:
    unsigned const N;
    double const rN;
    unsigned const SHIFT;
    unsigned const SZ;
    std::vector<double> cos;
    unsigned k, m;
    PART part;
    static double two_pi_q_divide_r(unsigned q, double r);
    void dev();
    double applyPrecise(double const*, unsigned, PART) const;
    void apply(unsigned, double const*, double*, double*) const;
    double devResult(double const*, double const*, double const*
            , unsigned*) const;
public:
    DFT(unsigned n) throw (std::invalid_argument);
    double apply(double const*, double*, double*, unsigned*) const;
    double applyReverse(double* samples, double const* re, double const* im) const;
    void printResult(std::ostream& os, double const* re, double const* im) const;
    void printInfo(std::ostream& os) const;
    double re(unsigned k, unsigned m) const;
    double im(unsigned k, unsigned m) const;
    double preciseRe(unsigned k, unsigned m) const;
    double preciseIm(unsigned k, unsigned m) const;
    double dev(unsigned* pk, unsigned* pm, double* tobe) const;
    double dev(unsigned k, unsigned m, double* tobe, PART) const;
    void code(std::ostream& os, std::string const& type) const;
}; //class DFT

#endif /* FOURIERE_H_ */
