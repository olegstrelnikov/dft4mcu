/**
 * Discrete Fouriere Transform (DFT) optimal calculator generator
 */

#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <iomanip>
#include <sstream>
#include "fouriere.h"

inline double DFT::two_pi_q_divide_r(unsigned q, double r) {
    return (q << 3) * std::atan(1.0) / r;
} //DFT::two_pi_q_divide_r()

/**
 * Deviation of the one DFT matrix entry
 * @param k entry's row
 * @param m entry's column
 * @param tobe precise value of the entry, deviation is calculated wrt it
 * @param part part of DFT matrix : Re or Im
 * @return deviation
 */
double DFT::dev(unsigned const k, unsigned const m, double* const tobe
        , PART const part) const {
    *tobe = (part == RE ? this->preciseRe(k, m) : this->preciseIm(k, m));
    return std::fabs( (part == RE ? this->re(k, m) : this->im(k, m)) - *tobe );
} //DFT::dev()

/**
 * Maximum deviation of the one DFT matrix entry
 * @param out pk where maximum deviation entry's row will be stored to
 * @param out pm where maximum deviation entry's column will be stored to
 * @param tobe precise value of the entry, deviation is calculated wrt it
 * @return deviation
 */
double DFT::dev(unsigned* const pk, unsigned* const pm
        , double* const tobe) const {
    *pk = this->k;
    *pm = this->m;
    return this->dev(*pk, *pm, tobe, this->part);
} //DFT::dev()

void DFT::dev() {
    double delta = 0.0;
    this->k = 0;
    this->m = 0;
    this->part = DFT::RE;
    for (unsigned k = 0; k < N; ++k) {
        for (unsigned m = 0; m < N; ++m) {
            if (std::fabs(this->re(k, m) - this->preciseRe(k, m)) > delta) {
                delta = std::fabs(this->re(k, m) - this->preciseRe(k, m));
                this->k = k;
                this->m = m;
                this->part = DFT::RE;
            }
            if (std::fabs(this->im(k, m) - this->preciseIm(k, m)) > delta) {
                delta = std::fabs(this->im(k, m) - this->preciseIm(k, m));
                this->k = k;
                this->m = m;
                this->part = DFT::IM;
            }
        }
    }
} //DFT::dev()

double DFT::re(unsigned const k, unsigned const m) const {
    unsigned icos = ((4*k*m) % (2*N)) >> SHIFT;
    if (icos > SZ) {
        icos = 2*SZ - icos;
    }
    unsigned Q = 4*(k*m % N);
    if (Q < N || Q > 3*N) {
        return this->cos[icos];
    } else if (Q > N && Q < 3*N) {
        return -this->cos[icos];
    } else if ( k == 0 || m == 0 ) {
        return 1.0;
    } else if ( k > N/2 ) {
        return this->re(N - k, m);
    } else if (2*k == N) {
        return ( (m & 1) ? -1.0 : 1.0 );
    }
    return 0.0;
} //DFT::re()

double DFT::im(unsigned const k, unsigned const m) const {
    unsigned icos = ((4*k*m) % (2*N)) >> SHIFT;
    if (icos > SZ) {
        icos = 2*SZ - icos;
    }
    unsigned isin = SZ - icos;
    unsigned Q = 4*(k*m % N);
    if (Q < 2*N) {
        return this->cos[isin];
    } else if (Q > 2*N) {
        return -this->cos[isin];
    } else if ( k > N/2 ) {
        return -this->im(N - k, m);
    }
    return 0.0;
} //DFT::im()

inline double DFT::preciseRe(unsigned const k, unsigned const m) const {
    return std::cos(2.0*std::atan(1.0)*4.0*k*m/rN);
} //DFT::preciseRe()

inline double DFT::preciseIm(unsigned const k, unsigned const m) const {
    return std::sin(2.0*std::atan(1.0)*4.0*k*m/rN);
} //DFT::preciseIm()

double DFT::applyPrecise(double const* const samples,
        unsigned const k, PART part) const {
    double y = 0.0;
    if (part == DFT::RE) {
        for (unsigned m = 0; m < N; ++m) {
            y += samples[m] * std::cos(DFT::two_pi_q_divide_r(k*m, rN));
        }
    } else {
        for (unsigned m = 0; m < N; ++m) {
            y += samples[m] * std::sin(DFT::two_pi_q_divide_r(k*m, rN));
        }
    }
    return y;
} //DFT::applyPrecise()

double DFT::devResult(double const* const samples, double const* const re
        , double const* const im, unsigned* const pk) const {
    double dev = 0.0, delta;
    for (unsigned k = 0; k < N; ++k) {
        delta = re[k] - applyPrecise(samples, k, DFT::RE);
        if (std::fabs(delta) > std::fabs(dev)) {
            dev = delta;
            if (pk) {
                *pk = k;
            }
        }
        delta =im[k] - applyPrecise(samples, k, DFT::IM);
        if (std::fabs(delta) > std::fabs(dev)) {
            dev = delta;
            if (pk) {
                *pk = k;
            }
        }
    } //for k
    return dev;
} //DFT::devResult()

static double rnd(double const from, double const to) {
    return from + std::rand() / static_cast<double>(RAND_MAX) * (to - from);
} //rnd()

void DFT::code(std::ostream& os, std::string const& type) const {
    unsigned static const COSINE_COLS = 3;
    os
        << "/**\n"
        << " * Calculate Discrete Fouriere Transform (DFT) of size " << this->N << "\n"
        << " * @param k no of DFT result's entry (has to be not more then " << (N/2) << ")\n"
        << " * @param samples array (0.." << (this->N - 1) << ") of DFT's source\n"
        << " * @param re where the real part of the result to store to\n"
        << " * @param im where the imaginary part of the result to store to\n"
        << " * To calculate DFT for k > " << (N/2) << ", please,"
        << " call calculateDFT(" << this->N << " - k, re, im)"
        << " and then negate im\n"
        << " */\n"
        ;
    os
        << "void calculateDFT" << this->N << "("
        << "unsigned k, "
        << type << " const* const samples, "
        << type << "* const re, " << type << "* const im"
        << ");\n\n"
        ;
    os
        << "void calculateDFT" << this->N << "("
        << "unsigned k, "
        << type << " const* const samples, "
        << type << "* const re, " << type << "* const im"
        << ") {\n"
        ;
    os
        << "    static const unsigned N = " << this->N << ";\n";
    os
        << "    static const unsigned SZ = " << this->SZ << ";\n";
    os
        << "    static const unsigned SHIFT = " << this->SHIFT << ";\n";
    os
        << "    static const " << type << " COSINE[" << this->cos.size() << "] = {\n";
    for (std::vector<double>::size_type i = 0; i <= this->cos.size() / COSINE_COLS; ++i) {
        for (size_t j = 0; j < COSINE_COLS && 3*i + j < this->cos.size(); ++j) {
            os
                << (j ? ", " : "") << std::setw(25) << std::setprecision(19)
                << this->cos[3*i + j] << (type == "float" ? "F" : "")
                ;
        }
        os
            << (i < this->cos.size() / COSINE_COLS ? "," : "") << "\n";
    }
    os
        << "    }; /* COSINE */\n";
    os
        << "    unsigned m, icos, isin, Q, md2;\n";
    os
		<< "    *re = samples[0];\n"
		<< "    *im = 0;\n"
		<< "    if (k == 0) {\n"
        << "        for (m = 1; m < N; ++m) {\n"
        << "            *re += samples[m];\n"
        << "        }\n"
		<< "    } else if (k < (N + 1)/2) {\n"
        << "        for (m = 1; m < N; ++m) {\n"
        << "            icos = ((4*k*m) % (2*N)) >> SHIFT;\n"
        << "            if (icos > SZ) {\n"
        << "                icos = 2*SZ - icos;\n"
        << "            }\n"
        << "            isin = SZ - icos;\n"
        << "            Q = 4*(k*m % N);\n"
        << "            if (Q < N || Q > 3*N) {\n"
        << "                *re += COSINE[icos] * samples[m];\n"
        << "            } else if (Q > N && Q < 3*N) {\n"
        << "                *re -= COSINE[icos] * samples[m];\n"
        << "            }\n"
        << "            if (Q < 2*N) {\n"
        << "                *im += COSINE[isin] * samples[m];\n"
        << "            } else if (Q > 2*N) {\n"
        << "                *im -= COSINE[isin] * samples[m];\n"
        << "            }\n"
        << "        }\n"
        << "    } else if (k == N/2) {\n"
        << "        for (md2 = 0; md2 < N/2; ++md2) {\n"
        << "            *re += samples[2*md2 + 2] - samples[2*md2 + 1];\n"
        << "        }\n"
        << "    }\n"
        ;
    os
        << "} /* calculateDFT" << this->N << "() */\n";
    os
        << "\n"
        << "/**\n"
        << " * Check calculateDFT" << this->N << "() function\n"
        << " * @return maximum deviation from precise value\n"
        << " */\n"
        ;
    os
        << type << " checkDFT" << this->N << "("
        << ");\n\n"
        ;
    os
        << type << " checkDFT" << this->N << "("
        << ") {\n"
        ;
    os
        << "    static const unsigned N = " << this->N << ";\n";
    os
        << "    static const " << type << " rN = N;\n";
    os
        << "    static const " << type << " DPHI = " << rnd(0.005, 5.0) << ";\n";
    os
        << "    static const " << type << " LVL = " << rnd(-1.0, 1.0) << ";\n";
    os
        << "    static const " << type << " MAG = " << rnd(0.5, 1.5) << ";\n";
    os
        << "    " << type << " samples[N], re, im, phi = "
        << rnd(-2.0, 2.0) << ", result = 0.0;\n";
    os
        << "    unsigned k, m;\n";
    os
        << "    for (k = 0; k < N; ++k, phi += DPHI) {\n"
        << "        samples[k] = LVL + MAG*std::sin(phi);\n"
        << "    }\n"
        << "    for (k = 0; k < (N + 1)/2; ++k) {\n"
        << "        calculateDFT" << N << "(k, samples, &re, &im);\n"
        << "        for (m = 0; m < N; ++m) {\n"
        << "            re -= samples[m] * std::cos(2.0*std::atan(1.0)*4.0*k*m/rN);\n"
        << "            im -= samples[m] * std::sin(2.0*std::atan(1.0)*4.0*k*m/rN);\n"
        << "        }\n"
        << "        if (re > result) {\n"
        << "            result = re;\n"
        << "        }\n"
        << "        if (im > result) {\n"
        << "            result = im;\n"
        << "        }\n"
        << "    }\n"
        << "    return result;\n"
        << "} /* checkDFT" << N << "() */\n"
        ;
}

void DFT::apply(unsigned k, double const* const samples, double* const re
        , double* const im) const {
    *re = samples[0];
    *im = 0;
    if (k == 0) {
        for (unsigned m = 1; m < N; ++m) {
            *re += samples[m];
        }
    } else if (k < (N + 1)/2) {
        for (unsigned m = 1; m < N; ++m) {
            unsigned icos = ((4*k*m) % (2*N)) >> SHIFT;
            if (icos > SZ) {
                icos = 2*SZ - icos;
            }
            unsigned isin = SZ - icos;
            unsigned Q = 4*(k*m % N);
            if (Q < N || Q > 3*N) {
                *re += this->cos[icos] * samples[m];
            } else if (Q > N && Q < 3*N) {
                *re -= this->cos[icos] * samples[m];
            }
            if (Q < 2*N) {
                *im += this->cos[isin] * samples[m];
            } else if (Q > 2*N) {
                *im -= this->cos[isin] * samples[m];
            }
        }
    } else if (k == N/2) {
        for (unsigned md2 = 0; md2 < N/2; ++md2) {
            *re += samples[2*md2 + 2] - samples[2*md2 + 1];
        }
    }
} //DFT::apply()

double DFT::apply(double const* const samples, double* const re
        , double* const im, unsigned* const pk) const {
    for (unsigned k = 0; k < N/2 + 1; ++k) {
        this->apply(k, samples, re + k, im + k);
    }
    for (unsigned k = N/2 + 1; k < N; ++k) {
        re[k] = re[N - k];
        im[k] = -im[N - k];
    }
    return devResult(samples, re, im, pk);
} //DFT::apply()

void DFT::printInfo(std::ostream& os) const {
    os << "DFT(N = " << N << "; cos.size() = " << cos.size()
        << ", max deviation "
        << (part == RE ? this->re(this->k, this->m) - this->preciseRe(this->k, this->m) : this->im(this->k, this->m) - this->preciseIm(this->k, this->m))
        << " is in entry ( " << this->k << ", " << this->m
        << " ) )\n";
} //DFT::printInfo()

void DFT::printResult(std::ostream& os, double const* re, double const* im) const {
    unsigned k = 0;
    for (double const* end = re + N; re < end; ++re, ++im, ++k) {
        os << std::setw(10) << k << " [" << std::setw(15) << *re << "] ["
                << std::setw(15) << *im << "]\n";
        ;
    }
} //DFT::printResult()

DFT::DFT(unsigned n) throw (std::invalid_argument)
        : N(n), rN(static_cast<double>(N))
        , SHIFT((N & 3 ? (N & 1 ? 0 : 1) : 2)), SZ(N >> SHIFT), cos(1 + SZ)
        , k(0), m(0), part(RE) {
    if (!N) {
        throw std::invalid_argument("Problem dimension must be greater then 0");
    }
    for (unsigned j = 0; j <= SZ; ++j) {
    	this->cos[j] = std::cos(DFT::two_pi_q_divide_r(j << SHIFT, 4.0*rN));
    }
    this->dev();
} //DFT::DFT()
