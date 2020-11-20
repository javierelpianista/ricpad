#ifndef RICPAD
#define RICPAD

#include <iostream>
#include <ginac/ginac.h>
#include <boost/multiprecision/mpfr.hpp>

namespace mp = boost::multiprecision;
namespace gi = GiNaC;

using namespace std;
using boost::multiprecision::mpfr_float;

extern map<string, string> options;
extern int Dmin, Dmax, d, digits;
extern gi::symbol x;
extern gi::ex potential;
extern mpfr_float E0;

mpfr_float gi_to_mpfr(
        const gi::numeric
        );

vector<mpfr_float> get_coefficients(
        const gi::ex, 
        const gi::ex &, 
        const int,
        const bool = false
        );
vector<mpfr_float> hankcoefs(
        const int,
        const vector<mpfr_float> &,
        const int,
        const mpfr_float
        );
mpfr_float hankdet(
        const int,
        const int,
        const vector<mpfr_float> &V,
        const int,
        const mpfr_float
        );
mpfr_float RPM_solve(
        const int,
        const int,
        const vector<mpfr_float> &,
        const int,
        const mpfr_float
        );
void read_input(
        string
        );
#endif
