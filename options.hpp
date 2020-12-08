#ifndef OPTIONS
#define OPTIONS

#include <map>
#include <string>
#include <iostream>
#include <boost/multiprecision/mpfr.hpp>

using namespace std;
using boost::multiprecision::mpfr_float;

class Options {
    public:

    map<string, int> ints;
    map<string, double> doubles;
    map<string, string> strings;
    map<string, mpfr_float> mpfrs;
    
    // For now we set the defaults on construction
    Options() {
        ints["use_rationals"] = 0;
        ints["use_complex"] = 1;

        ints["target_digits"] = -1;
        ints["Dmin"] = 3;
        ints["Dmax"] = -1;
        ints["d"] = 0;
        ints["digits"] = -1;
        ints["infinite_digits"] = 0;
        ints["log_nr"] = 0;

        ints["nr_max_iter"] = 100;

        strings["problem_type"] = "even";

        mpfrs["nr_step_size"] = mpfr_float("-1");
        mpfrs["nr_tolerance"] = mpfr_float("-1");
        mpfrs["E0I"] = mpfr_float(0);
    }

    ~Options() = default;

    void print() const {
        for ( auto const &x : ints ) 
            cout << x.first << " " << x.second << endl;
        for ( auto const &x : doubles ) 
            cout << x.first << " " << x.second << endl;
        for ( auto const &x : strings ) 
            cout << x.first << " " << x.second << endl;
        for ( auto const &x : mpfrs ) 
            cout << x.first << " " << x.second << endl;
    }
};
#endif
