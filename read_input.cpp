#include <map>
#include <string>
#include <fstream>

#include <boost/algorithm/string.hpp>
#include <boost/multiprecision/mpfr.hpp>

#include <options.hpp>
#include <problem.hpp>

using namespace std;
namespace mp = boost::multiprecision;

void read_input(const string filename, Options &opts, Problem &problem) {
    map<string, string> parsed_options;
    map<string, string>::iterator it;

    ifstream read_file(filename);
    string line;
    vector<string> tokens;

    while ( getline( read_file, line ) ) {
        boost::split(tokens, line, boost::is_any_of(" "));

        auto p = tokens.begin();
        string what = *p++;
        string value;

        if ( what == "D" ) {
            opts.ints["Dmin"] = stoi(*p++);
            //TODO program the option to not have Dmax
            opts.ints["Dmax"] = stoi(*p);
        } else if ( what == "d" ) {
            opts.ints["d"], stoi(*p);
        } else if ( what == "s" || what == "l" ) {
            opts.ints["s"] = stoi(*p);
        } else if ( what == "digits" ) {
            opts.ints["digits"] = stoi(*p);
        } else if ( what == "target_digits" ) {
            if ( *p == "max" ) {
                opts.ints["infinite_digits"] = 1;
            } else {
                opts.ints["target_digits"] = stoi(*p);
            }
        } else if ( what == "nr_max_iter" ) {
            opts.ints["nr_max_iter"] = stoi(*p);
        } else if ( what == "nr_tolerance" ) {
            opts.mpfrs["nr_tolerance"] = mpfr_float(*p);
        } else if ( what == "nr_diff_epsilon" ) {
            opts.mpfrs["nr_step_size"] = mpfr_float(*p);
        } else {
            while ( p != tokens.end() ) {
                value += *p++;
            }
            parsed_options[what] = value;
        }
    }

    // We will always use nr_tolerance = 10^(-digits + 180) and
    // nr_step_size = 10^(-digits + 100)

    // If target_digits is set, and the other relevant accuracy variables
    // are not, set NR tolerance to 10^(-(target_digits + 20)),  
    // nr_step_size to 10^(-target_digits + 100) and digits to 
    // target_digits + 200
    if ( opts.ints["target_digits"] != -1 ) {
        cout << "target_digits is: " << opts.ints["target_digits"] << endl;
        if ( opts.ints["digits"] == -1 ) 
            opts.ints["digits"] = opts.ints["target_digits"] + 200;

        mpfr_float::default_precision(opts.ints["digits"]);
        if ( opts.mpfrs["nr_tolerance"] == -1 ) {
            ostringstream oss;
            oss << "1E-" << opts.ints["target_digits"] + 20;
            opts.mpfrs["nr_tolerance"] = mpfr_float(oss.str());
        }
        if ( opts.mpfrs["nr_step_size"] == -1 ) {
            ostringstream oss;
            oss << "1E-" << opts.ints["target_digits"] + 100;
            opts.mpfrs["nr_step_size"] = mpfr_float(oss.str());
        }
    } 
    // If it hasn't been set, check if digits has been defined. If so, 
    // set nr_tolerance to 10^(-digits+180) and nr_step_size to 
    // 10^(-digits+100). 
    else if ( opts.ints["digits"] != -1 ) {
        if ( opts.ints["digits"] < 200 ) {
            cout << "digits should be at least 200" << endl;
            throw 1;
        }

        mpfr_float::default_precision(opts.ints["digits"]);

        opts.mpfrs["nr_tolerance"] = mpfr_float(
                "1E-" + to_string(opts.ints["digits"] - 180)
                );

        opts.mpfrs["nr_step_size"] = mpfr_float(
                "1E-" + to_string(opts.ints["digits"] - 100)
                );
    } 
    // Else, set target_digits to 40 and everything else accordingly
    else {
        opts.ints["digits"] = 240;
        mpfr_float::default_precision(opts.ints["digits"]);

        if ( opts.ints["infinite_digits"] == 0 ) 
            opts.ints["target_digits"] = 40;


        if ( opts.mpfrs["nr_tolerance"] == -1 ) 
            opts.mpfrs["nr_tolerance"] = mpfr_float("1E-60");

        if ( opts.mpfrs["nr_step_size"] == -1 )
            opts.mpfrs["nr_step_size"] = mpfr_float("1E-140");
    }


    // **********************************************************************
    // Now take care of options that define the problem or require 
    // initialization of previous variables
    // **********************************************************************
    
    it = parsed_options.find("E0");
    if ( it == parsed_options.end() ) {
        cout << "The value E0 should be given" << endl;
        throw 10;
    }
    
    it = parsed_options.find("var");
    if ( it == parsed_options.end() ) {
        cout << "The variable should be defined with var <VAR>" << endl;
        throw 10;
    } 
    
    it = parsed_options.find("pot");
    if ( it == parsed_options.end() ) {
        cout << "Potential should be defined with pot <POT>" << endl;
        throw 10;
    }

    // Set the precision of the output
    if ( opts.ints["target_digits"] > 0 ) {
        cout.precision(opts.ints["target_digits"]);
        }
    else {
        mpfr_float ndigits;
        ndigits = -log10(opts.mpfrs["nr_tolerance"]);
        cout.precision(static_cast<int>(ndigits));
    }
    opts.mpfrs["E0"] = mpfr_float(parsed_options["E0"]);
    problem.set_potential(parsed_options["pot"], parsed_options["var"]);
}
