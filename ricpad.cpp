#include <iostream>

#include <problem.hpp>
#include <options.hpp>
#include <read_input.hpp>
#include <ricpad.hpp>

using boost::multiprecision::mpfr_float;
using boost::multiprecision::mpc_complex;

using namespace std;

int main(int argc, char * argv[]) {
    Options opts;
    Problem problem;

    if ( argc == 1 ) {
        cout << "Usage: ricpad [input] [options]" << endl;
        return 1;
    }

    read_input(argv[1], opts, problem);

    if ( opts.ints.at("use_complex") ) {
        RPM_solve<mpc_complex>(problem, opts);
    } else {
        RPM_solve<mpfr_float>(problem, opts);
    }
}
