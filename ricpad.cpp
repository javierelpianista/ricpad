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
        return 2;
    }

    read_input(argv[1], opts, problem);

    opts.print();

    int result;
    if ( opts.ints.at("use_complex") ) {
        mpc_complex E;
        result = RPM_solve<mpc_complex>(problem, opts, E);
    } else {
        mpfr_float E;
        result = RPM_solve<mpfr_float>(problem, opts, E);
    }

    return result;
}
