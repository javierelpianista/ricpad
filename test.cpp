#include <ricpad.hpp>
#include <options.hpp>

namespace mp = boost::multiprecision;
namespace gi = GiNaC;

using namespace std;
using boost::multiprecision::mpfr_float;

int Dmin, Dmax, d, digits, s;
gi::symbol x;
gi::ex potential;
map<string, string> options;

int main() {
    Options opts;
    read_input("input", opts);

    options["problem_type"] = "even";
    options["tolerance"] = "1E-100";
    options["diff_epsilon"] = "1E-150";

    digits = 200;
    gi::Digits = digits;
    mpfr_float::default_precision(digits);

    cout << "Potential: " << potential << endl;
    cout << "Variable: " << x << endl;
    cout << "s: " << s << endl;

    vector<mpfr_float> Q = get_coefficients<mpfr_float>(potential, x, 2*Dmax);
    cout.precision(digits);

    for ( int i = Dmin; i<= Dmax; i++ ) {
        E0 = RPM_solve(i, d, Q, s, E0);
        cout << "D = " << i << " E0 = " << E0 << endl;
    }
    
    return 0;
}
