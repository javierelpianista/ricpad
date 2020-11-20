#include <ricpad.hpp>

int Dmin, Dmax, d, digits;
gi::symbol x;
gi::ex potential;
map<string, string> options;

int main() {
    read_input("input");

    options["problem_type"] = "even";
    options["tolerance"] = "1E-100";
    options["diff_epsilon"] = "1E-150";

    digits = 200;
    gi::Digits = digits;
    mpfr_float::default_precision(digits);

    cout << "Potential: " << potential << endl;
    cout << "Variable: " << x << endl;
    cout << potential.series(x == 0, 10) << endl;

    vector<mpfr_float> Q = get_coefficients(potential, x, 2*Dmax);
    cout.precision(digits);

    for ( int i = Dmin; i<= Dmax; i++ ) {
        E0 = RPM_solve(i, d, Q, 0, E0);
        cout << "D = " << i << " E0 = " << E0 << endl;
    }
    
    return 0;
}
