#include <options.hpp>
#include <problem.hpp>
#include <ricpad.hpp>

void read_input(const string, Options &, Problem &);

int main() {
    vector<mpfr_float> v;
    Options opts;
    Problem problem;

    read_input("input", opts, problem);

    mpfr_float E0 = opts.mpfrs["E0"];

    cout << "Options:" << endl;
    opts.print();
    cout << endl;

    RPM_solve<mpfr_float>(problem, opts);
}
