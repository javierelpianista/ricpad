#ifndef PROBLEM
#define PROBLEM

#include <map>
#include <string>

#include <ginac/ginac.h>
#include <boost/multiprecision/mpfr.hpp>

#include <options.hpp>

using namespace std;
namespace gi = GiNaC;

template <class T>
T gi_to_mp(
        const gi::numeric num
        ) 
{
    ostringstream oss;
    oss << num;

    return T(oss.str());
}


class Problem {
    gi::symbol x;
    gi::ex potential;

    public: 
        Problem() = default;
        ~Problem() = default;

        void set_potential(
                const string potential_str,
                const string variable_str
                ) {

            x = gi::symbol(variable_str);
            gi::symtab table;
            table[variable_str] = x;
            gi::parser reader(table);

            potential = reader(potential_str);
        };

        void print() {
            cout << potential << " (" << x << ")" << endl;
        };

        template <class T>
        vector<T> get_coefficients(
                int N, 
                const Options & options
                ) const
        {
            string problem_type;
            int use_rationals = 0;
            vector<T> data;
            data.reserve(N);

            problem_type = options.strings.at("problem_type");
            use_rationals = options.ints.at("use_rationals");

            gi::numeric to_push;
            gi::ex Q_series = gi::series_to_poly(
                    potential.series(x == 0, N+3)
                    );
            gi::ex tmp;

            int mult;

            if ( problem_type == "even" ) mult = 2;
            else if ( problem_type == "radial" ) mult = 1;

            for ( int i = 0; i <= N; i++ ) {
                if ( use_rationals == 1 ) {
                    to_push = gi::ex_to<gi::numeric>(Q_series.coeff(x, mult*i));
                } else {
                    tmp = Q_series.coeff(x, mult*i).evalf();
                    to_push = gi::ex_to<gi::numeric>(Q_series.coeff(x, mult*i).evalf());
                }
                data.push_back(gi_to_mp<T>(to_push));
            }

            return data;
        };

        template <class T>
        T get_neg_coeff(
                const Options & opts
                ) const
        {
            gi::ex coeff = potential.series(x, 1).coeff(x, -1);
            return gi_to_mp<T>(gi::ex_to<gi::numeric>(coeff));
        }
};
#endif
