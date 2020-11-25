#ifndef RICPAD
#define RICPAD

#include <string>
#include <iostream>
#include <ginac/ginac.h>
#include <boost/multiprecision/mpfr.hpp>
//#include <boost/multiprecision/mpc.hpp>

#include <problem.hpp>

namespace mp = boost::multiprecision;
namespace gi = GiNaC;

// Convert a number from a GiNaC::numeric to T

// Get the coefficients of the Taylor expansion of Q in a vector
template <class T>
vector<T> get_coefficients(
        const gi::ex Q, 
        const gi::ex &x, 
        const int N, 
        const bool use_rationals = false, 
        const string problem_type = "even"
        ) 
{
    vector<T> data;

    gi::numeric to_push;
    gi::ex Q_series = gi::series_to_poly(Q.series(x == 0, N+3));
    gi::ex tmp;
    T to_push_n;

    int mult;

    if ( problem_type == "even" ) mult = 2;

    for ( int i = 0; i <= N; i++ ) {
        if ( use_rationals ) {
            to_push = gi::ex_to<gi::numeric>(Q_series.coeff(x, mult*i));
        } else {
            tmp = Q_series.coeff(x, mult*i).evalf();
            to_push = gi::ex_to<gi::numeric>(Q_series.coeff(x, mult*i).evalf());
        }
        to_push_n = gi_to_mp<T>(to_push);
        data.push_back(to_push_n);
    }

    return data;
}


// Return a vector with the Hankel coefficients up to N
template <class T>
vector<T> hankcoefs(
        int N, 
        const Problem &problem, 
        int s,
        const T E0, 
        const Options & opts
        ) 
{
    vector<T> coefs;
    coefs.reserve(N+1);
    vector<T> v_coefs = problem.get_coefficients<T>(N, opts);
    coefs.push_back(E0 - v_coefs[0]);

    {
        T sum;
        for ( int j = 1; j <= N; j++ ) {
            sum = 0;
            for ( int k = 0; k <= j-1; k++ ) {
                sum += coefs[k]*coefs[j-k-1];
            }
            coefs.push_back((sum - v_coefs[j])/(2*j+2*s+1));
        }
    }

    return coefs;
}

// Return the hankel determinant with D and d, for Problem problem and parameter s,
// evaluated at point E0
template <class T>
T hankdet(
        const int D, 
        const int d, 
        const Problem &problem, 
        const int s,
        const T E0, 
        const Options & opts 
        ) {
    vector<T> coefs, coefsm1, coefsm2;

    if ( D == 0 ) {
        return 1;
    } else {
        coefsm1 = std::vector<T>(2*D, T(1));

        //for ( int i = d+1; i <= 2*D + d - 1; i++ ) 
        //    coefs.emplace_back(hankcoef(i, V, s, E0));

        coefs = hankcoefs<T>(2*D + d - 1, problem, s, E0, opts);
        coefs.erase(coefs.begin(), coefs.begin()+d+1);

        if ( D == 1 ) return coefs[0];

        for ( int j = 2; j <= D; j++ ) {
            coefsm2 = (vector<T>&&)(coefsm1);
            coefsm1 = (vector<T>&&)(coefs);

            for ( int k = 0; k <= 2*(D-j); k++ ) {
                coefs.emplace_back(
                        (coefsm1[k]*coefsm1[k+2] - 
                        coefsm1[k+1]*coefsm1[k+1]) /
                        coefsm2[k+2]
                        );
            }
        }
    }

    return coefs[0];
}


template <class T>
T RPM_find_root(
        const int D, 
        const int d, 
        const Problem &problem,
        const int s,
        const T E0,
        const Options & opts
        )
{
    T E = E0, Eold;
    T tol(opts.mpfrs.at("nr_tolerance"));
    T desv = tol + 1;
    T h(opts.mpfrs.at("nr_step_size"));
    T dH, H;

    int niter = 0;

    while ( desv > tol ) {
        Eold = E;
        H = hankdet<T>(D, d, problem, s, E, opts);
        dH = (hankdet<T>(D, d, problem, s, E + h, opts) - hankdet<T>(D, d, problem, s, E - h, opts));
        dH /= 2*h;

        E = Eold - H/dH;
        desv = mp::abs(Eold - E);

        if ( niter++ >= opts.ints.at("nr_max_iter") ) {
            cout << "Maximum number of NR iterations reached. Aborting." << endl;
            throw 2;
        }
    }

    return E;
}

template <class T>
T RPM_solve(const Problem & problem, const Options & opts) {
    T E = T(opts.mpfrs.at("E0"));
    T Eold = 0;
    T logdif = 0; 
    int d = opts.ints.at("d");
    int s = opts.ints.at("s");

    int accurate_digits = 0;

    for ( int D = opts.ints.at("Dmin"); D <= opts.ints.at("Dmax"); D++ ) {
        Eold = E;
        E = RPM_find_root<mpfr_float>(D, d, problem, s, E, opts);
        logdif = -log10(abs(E - Eold));
        cout << D << " " << E << endl;

        if ( 
                opts.ints.at("infinite_digits") == 0 &&
                static_cast<int>(logdif) >= opts.ints.at("target_digits") 
            ) 
            break;
    }

    return E;
}

#endif
