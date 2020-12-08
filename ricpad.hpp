#ifndef RICPAD
#define RICPAD

#include <string>
#include <iostream>
#include <fstream>
#include <ginac/ginac.h>
#include <boost/multiprecision/mpfr.hpp>
#include <boost/multiprecision/mpc.hpp>

#include <problem.hpp>

namespace mp = boost::multiprecision;
namespace gi = GiNaC;

using mp::mpfr_float;
using mp::mpc_complex;

// Assign the initial value to E0
template <class T>
void assign_E0( T & E0, const Options & opts ) {
}

template <>
void assign_E0<mpfr_float>( mpfr_float & E0, const Options & opts ) {
    E0 = opts.mpfrs.at("E0");
}

template <>
void assign_E0<mpc_complex>( mpc_complex & E0, const Options & opts ) {
    E0.real(opts.mpfrs.at("E0"));
    E0.imag(opts.mpfrs.at("E0I"));
}

// Return an adequate step size for each numeric type. If it's a complex type,
// then return (1+I)/sqrt(2)*h, if it is a real type, return h.
template <class T> 
T NR_step_size(const Options & opts) 
{
    T ans;
    cout << "NR_step_size not defined for the template argument chosen." 
        << endl;
    throw 2;
    return ans;
}

template <> 
mpfr_float NR_step_size<mpfr_float>(const Options & opts) {
    return opts.mpfrs.at("nr_step_size");
}

template <> 
mpc_complex NR_step_size<mpc_complex>(const Options & opts) {
    return opts.mpfrs.at("nr_step_size")*
        (mpc_complex(1,1))/mp::sqrt(mpfr_float(2));
}

namespace HankelCoefficients {
// Riccati coefficients for even problems defined by:
// [ -d^2/dx^2 + V(x) ] \psi(x) = E \psi(x).
// Where V(x) = V(-x), and therefore \psi(x) is either even or odd.
// s = 0 for even eigenfunctions, and s = 1 for odd eigenfunctions
// V(x) must not be singular at origin, meaning it can be expanded in 
// a Taylor series: V(x) = sum_{j=0}^{\infty} v_j x^{2j}
template <class T>
vector<T> symmetric(
        const int N,
        const int s,
        const T E,
        const vector<T>& v_coefs
        )
{ 
    vector<T> coefs;
    coefs.reserve(N+1);
    coefs.push_back(E - v_coefs[0]);

    T sum;

    for ( int j = 1; j <= N; j++ ) {
        sum = 0;
        for ( int k = 0; k <= j-1; k++ ) {
            sum += coefs[k]*coefs[j-k-1];
        }
        coefs.push_back((sum - v_coefs[j])/(2*j+2*s+1));
    }

    return coefs;
}

// Riccati coefficients for even problems defined by:
// [ -d^2/dx^2 + V(x) ] \psi(x) = E \psi(x).
// No symmetry is assumed for V(x) or \psi(x), and therefore both the initial
// coefficient f[0] = f0 and the energy E are required.
//
template <class T>
vector<T> asymmetric(
        const int N,
        const T E,
        const T f0, 
        const vector<T>& v_coefs
        )
{
    vector<T> coefs;
    coefs.reserve(N+1);
    coefs.push_back(f0);

    T sum;

    for ( int j = 1; j <= N; j++ ) {
        sum = 0;
        for ( int k = 0; k <= j-1; k++ ) {
            sum += coefs[k]*coefs[j-k-1];
        }
        if ( j == 1 ) sum += E;
        coefs.push_back((sum - v_coefs[j-1])/j);

    }
    
    return coefs;
}

} //namespace hankcoefs

// Return a vector with the Hankel coefficients up to N
template <class T>
vector<T> hankcoefs(
        const int N, 
        const Problem &problem, 
        const int s,
        const T E0, 
        const Options & opts
        ) 
{
    vector<T> coefs;
    coefs.reserve(N+1);

    vector<T> v_coefs = problem.get_coefficients<T>(N, opts);
    string problem_type = opts.strings.at("problem_type");

    if ( problem_type == "even" ) {
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
    } else if ( problem_type == "radial" ) {
        T vm1 = problem.get_neg_coeff<T>(opts);
        coefs.push_back(-vm1/(s+1));

        {
            T val;
            for ( int j = 1; j <= N; j++ ) {
                val = 0;
                for ( int k = 0; k <= j-1; k++ ) {
                    val += coefs[k]*coefs[j-k-1];
                }
                val -= 2*v_coefs[j-1];
                if ( j == 1 ) val += 2*E0;

                val = val/(2*s + j + 2);

                coefs.push_back(val);
            }
        }
    }

    return coefs;
}

// Return the hankel determinant with D and d, for Problem problem and 
// parameter s, evaluated at point E0
template <class T>
T hankdet(
        const int D, 
        const int d, 
        // Coefficients f[d+2]...f[2*D+d-2]. This vector is destroyed.
        std::vector<T>& coefs
        ) {
    vector<T> coefsm1, coefsm2;

    if ( D == 0 ) {
        return 1;
    } else {
        coefsm1 = std::vector<T>(2*D, T(1));

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


// Find the root of the Hankel determinant for D, d, s, and problem. Store it 
// in E. Return 0 if succeeded and 1 if the number of iterations was greater
// than opts.ints.at("nr_max_iter").
template <class T>
int RPM_find_root(
        const int D, 
        const int d, 
        const Problem &problem,
        const int s,
        const T E0,
        const Options & opts, 
        T & E
        )
{
    T Eold;
    mpfr_float tol(opts.mpfrs.at("nr_tolerance"));
    mpfr_float desv = tol + 1;
    T h;
    T dH, H;

    E = E0;

    int niter = 0;

    h = NR_step_size<T>(opts);

    while ( desv > tol ) {
        Eold = E;
        H = hankdet<T>(D, d, problem, s, E, opts);
        dH = hankdet<T>(D, d, problem, s, E + h, opts) - 
             hankdet<T>(D, d, problem, s, E - h, opts);
        dH /= 2*h;

        E = Eold - H/dH;
        desv = mp::abs(Eold - E);

        if ( opts.ints.at("log_nr") ) {
            ofstream f(opts.strings.at("log_file"), ios_base::app);
            f.precision(opts.ints.at("digits"));
            f << setw(8) << niter << " " << setw(opts.ints.at("digits") + 10);

            if ( opts.ints.at("use_complex") ) {
                f << left << E.real() << " " << E.imag() << endl;
            } else {
                f << left << E << endl;
            }
            f.close();
        }

        if ( niter++ >= opts.ints.at("nr_max_iter") ) {
            cout << "Maximum number of NR iterations reached." 
                << endl;
            return 1;
        }
    }

    return 0;
}

template <class T>
int RPM_solve(const Problem & problem, const Options & opts, T & E) {
    assign_E0(E, opts);

    T Eold = 0;
    T logdif = 0; 
    int d = opts.ints.at("d");
    int s = opts.ints.at("s");

    int accurate_digits = 0;

    for ( int D = opts.ints.at("Dmin"); D <= opts.ints.at("Dmax"); D++ ) {
        Eold = E;
        int result, newd = d, ntries = 0;

        result = RPM_find_root<T>(D, d, problem, s, Eold, opts, E);

        // If the iterations failed, try again with d+=1
        while ( result == 1 ) {
            cout << "Trying d = " << ++newd << endl;
            result = RPM_find_root<T>(D, newd, problem, s, Eold, opts, E);
            if ( ++ntries == 3 ) {
                cout << "Convergence failed with four different values of d."
                    << " Aborting." << endl;
                return 1;
            }
        }

        logdif = -log10(abs(E - Eold));

        cout << setw(6) << D << " " << setw(opts.ints.at("target_digits") + 10); 
        if ( opts.ints.at("use_complex") ) {
            cout << left << E.real() << " " << E.imag() << endl;
        } else {
            cout << left << E << endl;
        }

        if ( 
                opts.ints.at("infinite_digits") == 0 &&
                static_cast<int>(logdif) >= opts.ints.at("target_digits") 
            ) 
            break;

        // Output the same to the log file
        if ( opts.ints.at("log_nr") ) {
            ofstream f(opts.strings.at("log_file"), ios_base::app);
            f.precision(opts.ints.at("digits"));
            f << endl;
            f << setw(8) << "D = " + to_string(D) << " " << setw(opts.ints.at("digits") + 10);

            if ( opts.ints.at("use_complex") ) {
                f << left << E.real() << " " << E.imag() << endl << endl;
            } else {
                f << left << E << endl << endl;
            }
            f.close();
        }
    }

    return 0;
}

#endif
