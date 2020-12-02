#include <iostream>

#include <boost/multiprecision/mpfr.hpp>
#include <boost/multiprecision/mpc.hpp>
#include <boost/multiprecision/eigen.hpp>
#include <Eigen/Dense>

#include <problem.hpp>
#include <ricpad.hpp>
#include <options.hpp>
#include <read_input.hpp>

using namespace std;
using boost::multiprecision::mpfr_float;
using boost::multiprecision::mpc_complex;
using C = mpc_complex;

namespace mp = boost::multiprecision;

//std::function<mpfr_float(mpfr_float&)> test = [](mpfr_float & x){ return x*x; };


/*
template <typename C, typename R>
void NR_solve( 
        std::vector<std::function<C(C&)>> f, 
        std::vector<C> &x, 
        R tol, 
        C h 
        ) 
{
    R desv = tol + 1;
    C xold, fval, dfval, xp, xm;
    int niter, max_iter = 100;

    while ( desv > tol ) {
        xold = x[0];
        fval = f[0](x[0]);
        xp = x[0] + h;
        xm = x[0] - h;

        dfval = (f[0](xp) - f[0](xm))/(2*h);

        x[0] = xold - fval/dfval;
        desv = mp::abs(x[0] - xold);

        if ( niter++ >= max_iter ) {
            cout << "Maximum number of NR iterations reached." << endl;
            exit(1);
        }
    }
}
*/

template <
    typename C, // complex number type
    typename R, // real number type
    int N       // Number of equations and unknowns
>

int NR_solve( 
        std::vector<std::function<C(vector<C>&)>> f, 
        std::vector<C>& x, 
        R tol, 
        C h 
        ) 
{
    Eigen::Matrix<C, N, N> jacobian, inv_jacobian;

    if ( f.size() != N ) {
        cout << "The size of f is wrong." << endl;
        return 1;
    } else if ( x.size() != N ) {
        cout << "The size of x is wrong." << endl;
    }

    vector<C> xp(x), xm(x), row;

    for ( int i = 0; i < x.size(); i++ ) {
        for ( int j = 0; j < x.size(); j++ ) {
            xp[j] = x[j] + h;
            xm[j] = x[j] - h;

            C val = (f[i](xp) - f[i](xm))/(2*h);

            xp[j] = x[j];
            xm[j] = x[j];

            jacobian(i, j) = move(val);
        }
    }

    inv_jacobian = jacobian.inverse();

    for ( int i = 0; i < N; i++ ) {
        for ( int j = 0; j < N; j++ ) {
            cout << jacobian(i, j).real() << " ";
        }
        cout << endl;
    }

    for ( int i = 0; i < N; i++ ) {
        for ( int j = 0; j < N; j++ ) {
            cout << inv_jacobian(i, j).real() << " ";
        }
        cout << endl;
    }


    return 0;
}


int main(int argc, char * argv[]) {
    mpfr_float::default_precision(1000);
    mpc_complex::default_precision(1000);

    mpfr_float x = 4, tol("1E-800");
    mpc_complex z0(4), h("(1E-900,1E-900)");

    std::vector<mpc_complex> z;
    z.emplace_back(mpc_complex(4));
    z.emplace_back(mpc_complex(6));

    std::vector<std::function<C(vector<C>&)>> f;

    std::function<C(vector<C>&)> test1 = 
        [](vector<C>& x){ return x[0]*x[0] + x[1] - 2; };

    std::function<C(vector<C>&)> test2 = 
        [](vector<C>& x){ return 2*x[1]*x[1] -2*x[0] - 3; };

    f.push_back(move(test1));
    f.push_back(move(test2));

    NR_solve<mpc_complex, mpfr_float, 2>(f, z, tol, h);

    return 0;
}
