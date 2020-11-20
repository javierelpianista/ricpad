#include <ricpad.hpp>

mpfr_float gi_to_mpfr(
        const gi::numeric num
        ) 
{
    ostringstream oss;

    oss << num;

    return mpfr_float(oss.str());
}

vector<mpfr_float> get_coefficients(
        const gi::ex Q, 
        const gi::ex &x, 
        const int N, 
        const bool use_rationals
        ) 
{
    vector<mpfr_float> data;

    gi::numeric to_push;
    gi::ex Q_series = gi::series_to_poly(Q.series(x == 0, N+3));
    gi::ex tmp;
    mpfr_float to_push_n;

    int mult;

    if ( options["problem_type"] == "even" ) mult = 2;

    for ( int i = 0; i <= N; i++ ) {
        if ( use_rationals ) {
            to_push = gi::ex_to<gi::numeric>(Q_series.coeff(x, mult*i));
        } else {
            tmp = Q_series.coeff(x, mult*i).evalf();
            to_push = gi::ex_to<gi::numeric>(Q_series.coeff(x, mult*i).evalf());
        }
        to_push_n = gi_to_mpfr(to_push);
        data.push_back(to_push_n);
    }

    return data;
}

// Return a vector with the Hankel coefficients up to N
vector<mpfr_float> hankcoefs(
        const int N, 
        const vector<mpfr_float> &V, 
        const int s,
        const mpfr_float E0
        ) 
{
    vector<mpfr_float> coefs;
    coefs.reserve(N+1);
    coefs.push_back(E0 - V[0]);

    {
        mpfr_float sum;
        for ( int j = 1; j <= N; j++ ) {
            sum = 0;
            for ( int k = 0; k <= j-1; k++ ) {
                sum += coefs[k]*coefs[j-k-1];
            }
            coefs.push_back((sum - V[j])/(2*j+s+1));
        }
    }

    return coefs;
}

mpfr_float hankdet(
        const int D, 
        const int d, 
        const vector<mpfr_float> &V, 
        const int s,
        const mpfr_float E0) {
    vector<mpfr_float> coefs, coefsm1, coefsm2;

    if ( D == 0 ) {
        return 1;
    } else {
        coefsm1 = std::vector<mpfr_float>(2*D, mpfr_float(1));

        //for ( int i = d+1; i <= 2*D + d - 1; i++ ) 
        //    coefs.emplace_back(hankcoef(i, V, s, E0));

        coefs = hankcoefs(2*D + d - 1, V, s, E0);
        coefs.erase(coefs.begin(), coefs.begin()+d+1);

        if ( D == 1 ) return coefs[0];

        for ( int j = 2; j <= D; j++ ) {
            coefsm2 = (vector<mpfr_float>&&)(coefsm1);
            coefsm1 = (vector<mpfr_float>&&)(coefs);

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

mpfr_float RPM_solve(
        const int D, 
        const int d, 
        const vector<mpfr_float> &V, 
        const int s,
        const mpfr_float E0 
        )
{
    mpfr_float E = E0, Eold;
    mpfr_float tol(options["tolerance"]);
    mpfr_float desv = tol + 1;
    mpfr_float h(options["diff_epsilon"]);
    mpfr_float dH, H;

    while ( desv > tol ) {
        Eold = E;
        H = hankdet(D, d, V, s, E);
        dH = (hankdet(D, d, V, s, E + h) - hankdet(D, d, V, s, E - h));
        dH /= 2*h;

        E = Eold - H/dH;
        desv = mp::abs(Eold - E);
    }

    return E;
}

