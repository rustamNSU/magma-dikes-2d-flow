#pragma once
#include <cmath>

inline double liquidus_temperature(double p){
    constexpr double al = 1205.7; 
    constexpr double bl = 6.0;    // bl = 6 for p in kbar (1e8 Pa)
    constexpr double cl = 285.7;
    constexpr double dl = 200.0;
    constexpr double el = 0.7;    // el = 0.7 for p in kbar (1e8 Pa)
    constexpr double fl = 11.0;   // fl = 6 for p in kbar (1e8 Pa)
    constexpr double xh2od = 1.0;
    double pkbar = p * 1e-8;
    return al + bl*pkbar - xh2od*(cl + fl*pkbar - dl/(pkbar + el));
}


inline double solidus_temperature(double p){
    constexpr double as = 854.0896; 
    constexpr double bs = 6.0;     // bs =6 for p in kbar (1e8 Pa)
    constexpr double cs = 224.0896;
    constexpr double ds = 80.0;
    constexpr double es = 0.357;   // es for p in kbar (1e8 Pa)
    constexpr double fs = 6.0;     // fs for p in kbar (1e8 Pa)
    constexpr double xh2od = 1.0;
    double pkbar = p * 1e-8;
    return as + bs*pkbar - xh2od*(cs + fs*pkbar - ds/(pkbar + es));
}


inline double beta_equilibrium(double p, double T){
    double aF = -4.974;
    double bF = 28.623;
    double cF = -52.708;
    double dF = 34.816;
    double Tl = liquidus_temperature(p);
    double Ts = solidus_temperature(p);
    double x = (T-Ts) / (Tl - Ts);
    return 1.0 / (1.0 + std::exp(aF + bF*x + cF*x*x + dF*x*x*x));
}