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


/* Correction for beta change in magma assembling (beta - betach)*/
inline double phi_chamber(double beta, double betach = 0.45){
    double dbeta = beta - betach;
    return std::exp(4.33*dbeta + 10.48*dbeta*dbeta);
}


inline double theta_coef(double beta, double betach = 0.45){
    double phi = beta < 0.25 ? 0.63 : phi_chamber(beta, betach);
    constexpr double eps = 0.999916;
    constexpr double bstar = 0.673;
    constexpr double gamma = 3.98937;
    constexpr double delta = 16.9386;
    // double alpha = 0.5*std::sqrt(M_PI)/eps/bstar;
    constexpr double alpha = 0.88622692545/eps/bstar;
    double fbx = std::erf(alpha * beta * (1.0 + std::pow(beta / bstar, gamma)));
    double theta = phi * (1 + std::pow(beta / bstar, delta)) * std::pow(1.0 - eps*fbx, -2.5*bstar);
    return std::max(1.0, theta);
}
// function theta=Theta2(bx)
//     be0=0.45;
//     teta1=exp(4.33.*(bx-be0)+10.48.*(bx-be0).*(bx-be0));
//     teta1(bx<0.25)=0.63;
//     % teta1(:) = 1.0;
//     eps= 0.999916;
//     bp = 0.673;
//     gamma = 3.98937;
//     delta = 16.9386;
//     Fbx=erf((sqrt(pi)/2/eps/bp).*bx.*(1+(bx./bp).^gamma));
//     teta2 = (1+(bx./bp).^delta).*(1-eps.*Fbx).^(-2.5*bp);
//     theta=max(1.,teta1.*teta2);
// end