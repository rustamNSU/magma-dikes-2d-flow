#include "ViscosityModels.hpp"
#include <algorithm>
#include <numeric>


void GiordanoViscosity::setComposition(const std::array<double, 11>& composition_){
    this->composition = composition_;

    /* Normalize composition */
    double sum = std::accumulate(composition.begin(), composition.end(), 0.0);
    for (int i = 0; i < composition.size(); ++i){
        composition[i] = composition[i] / sum;
        base_mole[i] = composition[i] / mw_compositions[i];
    }
}


double GiordanoViscosity::calculateViscosity(double wt_h2o, double temperature) const{
    constexpr double A = -4.55;
    
    /* Mole fraction */
    double mp_h2o = wt_h2o / mw_h2o;
    auto wtn = base_mole;
    double sum = 0;
    for (auto &e : wtn){
        e = e * (1.0 - wt_h2o);
        sum += e;
    }
    sum += mp_h2o;
    double xsio2  = 100.0 * wtn[0] / sum;
    double xtio2  = 100.0 * wtn[1] / sum;
    double xal2o3 = 100.0 * wtn[2] / sum;
    double xfeo   = 100.0 * wtn[3] / sum;
    double xmno   = 100.0 * wtn[4] / sum;
    double xmgo   = 100.0 * wtn[5] / sum;
    double xcao   = 100.0 * wtn[6] / sum;
    double xna2o  = 100.0 * wtn[7] / sum;
    double xk2o   = 100.0 * wtn[8] / sum;
    double xp2o5  = 100.0 * wtn[9] / sum;
    double xf2o1  = 100.0 * wtn[10] / sum;
    double xh2o   = 100.0 * mp_h2o / sum;
    std::array<double, 10> bmf;
    std::array<double, 7> cmf;
    bmf[0] = xsio2 + xtio2; // b1  
    bmf[1] = xal2o3; // b2  
    bmf[2] = xfeo + xmno + xp2o5; // b3  
    bmf[3] = xmgo; // b4  
    bmf[4] = xcao; // b5  
    bmf[5] = xna2o + xh2o + xf2o1; // b6  
    bmf[6] = xh2o + xf2o1 + std::log(1.0 + xh2o); // b7  
    bmf[7] = (xsio2 + xtio2) * (xfeo + xmno + xmgo); // b11 
    bmf[8] = (xsio2 + xtio2 + xal2o3 + xp2o5) * (xna2o + xk2o + xh2o); // b12 
    bmf[9] = (xal2o3) * (xna2o + xk2o);  // b13 
    
    cmf[0] = xsio2; // b1  
    cmf[1] = xtio2 + xal2o3; // b2  
    cmf[2] = xfeo + xmno + xmgo; // b3  
    cmf[3] = xcao; // b4
    cmf[4] = xna2o + xk2o; // b5  
    cmf[5] = std::log(1.0 + xh2o + xf2o1); // b6  
    cmf[6] = (xal2o3 + xfeo + xmno + xmgo + xcao - xp2o5) * (xna2o + xk2o + xh2o + xf2o1); // b7  

    double B = std::inner_product(bcoef.begin(), bcoef.end(), bmf.begin(), 0.0);
    double C = std::inner_product(ccoef.begin(), ccoef.end(), cmf.begin(), 0.0);
    double Tg = B / (12.0 - A) + C; // Glass transition temperature
    double F = B / (Tg * (1.0 - C / Tg) * (1.0 - C / Tg)); // The fragility of silicate melts
    double T = temperature + 273.15;
    double viscosity = std::pow(10.0, A + B / (T - C));
    return viscosity;
}