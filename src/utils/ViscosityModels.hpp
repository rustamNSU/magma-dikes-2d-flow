#pragma once
#include <string>
#include <cmath>
#include <array>


/**
 * @brief Vogel–Fulcher–Tammann equation for temperature dependent viscosity with constant parameters
 * 
 * @param T [K]
 * @param A 
 * @param B 
 * @param C [K]
 * @return double 
 */
inline double VFT_constant_viscosity(double T, double A, double B, double C){
    constexpr double ln10 = 2.302585092994045684;
    return std::exp(ln10 * (A + B / (T - C)));
}


class GiordanoViscosity{
    private:
        double sio2;
        double tio2;
        double al2o3;
        double feot; // FeO(T)
        double mno;  // MnO
        double mgo;  // MgO
        double cao;  // CaO
        double na2o; // Na2O
        double k2o;  // K2O
        double p2o5; // P2O5
        double f2o1; // F2O-1

        std::array<double, 11> composition; // Without h2o
        std::array<double, 11> mw_compositions = {
            60.0843, 79.8658, 101.961276, 71.8444, 70.937449, 40.3044, 56.0774, 61.97894, 94.1960, 141.9446, 18.9984
        }; // Molar weights of oxides (g/mol)
        std::array<double, 11> base_mole;
        double mw_h2o = 18.01528; // Molecular weight of h2o (g/mol)

        std::array<double, 10> bcoef = {
            159.56, -173.34, 72.13, 75.69, -38.98, -84.08, 141.54, -2.43, -0.91, 17.62
        }; // b1, ..., b7, b11, b12, b13
        std::array<double, 7> ccoef = {
            2.75, 15.72, 8.32, 10.2, -12.29, -99.54, 0.3
        }; // c1, ..., c6, b11


    public:
        GiordanoViscosity() = default;
        /**
         * @brief Set the Composition object
         * 
         * @param composition_ oxides without h2o in %
         */
        void setComposition(const std::array<double, 11>& composition_);


        /**
         * @brief 
         * @param wt_h2o weight fraction of water in the melt (not % !!!)
         * @param temperature melt temperature [C]
         * @return 
         */
        double calculateViscosity(double wt_h2o, double temperature) const;
};