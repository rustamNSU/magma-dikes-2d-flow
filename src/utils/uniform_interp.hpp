#pragma once
#include <vector>

class UniformInterpolation1d{
    private:
        std::vector<double> X;
        std::vector<double> F;
        double epsilon = 1e-10;
    
    public:
        UniformInterpolation1d(const std::vector<double>& X, const std::vector<double>& F);
        double getValue(double x) const;
};


class UniformInterpolation2d{
    private:
        std::vector<double> X;
        std::vector<double> Y;
        std::vector<double> F;
        int nx;
        int ny;
        double epsilon = 1e-10;
        inline int ID(int ix, int iy) const{
            return ix*ny + iy;
        }
    
    public:
        UniformInterpolation2d(
            const std::vector<double>& X, 
            const std::vector<double>& Y, 
            const std::vector<double>& F
        );
        double getValue(double x, double y) const;
};


class UniformInterpolation3d{
    private:
        std::vector<double> X;
        std::vector<double> Y;
        std::vector<double> Z;
        std::vector<double> F;
        int nx;
        int ny;
        int nz;
        inline int ID(int ix, int iy, int iz) const{
            return (ix*ny + iy)*nz + iz;
        }
        double epsilon = 1e-10;
    
    public:
        UniformInterpolation3d(
            const std::vector<double>& X, 
            const std::vector<double>& Y, 
            const std::vector<double>& Z, 
            const std::vector<double>& F
        );
        double getValue(double x, double y, double z) const;
};
