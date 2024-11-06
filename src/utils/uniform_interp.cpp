#include "uniform_interp.hpp"
#include <stdexcept>
#include <cmath>
#include <tuple>


inline bool is_uniform_data(const std::vector<double>& x){

}


inline std::tuple<int, double> get_id(double x, const std::vector<double>& X){
    int N = X.size();
    int id = (int)std::floor((x - X[0]) / (X[1] - X[0]));
    if (id < 0) return std::make_tuple(0, X[0]);
    if (id > X.size() - 2) return std::make_tuple(static_cast<int>(X.size() - 2), X[N-1]);
    return std::make_tuple(id, x);
}


inline double linear_interp(
    double x, 
    double x1, 
    double x2,
    double y1,
    double y2
){
    double wt2 = (x - x1) / (x2 - x1);
    double wt1 = 1.0 - wt2;
    return wt1*y1 + wt2*y2;
}


UniformInterpolation1d::UniformInterpolation1d(
    const std::vector<double>& X, 
    const std::vector<double>& F
) : X(X), F(F)
{
    if (X.size() != F.size()){
        throw std::invalid_argument("x and y must be the same size!\n");
    }
    double dx = X[1] - X[0];
    for (int i = 0; i < X.size()-1; i++){
        double dxx = X[i+1] - X[i];
        if (std::abs(dxx - dx) > epsilon){
            throw std::invalid_argument("x must be uniform spaced!\n");
        }
    }
}


double UniformInterpolation1d::getValue(double x) const{
    auto [id, xx] = get_id(x, X);
    return linear_interp(xx, X[id], X[id+1], F[id], F[id+1]);
}


UniformInterpolation2d::UniformInterpolation2d(
    const std::vector<double>& X, 
    const std::vector<double>& Y, 
    const std::vector<double>& F
) : X(X), Y(Y), F(F)
{
    if (X.size() * Y.size() != F.size()){
        throw std::invalid_argument("x * y must be the same size as F!\n");
    }
    nx = X.size();
    ny = Y.size();
}


double UniformInterpolation2d::getValue(double x, double y) const{
    auto [idx, xx] = get_id(x, X);
    auto [idy, yy] = get_id(y, Y);
    
    double f1 = linear_interp(xx, X[idx], X[idx+1], F[ID(idx, idy)], F[ID(idx+1, idy)]);
    double f2 = linear_interp(xx, X[idx], X[idx+1], F[ID(idx, idy+1)], F[ID(idx+1, idy+1)]);
    return linear_interp(yy, Y[idy], Y[idy+1], f1, f2);
}


UniformInterpolation3d::UniformInterpolation3d(
    const std::vector<double>& X, 
    const std::vector<double>& Y, 
    const std::vector<double>& Z, 
    const std::vector<double>& F
) : X(X), Y(Y), Z(Z), F(F){
    if (X.size() * Y.size() * Z.size() != F.size()){
        throw std::invalid_argument("x * y must be the same size as F!\n");
    }
    nx = X.size();
    ny = Y.size();
    nz = Z.size();
}


double UniformInterpolation3d::getValue(double x, double y, double z) const{
    auto [idx, xx] = get_id(x, X);
    auto [idy, yy] = get_id(y, Y);
    auto [idz, zz] = get_id(z, Z);
    
    double f11 = linear_interp(xx, X[idx], X[idx+1], F[ID(idx, idy, idz)], F[ID(idx+1, idy, idz)]);
    double f21 = linear_interp(xx, X[idx], X[idx+1], F[ID(idx, idy+1, idz)], F[ID(idx+1, idy+1, idz)]);
    double f12 = linear_interp(xx, X[idx], X[idx+1], F[ID(idx, idy, idz+1)], F[ID(idx+1, idy, idz+1)]);
    double f22 = linear_interp(xx, X[idx], X[idx+1], F[ID(idx, idy+1, idz+1)], F[ID(idx+1, idy+1, idz+1)]);
    double f1 = linear_interp(yy, Y[idy], Y[idy+1], f11, f21);
    double f2 = linear_interp(yy, Y[idy], Y[idy+1], f12, f22);
    return linear_interp(zz, Z[idz], Z[idz+1], f1, f2);
}