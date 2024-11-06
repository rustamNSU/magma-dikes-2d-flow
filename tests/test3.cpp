#include <iostream>
#include <cmath>
#include <vector>
#include <tuple>
#include <fmt/ranges.h>
#include "uniform_interp.hpp"
#include <iostream>



std::vector<double> generateX(int N, double xmin = 0.0, double xmax = 1.0){
    double dx = (xmax - xmin) / (N-1);
    std::vector<double> x(N);
    for (int i = 0; i < N; i++){
        x[i] = xmin + i*dx;
    }
    return x;
}


template<typename Func>
std::vector<double> generateF(Func f, const std::vector<double>& x){
    std::vector<double> y(x.size());
    for (int i = 0; i < x.size(); i++){
        y[i] = f(x[i]);
    }
    return y;
}


template<typename Func>
std::vector<double> generateF(Func func, const std::vector<double>& x, const std::vector<double>& y){
    std::vector<double> f(x.size()*y.size());
    int ny = y.size();
    for (int ix = 0; ix < x.size(); ix++){
        for (int iy = 0; iy < y.size(); iy++){
            f[ix*ny + iy] = func(x[ix], y[iy]);
        }
    }
    return f;
}


template<typename Func>
std::vector<double> generateF(
    Func func, 
    const std::vector<double>& x, 
    const std::vector<double>& y, 
    const std::vector<double>& z
){
    std::vector<double> f(x.size()*y.size()*z.size());
    int ny = y.size();
    int nz = z.size();
    for (int ix = 0; ix < x.size(); ix++){
        for (int iy = 0; iy < y.size(); iy++){
            for (int iz = 0; iz < z.size(); iz++){
                f[(ix*ny + iy)*nz + iz] = func(x[ix], y[iy], z[iz]);
            }
        }
    }
    return f;
}


int main(int argc, char ** argv){
    std::vector<double> X, Y, Z, F, xx, yy, zz, ff;
    auto func1d = [](double x){return std::sin(2 * M_PI * x);};
    X = generateX(11, 0.0, 1.0);
    F = generateF(func1d, X);
    UniformInterpolation1d interp1d(X, F);

    xx = generateX(29, -0.1, 1.1);
    ff.resize(xx.size());
    for (int i = 0; i < xx.size(); i++){
        ff[i] = interp1d.getValue(xx[i]);
    }
    fmt::print("xx = {}\n", xx);
    fmt::print("ff = {}\n", ff);
    fmt::print("\n\n\n2d-interp\n");

    /* Bilinear */
    X = generateX(11, -1.0, 1.0);
    Y = generateX(11, -1.0, 1.0);
    auto func2d = [](double x, double y){return x*x+y*y*y;};
    F = generateF(func2d, X, Y);
    UniformInterpolation2d interp2d(X, Y, F);
    xx = generateX(33, -1.1, 1.1);
    yy = generateX(19, -1.1, 1.1);
    ff.resize(xx.size()*yy.size());
    int id = 0;
    for (int ix = 0; ix < xx.size(); ix++){
        for (int iy = 0; iy < yy.size(); iy++){
            ff[id] = interp2d.getValue(xx[ix], yy[iy]);
            id++;
        }
    }
    fmt::print("xx = {}\n", xx);
    fmt::print("yy = {}\n", yy);
    fmt::print("ff = {}\n", ff);

    fmt::print("\n\n\n3d-interp\n");
    /* Trilinear */
    X = generateX(101, 0.0, 2.0);
    Y = generateX(101, 0.0, 1.0);
    Z = generateX(101, -1.0, 1.0);
    auto func3d = [](double x, double y, double z){return x + y*y + z*z*z;};
    F = generateF(func3d, X, Y, Z);
    UniformInterpolation3d interp3d(X, Y, Z, F);
    xx = generateX(200, 0.0, 2.0);
    yy = generateX(200, 0.0, 1.0);
    zz = generateX(200, -1.0, 1.0);
    id = 0;
    double error = 0;
    for (int ix = 0; ix < xx.size(); ix++){
        for (int iy = 0; iy < yy.size(); iy++){
            for (int iz = 0; iz < zz.size(); iz++){
                double err = std::abs(interp3d.getValue(xx[ix], yy[iy], zz[iz]) - func3d(xx[ix], yy[iy], zz[iz]));
                error = err > error ? err : error;
                id++;
            }
        }
    }
    fmt::print("error = {}\n", error);
    return 0;
}