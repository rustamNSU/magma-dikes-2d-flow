#include "Mesh.hpp"
#include "Elasticity.hpp"
#include "DikeData.hpp"
#include <iostream>

int main(){
    Mesh mesh(10, -1.0, 1.0);
    Elasticity elasticity(1.0, 0.2, &mesh);
    DikeData dike(&mesh);
    std::cout << mesh.getx() << std::endl;
    std::cout << mesh.getxl() << std::endl;
    std::cout << mesh.getxr() << std::endl;
    std::cout << elasticity.get_matrix() << std::endl;
    return 0;
}