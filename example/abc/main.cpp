// Std includes
#include <cmath>
#include <iostream>
// Thirdparties includes
#include <Eigen/Dense>
// Lib includes
#include "fl0w/abc.h"

const unsigned int DIM = 3;

using TypeScalar = double;
using TypeVector = Eigen::Matrix<TypeScalar, DIM, 1>;
using TypeMatrix = Eigen::Matrix<TypeScalar, DIM, DIM>;
template<typename... Args>
using TypeRef = Eigen::Ref<Args...>;
using TypeFlow = fl0w::ABC<TypeVector, TypeMatrix, TypeRef>;

void print(const TypeFlow& flow, const TypeVector& x, const TypeScalar& t) {
    std::cout << std::endl;
    std::cout << "flow.getVelocity(" << x.transpose() << ", " << t << ") -> " << flow.getVelocity(x, t).transpose() << std::endl;
    std::cout << "flow.getVorticity(" << x.transpose() << ", " << t << ") -> " << flow.getVorticity(x, t).transpose() << std::endl;
    std::cout << "flow.getAcceleration(" << x.transpose() << ", " << t << ") -> " << flow.getAcceleration(x, t).transpose() << std::endl;
    std::cout << std::endl;
}

int main () { 
    TypeFlow flow;
    TypeVector x;
    double t;
    // Init
    x << 0.0, 0.0, 0.0;
    t = 0.0;
    print(flow, x, t);
}
