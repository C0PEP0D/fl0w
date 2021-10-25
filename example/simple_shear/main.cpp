// Std includes
#include <cmath>
#include <iostream>
// Thirdparties includes
#include <Eigen/Dense>
// Lib includes
#include "fl0w/simple_shear.h"

const unsigned int DIM = 3;

using TypeScalar = double;
using TypeVector = Eigen::Matrix<TypeScalar, DIM, 1>;
using TypeMatrix = Eigen::Matrix<TypeScalar, DIM, DIM>;
template<typename... Args>
using TypeRef = Eigen::Ref<Args...>;
using TypeFlow = fl0w::SimpleShear<TypeVector, TypeMatrix, TypeRef>;

void test(const TypeFlow& flow, const TypeVector& x, const TypeScalar& t) {
    // Expected values
    TypeVector xW;
    xW << 0.0,
          0.0,
         -0.5;
    TypeMatrix xS;
    xS << 0.0, 0.5, 0.0,
          0.5, 0.0, 0.0,
          0.0, 0.0, 0.0;
    // Output
    TypeVector w = flow.getVorticity(x, t);
    TypeMatrix S = flow.getStrain(x, t);
    std::cout << std::endl;
    std::cout << "flow.getVorticity(" << x.transpose() << ", " << t << ") = " << std::endl; 
    std::cout << w << std::endl;
    std::cout << "flow.getStrain(" << x.transpose() << ", " << t << ") = " << std::endl; 
    std::cout << S << std::endl;
    // Test Output
    std::cout << std::endl;
    std::cout << "Test vorticity succeeded : " << (w == xW) << std::endl; 
    std::cout << "Test strain succeeded : " << (S == xS) << std::endl; 
    std::cout << std::endl;
    // Small temporary additional computations
    std::cout << "flow.getJacobian(" << x.transpose() << ", " << t << ") = " << std::endl; 
    std::cout << flow.getJacobian(x, t) << std::endl;
    TypeVector dir;
    dir << 1.0,
           0.0,
           0.0;
    std::cout << "flow.getJacobian(" << x.transpose() << ", " << t << ") * " << dir.transpose() << " = " << std::endl; 
    std::cout << flow.getJacobian(x, t).transpose() * dir  << std::endl;
}


void print(const TypeFlow& flow, const TypeVector& x, const TypeScalar& t) {
    std::cout << std::endl;
    std::cout << "flow.getVelocity(" << x.transpose() << ", " << t << ") -> " << flow.getVelocity(x, t).transpose() << std::endl;
    std::cout << "flow.getVorticity(" << x.transpose() << ", " << t << ") -> " << flow.getVorticity(x, t).transpose() << std::endl;
    std::cout << "flow.getAcceleration(" << x.transpose() << ", " << t << ") -> " << flow.getAcceleration(x, t).transpose() << std::endl;
    std::cout << std::endl;
}

int main () { 
    TypeFlow flow;
    flow.create(1.0);
    TypeVector x;
    double t;
    // Init
    x << 0.0, 0.0, 0.0;
    t = 0.0;
    test(flow, x, t);
}
