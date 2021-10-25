// Std includes
#include <vector>
#include <iostream>
// Thirdparties includes
#include <Eigen/Dense>
// Lib includes
#include "fl0w/kinematic.h"

const unsigned int DIM = 3;

using TypeScalar = double;
using TypeVector = Eigen::Matrix<TypeScalar, DIM, 1>;
using TypeMatrix = Eigen::Matrix<TypeScalar, DIM, DIM>;
template<typename... Args>
using TypeRef = Eigen::Ref<Args...>;
using TypeFlow = fl0w::kinematic::Kinematic<TypeVector, TypeMatrix, TypeRef, std::vector>;

void print(const TypeFlow& flow, const TypeVector& x, const TypeScalar& t) {
    std::cout << std::endl;
    std::cout << "flow.getVelocity(" << x.transpose() << ", " << t << ") -> " << flow.getVelocity(x, t).transpose() << std::endl;
    std::cout << "flow.getVorticity(" << x.transpose() << ", " << t << ") -> " << flow.getVorticity(x, t).transpose() << std::endl;
    std::cout << "flow.getAcceleration(" << x.transpose() << ", " << t << ") -> " << flow.getAcceleration(x, t).transpose() << std::endl;
    std::cout << std::endl;
}

int main () { 
    TypeFlow flow;
    flow.init();
    std::cout << "k " << (*flow.sK)[0] << std::endl;
    std::cout << "k " << (*flow.sK)[1] << std::endl;
    std::cout << "k " << (*flow.sK)[2] << std::endl;
    std::cout << "k " << (*flow.sK)[3] << std::endl;
    std::cout << "n0 " << (*flow.sK)[0].norm() << std::endl;
    std::cout << "nl " << (*flow.sK)[(*flow.sK).size()-1].norm() << std::endl;
    std::cout << "a " << flow.a[0] << std::endl;
    std::cout << "b " << flow.b[0] << std::endl;
    std::cout << "o " << flow.omega[0] << std::endl;
    TypeVector x;
    double t;
    // Init
    x << 0.0, 0.0, 0.0;
    t = 0.0;
    print(flow, x, t);
    //x << 10.0, -10.0, 40.0;
    //t = 10.0;
    //print(flow, x, t);
    //x << 1.0, -3.0, -1000.0;
    //t = 40.0;
    //print(flow, x, t);
}
