// Std includes
#include <cmath>
#include <iostream>
// Thirdparties includes
#include <Eigen/Dense>
// Lib includes
#include "fl0w/analytic/simple_shear.h"
#include "fl0w/analytic/taylor_green_vortex.h"
#include "fl0w/analytic/uniform.h"

const unsigned int DIM = 2;

using tScalar = double;
using tSpaceVector = Eigen::Matrix<tScalar, DIM, 1>;
using tSpaceMatrix = Eigen::Matrix<tScalar, DIM, DIM>;
template<typename... Args>
using tView = Eigen::Map<Args...>;

int main () {
	// init
	const tSpaceVector x({0.0, 0.0});
	const double t = 0.0;
	// uniform
	std::cout << "\n";
	std::cout << "uniform\n";
	std::cout << "\n";
    std::cout << "velocity( (" << x.transpose() << ")T, " << t << " )T -> " << fl0w::analytic::FlowUniform<tSpaceVector, tSpaceMatrix, tView, 0.0, 0.0>::getVelocity(x.data(), t).transpose() << "\n";
    std::cout << "velocityGradients( (" << x.transpose() << ")T, " << t << " ) -> \n" << fl0w::analytic::FlowUniform<tSpaceVector, tSpaceMatrix, tView, 0.0, 0.0>::getVelocityGradients(x.data(), t) << "\n";
    std::cout << "acceleration( (" << x.transpose() << ")T, " << t << " )T -> " << fl0w::analytic::FlowUniform<tSpaceVector, tSpaceMatrix, tView, 0.0, 0.0>::getAcceleration(x.data(), t).transpose() << "\n";
    std::cout << std::endl;
    // taylor green vortex
   	std::cout << "\n";
   	std::cout << "taylor green vortex\n";
   	std::cout << "\n";
    std::cout << "velocity(" << x.transpose() << ", " << t << " )T -> " << fl0w::analytic::FlowTaylorGreenVortex<tSpaceVector, tSpaceMatrix, tView>::getVelocity(x.data(), t).transpose() << "\n";
    std::cout << "velocityGradients(" << x.transpose() << ", " << t << " ) -> \n" << fl0w::analytic::FlowTaylorGreenVortex<tSpaceVector, tSpaceMatrix, tView>::getVelocityGradients(x.data(), t) << "\n";
    std::cout << "acceleration( (" << x.transpose() << ")T, " << t << " )T -> " << fl0w::analytic::FlowTaylorGreenVortex<tSpaceVector, tSpaceMatrix, tView>::getAcceleration(x.data(), t).transpose() << "\n";
    std::cout << std::endl;
    // simple shear
   	std::cout << "\n";
   	std::cout << "simple shear\n";
   	std::cout << "\n";
    std::cout << "velocity(" << x.transpose() << ", " << t << " )T -> " << fl0w::analytic::FlowSimpleShear<tSpaceVector, tSpaceMatrix, tView, 1.0>::getVelocity(x.data(), t).transpose() << "\n";
    std::cout << "velocityGradients(" << x.transpose() << ", " << t << " ) -> \n" << fl0w::analytic::FlowSimpleShear<tSpaceVector, tSpaceMatrix, tView, 1.0>::getVelocityGradients(x.data(), t) << "\n";
    std::cout << "acceleration( (" << x.transpose() << ")T, " << t << " )T -> " << fl0w::analytic::FlowSimpleShear<tSpaceVector, tSpaceMatrix, tView, 1.0>::getAcceleration(x.data(), t).transpose() << "\n";
    std::cout << std::endl;
}
