// Std includes
#include <cmath>
#include <iostream>
// Thirdparties includes
#include <Eigen/Dense>
// Lib includes
#include "fl0w/jhtdb.h"

const unsigned int DIM = 3;

using tScalar = double;
using tSpaceVector = Eigen::Matrix<tScalar, DIM, 1>;
using tSpaceMatrix = Eigen::Matrix<tScalar, DIM, DIM>;
template<typename... Args>
using tView = Eigen::Map<Args...>;

int main () {
	// init
	const tSpaceVector x({0.0, 0.0, 0.0});
	const double t = 0.0;
	// uniform
	std::cout << "\n";
	std::cout << "jhtdb\n";
	std::cout << "\n";
    std::cout << "velocity( (" << x.transpose() << ")T, " << t << " )T -> " << fl0w::Jhtdb<tSpaceVector, tSpaceMatrix, tView>::Isotropic::getVelocity(x.data(), t).transpose() << "\n";
    std::cout << "velocityGradients( (" << x.transpose() << ")T, " << t << " ) -> \n" << fl0w::Jhtdb<tSpaceVector, tSpaceMatrix, tView>::Isotropic::getVelocityGradients(x.data(), t) << "\n";
    std::cout << std::endl;
}
