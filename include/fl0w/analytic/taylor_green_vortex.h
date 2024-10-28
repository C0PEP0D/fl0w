#ifndef FL0W_TAYLOR_GREEN_VORTEX_H
#define FL0W_TAYLOR_GREEN_VORTEX_H
#pragma once

// std includes
#include <cmath> // sin, cos
// module includes
#include "fl0w/flow.h"

namespace fl0w {

namespace analytic {

template<typename _tSpaceVector, typename _tSpaceMatrix, template<typename...> class _tView>
class FlowTaylorGreenVortex : public Flow<_tSpaceVector, _tSpaceMatrix, _tView> {
    public:
        using tBase = Flow<_tSpaceVector, _tSpaceMatrix, _tView>;
        using typename tBase::tSpaceVector;
        using typename tBase::tSpaceMatrix;
        template<typename... Args> using tView = typename tBase::tView<Args...>;
    public:
        static tSpaceVector getVelocity(const double* pX, const double t) {
            // TODO: BACK TO PREVIOUS DEFINITION
            return 0.5 * tSpaceVector({
                -std::cos(pX[0]) * std::sin(pX[1]), 
                std::sin(pX[0]) * std::cos(pX[1])
            });
            // return tSpaceVector({
                    // std::cos(pX[0]) * std::sin(pX[1]),
                    // -std::sin(pX[0]) * std::cos(pX[1])
            // });
        };
        static tSpaceMatrix getVelocityGradients(const double* pX, const double t) {
            // TODO: BACK TO PREVIOUS DEFINITION
            tSpaceMatrix velocityGradients = tSpaceMatrix::Zero();
            velocityGradients(0,0) = std::sin(pX[0]) * std::sin(pX[1]); velocityGradients(0,1) = -std::cos(pX[0]) * std::cos(pX[1]);
            velocityGradients(1,0) = std::cos(pX[0]) * std::cos(pX[1]); velocityGradients(1,1) = -std::sin(pX[0]) * std::sin(pX[1]);
            return 0.5 * velocityGradients;
            // velocityGradients(0,0) = -std::sin(pX[0]) * std::sin(pX[1]); velocityGradients(0,1) = std::cos(pX[0]) * std::cos(pX[1]);
            // velocityGradients(1,0) = -std::cos(pX[0]) * std::cos(pX[1]); velocityGradients(1,1) = std::sin(pX[0]) * std::sin(pX[1]);
            // return velocityGradients;
        };
        static tSpaceVector getAcceleration(const double* pX, const double t) {
            return tSpaceVector::Zero();
        };
};

}

}

#endif
