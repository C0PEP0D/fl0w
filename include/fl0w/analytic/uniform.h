#ifndef FL0W_UNIFORM_H
#define FL0W_UNIFORM_H
#pragma once

#include "fl0w/flow.h"
#include <array>

namespace fl0w {

namespace analytic {

template<typename _tSpaceVector, typename _tSpaceMatrix, template<typename...> class _tView, std::array<double, _tSpaceVector::SizeAtCompileTime> Velocity>
class FlowUniform : public Flow<_tSpaceVector, _tSpaceMatrix, _tView> {
    public:
        using tBase = Flow<_tSpaceVector, _tSpaceMatrix, _tView>;
        using typename tBase::tSpaceVector;
        using typename tBase::tSpaceMatrix;
        template<typename... Args> using tView = typename tBase::tView<Args...>;
    public:
        static tSpaceVector getVelocity(const double* pX, const double t) {
            return tView<const tSpaceVector>(Velocity.data());
        };
        static tSpaceMatrix getVelocityGradients(const double* pX, const double t) {
            return tSpaceMatrix::Zero();
        };
        static tSpaceVector getAcceleration(const double* pX, const double t) {
            return tSpaceVector::Zero();
        };
};

}

}

#endif
