#ifndef FL0W_TAYLOR_GREEN_VORTEX_H
#define FL0W_TAYLOR_GREEN_VORTEX_H
#pragma once

// std includes
#include <cmath> // sin, cos
// module includes
#include "fl0w/flow.h"

namespace fl0w {

template<typename TypeVector, typename TypeMatrix, template<typename...> class TypeRef>
class TaylorGreenVortex : public Flow<TypeVector, TypeMatrix, TypeRef> {
    public:
        TaylorGreenVortex();

        TypeVector getVelocity(const TypeRef<const TypeVector>& x, const double& t) const override;
        TypeMatrix getJacobian(const TypeRef<const TypeVector>& x, const double& t) const override;
        TypeVector getAcceleration(const TypeRef<const TypeVector>& x, const double& t) const override;
};

// TaylorGreenVortex class

template<typename TypeVector, typename TypeMatrix, template<typename...> class TypeRef>
TaylorGreenVortex<TypeVector, TypeMatrix, TypeRef>::TaylorGreenVortex() : Flow<TypeVector, TypeMatrix, TypeRef>::Flow() { 

}

template<typename TypeVector, typename TypeMatrix, template<typename...> class TypeRef>
TypeVector TaylorGreenVortex<TypeVector, TypeMatrix, TypeRef>::getVelocity(const TypeRef<const TypeVector>& x, const double& t) const {
    TypeVector u = TypeVector::Zero();
    u(0) = std::cos(x(0)) * std::sin(x(1));
    u(1) = -std::sin(x(0)) * std::cos(x(1));
    return u;
}

template<typename TypeVector, typename TypeMatrix, template<typename...> class TypeRef>
TypeMatrix TaylorGreenVortex<TypeVector, TypeMatrix, TypeRef>::getJacobian(const TypeRef<const TypeVector>& x, const double& t) const {
    TypeMatrix J = TypeMatrix::Zero();
    J(0,0) = -std::sin(x(0)) * std::sin(x(1)); J(0,1) = std::cos(x[0]) * std::cos(x[1]);
    J(1,0) = -std::cos(x[0]) * std::cos(x(1)); J(1,1) = std::sin(x(0)) * std::sin(x(1));
    return J;
}

template<typename TypeVector, typename TypeMatrix, template<typename...> class TypeRef>
TypeVector TaylorGreenVortex<TypeVector, TypeMatrix, TypeRef>::getAcceleration(const TypeRef<const TypeVector>& x, const double& t) const {
    return TypeVector::Zero();
}

}

#endif
