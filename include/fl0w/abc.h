#ifndef FL0W_ABC_H
#define FL0W_ABC_H
#pragma once

// std includes
#include <cmath> // sin, cos
// module includes
#include "fl0w/flow.h"

namespace fl0w {

template<typename TypeVector, typename TypeMatrix, template<typename...> class TypeRef>
class ABC : public Flow<TypeVector, TypeMatrix, TypeRef> {
    public:
        ABC();

        TypeVector getVelocity(const TypeRef<const TypeVector>& x, const double& t) const override;
        TypeMatrix getJacobian(const TypeRef<const TypeVector>& x, const double& t) const override;
        TypeVector getAcceleration(const TypeRef<const TypeVector>& x, const double& t) const override;
};

// ABC class

template<typename TypeVector, typename TypeMatrix, template<typename...> class TypeRef>
ABC<TypeVector, TypeMatrix, TypeRef>::ABC() : Flow<TypeVector, TypeMatrix, TypeRef>::Flow() { 

}

template<typename TypeVector, typename TypeMatrix, template<typename...> class TypeRef>
TypeVector ABC<TypeVector, TypeMatrix, TypeRef>::getVelocity(const TypeRef<const TypeVector>& x, const double& t) const {
    TypeVector u;
    u(0) = std::sin(x(2)) + std::cos(x(1));
    u(1) = std::sin(x(0)) + std::cos(x(2));
    u(2) = std::sin(x(1)) + std::cos(x(0));
    return u;
}

template<typename TypeVector, typename TypeMatrix, template<typename...> class TypeRef>
TypeMatrix ABC<TypeVector, TypeMatrix, TypeRef>::getJacobian(const TypeRef<const TypeVector>& x, const double& t) const {
    TypeMatrix J;
    J(0,0) = 0.0;               J(0,1) = -std::sin(x[1]);   J(0,2) = std::cos(x[2]);
    J(1,0) = std::cos(x[0]);    J(1,1) = 0.0;               J(1,2) = -std::sin(x[2]);
    J(2,0) = -std::sin(x[0]);  J(2,1) = std::cos(x[1]);    J(2,2) = 0.0;
    return J;
}

template<typename TypeVector, typename TypeMatrix, template<typename...> class TypeRef>
TypeVector ABC<TypeVector, TypeMatrix, TypeRef>::getAcceleration(const TypeRef<const TypeVector>& x, const double& t) const {
    return TypeVector::Zero();
}

}

#endif
