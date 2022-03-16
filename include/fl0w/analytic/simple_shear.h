#ifndef FL0W_SIMPLE_SHEAR_H
#define FL0W_SIMPLE_SHEAR_H
#pragma once

#include "fl0w/flow.h"
#include <limits>

namespace fl0w {

namespace analytic {

template<typename TypeVector, typename TypeMatrix, template<typename...> class TypeRef>
class SimpleShear : public Flow<TypeVector, TypeMatrix, TypeRef> {
    public:
        SimpleShear();

        void create(const double& gamma);
        TypeVector getVelocity(const TypeRef<const TypeVector>& x, const double& t) const override;
        TypeMatrix getVelocityGradients(const TypeRef<const TypeVector>& x, const double& t) const override;
        TypeVector getAcceleration(const TypeRef<const TypeVector>& x, const double& t) const override;

        // Attributes :
        double gamma;
};

// SimpleShear class

template<typename TypeVector, typename TypeMatrix, template<typename...> class TypeRef>
SimpleShear<TypeVector, TypeMatrix, TypeRef>::SimpleShear() : Flow<TypeVector, TypeMatrix, TypeRef>::Flow(), gamma(std::numeric_limits<double>::quiet_NaN()) {

}

template<typename TypeVector, typename TypeMatrix, template<typename...> class TypeRef>
void SimpleShear<TypeVector, TypeMatrix, TypeRef>::create(const double& p_gamma) {
    gamma = p_gamma;
}

template<typename TypeVector, typename TypeMatrix, template<typename...> class TypeRef>
TypeVector SimpleShear<TypeVector, TypeMatrix, TypeRef>::getVelocity(const TypeRef<const TypeVector>& x, const double& t) const {
    TypeVector u; u.fill(0.0);
    u(0) = gamma * x(1);
    return u;
}

template<typename TypeVector, typename TypeMatrix, template<typename...> class TypeRef>
TypeMatrix SimpleShear<TypeVector, TypeMatrix, TypeRef>::getVelocityGradients(const TypeRef<const TypeVector>& x, const double& t) const {
    TypeMatrix J; J.fill(0.0);
    J(0,1) = gamma;
    return J;
}

template<typename TypeVector, typename TypeMatrix, template<typename...> class TypeRef>
TypeVector SimpleShear<TypeVector, TypeMatrix, TypeRef>::getAcceleration(const TypeRef<const TypeVector>& x, const double& t) const {
    TypeVector a; a.fill(0.0);
    return a;
}

}

}

#endif
