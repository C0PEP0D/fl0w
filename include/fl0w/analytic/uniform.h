#ifndef FL0W_ABC_H
#define FL0W_ABC_H
#pragma once

#include "fl0w/flow.h"

namespace fl0w {

namespace analytic {

template<typename TypeVector, typename TypeMatrix, template<typename...> class TypeRef>
class Uniform : public Flow<TypeVector, TypeMatrix, TypeRef> {
    public:
        Uniform();

        TypeVector getVelocity(const TypeRef<const TypeVector>& x, const double& t) const override;
        TypeMatrix getVelocityGradients(const TypeRef<const TypeVector>& x, const double& t) const override;
        TypeVector getAcceleration(const TypeRef<const TypeVector>& x, const double& t) const override;
    public:
        TypeVector u;
};

// ABC class

template<typename TypeVector, typename TypeMatrix, template<typename...> class TypeRef>
Uniform<TypeVector, TypeMatrix, TypeRef>::Uniform() : Flow<TypeVector, TypeMatrix, TypeRef>::Flow() : u(TypeVector::Zero()) { 
}

template<typename TypeVector, typename TypeMatrix, template<typename...> class TypeRef>
TypeVector Uniform<TypeVector, TypeMatrix, TypeRef>::getVelocity(const TypeRef<const TypeVector>& x, const double& t) const {
    return u;
}

template<typename TypeVector, typename TypeMatrix, template<typename...> class TypeRef>
TypeMatrix Uniform<TypeVector, TypeMatrix, TypeRef>::getVelocityGradients(const TypeRef<const TypeVector>& x, const double& t) const {
    return TypeMatrix::Zero();
}

template<typename TypeVector, typename TypeMatrix, template<typename...> class TypeRef>
TypeVector Uniform<TypeVector, TypeMatrix, TypeRef>::getAcceleration(const TypeRef<const TypeVector>& x, const double& t) const {
    return TypeVector::Zero();
}

}

}

#endif
