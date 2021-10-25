#ifndef FL0W_FLOW_H
#define FL0W_FLOW_H
#pragma once

// std includes
#include <cstddef> // size_t

namespace fl0w {

template<typename TypeVector, typename TypeMatrix, template<typename...> class TypeRef>
class Flow {
    public:
        Flow();
        
        virtual TypeVector getVelocity(const TypeRef<const TypeVector>& x, const double& t) const = 0;
        virtual TypeMatrix getJacobian(const TypeRef<const TypeVector>& x, const double& t) const = 0;
        virtual TypeVector getAcceleration(const TypeRef<const TypeVector>& x, const double& t) const = 0;
        
        TypeVector getVorticity(const TypeRef<const TypeVector>& x, const double& t) const;
        TypeMatrix getStrain(const TypeRef<const TypeVector>& x, const double& t) const;
};

// Flow class

template<typename TypeVector, typename TypeMatrix, template<typename...> class TypeRef>
Flow<TypeVector, TypeMatrix, TypeRef>::Flow() {

}

template<typename TypeVector, typename TypeMatrix, template<typename...> class TypeRef>
TypeVector Flow<TypeVector, TypeMatrix, TypeRef>::getVorticity(const TypeRef<const TypeVector>& x, const double& t) const {
    // Get jacobian
    TypeMatrix J = getJacobian(x, t);
    // Get vorticity from jacobian
    TypeVector w;
    for(std::size_t i = 0; i < w.size(); i++) {
        w(i) = 0.5 * ( J((i+2) % w.size(), (i+1) % w.size()) - J((i+1) % w.size(), (i+2) % w.size()) );
    }
    return w;
}

template<typename TypeVector, typename TypeMatrix, template<typename...> class TypeRef>
TypeMatrix Flow<TypeVector, TypeMatrix, TypeRef>::getStrain(const TypeRef<const TypeVector>& x, const double& t) const {
    // Get jacobian
    TypeMatrix J = getJacobian(x, t);
    // Get strain from jacobian
    TypeMatrix S;
    for(std::size_t i = 0; i < S.cols(); i++) {
        for(std::size_t j = 0; j < S.cols(); j++) {
            S(i, j) = 0.5 * ( J(i, j) + J(j, i) );
        }
    }
    return S;
}

}

#endif
