#ifndef FLOW_H
#define FLOW_H
#pragma once

#include "fl0w/flow.h"
#include <cmath>
#include <limits>

namespace fl0w {

namespace analytic {

template<typename TypeVector, typename TypeMatrix, template<typename...> class TypeRef>
class BatchelorVortex : public Flow<TypeVector, TypeMatrix, TypeRef> {
    public:
        BatchelorVortex();

        void create(const double& q, const double& r1, const double& r2);
        TypeVector getVelocity(const TypeRef<const TypeVector>& x, const double& t) const override;
        TypeMatrix getVelocityGradients(const TypeRef<const TypeVector>& x, const double& t) const override;
        TypeVector getAcceleration(const TypeRef<const TypeVector>& x, const double& t) const override;

        // Attributes :
        
        double q;
        double r1;
        double r2;
};
// BatchelorVortex class

template<typename TypeVector, typename TypeMatrix, template<typename...> class TypeRef>
BatchelorVortex<TypeVector, TypeMatrix, TypeRef>::BatchelorVortex() : Flow<TypeVector, TypeMatrix, TypeRef>::Flow(), 
    q(std::numeric_limits<double>::quiet_NaN()),
    r1(std::numeric_limits<double>::quiet_NaN()),
    r2(std::numeric_limits<double>::quiet_NaN()) {

}

template<typename TypeVector, typename TypeMatrix, template<typename...> class TypeRef>
void BatchelorVortex<TypeVector, TypeMatrix, TypeRef>::create(const double& p_q, const double& p_r1, const double& p_r2) {
    q = p_q;
    r1 = p_r1;
    r2 = p_r2;
}

template<typename TypeVector, typename TypeMatrix, template<typename...> class TypeRef>
TypeVector BatchelorVortex<TypeVector, TypeMatrix, TypeRef>::getVelocity(const TypeRef<const TypeVector>& x, const double& t) const {
    TypeVector vr; vr = x; vr(2) = 0.0;
    const double r = vr.norm();
    TypeVector er;
    TypeVector et;
    if (r > 0.0) {
        // Compute cylindric base
        er = vr/r;
        et = er;
        // Swap coordinates to get the normal vector
        et(0) = -er(1); et(1) = er(0);
    }
    TypeVector ez; ez.fill(0.0); ez(2) = 1.0;
    
    if (r > 0.0) {
        return (1.0 - std::exp(-std::pow((r/r1), 2))) / (r / r1) * et + q * exp(-std::pow(r/r2, 2)) * ez;
    } else {
        return q * ez;
    }
}

template<typename TypeVector, typename TypeMatrix, template<typename...> class TypeRef>
TypeVector BatchelorVortex<TypeVector, TypeMatrix, TypeRef>::getAcceleration(const TypeRef<const TypeVector>& x, const double& t) const {
    TypeVector a; a.fill(0.0);
    return a;
}

template<typename TypeVector, typename TypeMatrix, template<typename...> class TypeRef>
TypeMatrix BatchelorVortex<TypeVector, TypeMatrix, TypeRef>::getVelocityGradients(const TypeRef<const TypeVector>& x, const double& t) const {
    TypeVector vr; vr = x; vr(2) = 0.0;
    const double r = vr.norm();
    TypeVector er;
    TypeVector et;
    if (r > 0.0) {
        // Compute cylindric base
        er = vr/r;
        et = er;
        // Swap coordinates to get the normal vector
        et(0) = -er(1); et(1) = er(0);
    }
    TypeVector ez; ez.fill(0.0); ez(2) = 1.0;

    TypeVector dr = ((2.0 * std::pow(r, 2) - std::pow(r1, 2)) * std::exp(-std::pow(r/r1, 2)) - std::pow(r1, 2)) / (r1 * std::pow(r, 2)) * et + -q * r/std::pow(r2, 2) * std::exp(-std::pow(r/r2, 2)) * ez;
    TypeMatrix J;
    if (r > 0.0) {
        for(std::size_t i = 0; i < x.size(); i++) {
            for(std::size_t j = 0; j < x.size(); j++) {
                J(i, j) = dr(i) * x(j)/r;
            }
        }
    } else {
        J.fill(0.0);
    }
    return J;
}

}

}

#endif
