#ifndef FLOW_H
#define FLOW_H
#pragma once

#include "fl0w/flow.h"
#include <cmath>
#include <limits>
#include <memory>

namespace fl0w {

template<typename TypeVector, typename TypeMatrix, template<typename...> class TypeRef>
class Turbulence : public Flow<TypeVector, TypeMatrix, TypeRef> {
    public:
        Turbulence();
        virtual void createFromDissipationRate(const double& nu, const double& eps, const double& l);
        virtual void createFromCharacteristicLengthsAndTime(const double& eta, const double& tauEta, const double& l);
        
        double nu; // Flow kinematic viscosity
        std::shared_ptr<double> sEps; // Turbulence dissipation rate
        double eta; // Characteristic length of small eddys
        double l; // Characteristic length of big eddys
        double tauEta; // Characteristic time of small eddys
        double uEta; // Characteristic velocity of small eddys
};
// Turbulence class

template<typename TypeVector, typename TypeMatrix, template<typename...> class TypeRef>
Turbulence<TypeVector, TypeMatrix, TypeRef>::Turbulence() : Flow<TypeVector, TypeMatrix, TypeRef>::Flow(),
    nu(std::numeric_limits<double>::quiet_NaN()), 
    eta(std::numeric_limits<double>::quiet_NaN()),
    l(std::numeric_limits<double>::quiet_NaN()),
    tauEta(std::numeric_limits<double>::quiet_NaN()),
    uEta(std::numeric_limits<double>::quiet_NaN()) 
{
    sEps = std::make_shared<double>(std::numeric_limits<double>::quiet_NaN());
    createFromDissipationRate(1e-6, 1e-5, 1e-2);
}

template<typename TypeVector, typename TypeMatrix, template<typename...> class TypeRef>
void Turbulence<TypeVector, TypeMatrix, TypeRef>::createFromDissipationRate(const double& p_nu, const double& p_eps, const double& p_l) {
    nu = p_nu;
    *sEps = p_eps;
    l = p_l;
    
    eta = std::pow(std::pow(nu, 3)/(*sEps), 0.25); // Small eddys length scale
    tauEta = std::pow(nu/(*sEps), 0.5);
    uEta = eta/tauEta;
}

template<typename TypeVector, typename TypeMatrix, template<typename...> class TypeRef>
void Turbulence<TypeVector, TypeMatrix, TypeRef>::createFromCharacteristicLengthsAndTime(const double& p_eta, const double& p_tauEta, const double& p_l) {
    eta = p_eta;
    tauEta = p_tauEta;
    l = p_l;

    uEta = eta/tauEta;
    nu = std::pow(uEta, 2)/tauEta;
    *sEps = nu/std::pow(tauEta, 2);
}

}

#endif
