#ifndef FL0W_KINEMATIC_H
#define FL0W_KINEMATIC_H
#pragma once

// Include Standards
#include <memory>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <iostream>
// Include from fl0w
#include "fl0w/turbulence.h"
#include "fl0w/randgen.h"

namespace fl0w {

namespace kinematic {
        
class Spectrum {
    public:
        inline Spectrum();
        inline virtual const double operator()(const double& k) const = 0;
};

class SpectrumKolmogorov : public Spectrum {
    public:
        inline SpectrumKolmogorov();
        inline const double operator()(const double& k) const override;

        std::shared_ptr<double> sK0;
        std::shared_ptr<double> sEps;
};

template<typename TypeGenerated, template<typename...> class TypeContainer>
class Generator {
    public:
        Generator();
        virtual const TypeContainer<TypeGenerated>& gen(TypeContainer<TypeGenerated>& toGen) const = 0;
    public:
        std::shared_ptr<std::size_t> sNk;
        std::shared_ptr<randgen::Generator<double>> sRandomGenerator;
};

template<typename TypeGenerated, template<typename...> class TypeContainer>
class GeneratorWaveVector : public Generator<TypeGenerated, TypeContainer> {
    public:
        GeneratorWaveVector();
        const TypeContainer<TypeGenerated>& gen(TypeContainer<TypeGenerated>& k) const override = 0;
    public:
        double inf;
        double sup;
        using Generator<TypeGenerated, TypeContainer>::sNk;
        using Generator<TypeGenerated, TypeContainer>::sRandomGenerator;
};

template<typename TypeGenerated, template<typename...> class TypeContainer>
class GeneratorWaveVectorLinear : public GeneratorWaveVector<TypeGenerated, TypeContainer> {
    public:
        using GeneratorWaveVector<TypeGenerated, TypeContainer>::GeneratorWaveVector;
        const TypeContainer<TypeGenerated>& gen(TypeContainer<TypeGenerated>& k) const override;
    
    public:
        using GeneratorWaveVector<TypeGenerated, TypeContainer>::sNk;
        using GeneratorWaveVector<TypeGenerated, TypeContainer>::inf;
        using GeneratorWaveVector<TypeGenerated, TypeContainer>::sup;
        using GeneratorWaveVector<TypeGenerated, TypeContainer>::sRandomGenerator;
};

template<typename TypeGenerated, template<typename...> class TypeContainer>
class GeneratorWaveVectorGeometric : public GeneratorWaveVector<TypeGenerated, TypeContainer> {
    public:
        using GeneratorWaveVector<TypeGenerated, TypeContainer>::GeneratorWaveVector;
        const TypeContainer<TypeGenerated>& gen(TypeContainer<TypeGenerated>& k) const override;
    
    public:
        using GeneratorWaveVector<TypeGenerated, TypeContainer>::sNk;
        using GeneratorWaveVector<TypeGenerated, TypeContainer>::inf;
        using GeneratorWaveVector<TypeGenerated, TypeContainer>::sup;
        using GeneratorWaveVector<TypeGenerated, TypeContainer>::sRandomGenerator;
};

template<typename TypeGenerated, template<typename...> class TypeContainer>
class GeneratorAmplitude : public Generator<TypeGenerated, TypeContainer> {
    public:
        GeneratorAmplitude();
        const TypeContainer<TypeGenerated>& gen(TypeContainer<TypeGenerated>& a) const override;
    public:
        using Generator<TypeGenerated, TypeContainer>::sNk;
        using Generator<TypeGenerated, TypeContainer>::sRandomGenerator;
    public:
        std::shared_ptr<TypeContainer<TypeGenerated>> sK;
        std::shared_ptr<Spectrum> sE;
};

template<typename TypeVector, template<typename...> class TypeContainer>
class GeneratorOmega : public Generator<double, TypeContainer> {
    public:
        GeneratorOmega();
        const TypeContainer<double>& gen(TypeContainer<double>& omega) const override = 0;
    public:
        std::shared_ptr<TypeContainer<TypeVector>> sK;
        std::shared_ptr<Spectrum> sE;
        using Generator<double, TypeContainer>::sNk;
        using Generator<double, TypeContainer>::sRandomGenerator;
};

template<typename TypeVector, template<typename...> class TypeContainer>
class GeneratorOmegaLewis2000 : public GeneratorOmega<TypeVector, TypeContainer> {
    public:
        GeneratorOmegaLewis2000();
        const TypeContainer<double>& gen(TypeContainer<double>& omega) const override;
    public:
        using GeneratorOmega<TypeVector, TypeContainer>::sK;
        using GeneratorOmega<TypeVector, TypeContainer>::sE;
        using GeneratorOmega<TypeVector, TypeContainer>::sNk;
        using GeneratorOmega<TypeVector, TypeContainer>::sRandomGenerator;
    public:
        double lambda;
};

template<typename TypeVector, typename TypeMatrix, template<typename...> class TypeRef, template<typename...> class TypeContainer>
class Kinematic : public Turbulence<TypeVector, TypeMatrix, TypeRef> {
    public:
        // Methods :
        Kinematic();
        void createFromDissipationRate(const double& nu, const double& eps, const double& l) override;
        void createFromCharacteristicLengthsAndTime(const double& eta, const double& tauEta, const double& l) override;
        
        void init();
        
        TypeVector getVelocity(const TypeRef<const TypeVector>& x, const double& t) const override;
        TypeVector getAcceleration(const TypeRef<const TypeVector>& x, const double& t) const override;
        TypeMatrix getVelocityGradients(const TypeRef<const TypeVector>& x, const double& t) const override;

        // Attributes :
        std::shared_ptr<randgen::Generator<double>> sRandomGenerator;

        std::shared_ptr<double> sK0;
        std::shared_ptr<std::size_t> sNk; // Number of Modes
        TypeContainer<double> periodicBox; // Periodic box
        std::shared_ptr<SpectrumKolmogorov> sE;
        std::shared_ptr<GeneratorWaveVector<TypeVector, TypeContainer>> sKGen;
        std::shared_ptr<TypeContainer<TypeVector>> sK;
        std::shared_ptr<GeneratorAmplitude<TypeVector, TypeContainer>> sAGen;
        TypeContainer<TypeVector> a;
        TypeContainer<TypeVector> b;
        std::shared_ptr<GeneratorOmegaLewis2000<TypeVector, TypeContainer>> sOGen;
        TypeContainer<double> omega;
        
        struct Data {enum size_t{k, a, b, omega};};
        TypeContainer<std::tuple<TypeVector, TypeVector, TypeVector, double>> data; // k, a, b, omega
        
        // Using Attributes
        using Turbulence<TypeVector, TypeMatrix, TypeRef>::nu; // Flow kinematic viscosity
        using Turbulence<TypeVector, TypeMatrix, TypeRef>::sEps; // Turbulence dissipation rate
        using Turbulence<TypeVector, TypeMatrix, TypeRef>::eta; // Characteristic length of small eddys
        using Turbulence<TypeVector, TypeMatrix, TypeRef>::l; // Characteristic length of big eddys
        using Turbulence<TypeVector, TypeMatrix, TypeRef>::tauEta; // Characteristic time of small eddys
        using Turbulence<TypeVector, TypeMatrix, TypeRef>::uEta; // Characteristic velocity of small eddys
};

// Kinematic Class

template<typename TypeVector, typename TypeMatrix, template<typename...> class TypeRef, template<typename...> class TypeContainer>
Kinematic<TypeVector, TypeMatrix, TypeRef, TypeContainer>::Kinematic() : Turbulence<TypeVector, TypeMatrix, TypeRef>::Turbulence() {
    // Basic Constants
    sK0 = std::make_shared<double>(1.6);
    sNk = std::make_shared<std::size_t>(100);
    sRandomGenerator = std::make_shared<randgen::Generator<double>>();
    // Wave Vector Generator
    sKGen = std::make_shared<GeneratorWaveVectorLinear<TypeVector, TypeContainer>>();
    sKGen->sRandomGenerator = sRandomGenerator;
    sKGen->sNk = sNk;
    // Wave Vectors
    sK = std::make_shared<TypeContainer<TypeVector>>(*sNk);
    // Create other turbulent constants
    createFromDissipationRate(1e-6, 1e-5, 1e-2);
    // Spectrum
    sE = std::make_shared<SpectrumKolmogorov>();
    sE->sK0 = sK0;
    sE->sEps = sEps;
    // Amplitude generator
    sAGen = std::make_shared<GeneratorAmplitude<TypeVector, TypeContainer>>();
    sAGen->sRandomGenerator = sRandomGenerator;
    sAGen->sNk = sNk;
    sAGen->sK = sK;
    sAGen->sE = sE;
    // Omega generator
    sOGen = std::make_shared<GeneratorOmegaLewis2000<TypeVector, TypeContainer>>();
    sOGen->sRandomGenerator = sRandomGenerator;
    sOGen->sNk = sNk;
    sOGen->sK = sK;
    sOGen->sE = sE;
    sOGen->lambda = 0.4;
}

template<typename TypeVector, typename TypeMatrix, template<typename...> class TypeRef, template<typename...> class TypeContainer>
void Kinematic<TypeVector, TypeMatrix, TypeRef, TypeContainer>::createFromDissipationRate(const double& p_nu, const double& p_eps, const double& p_l) {
    Turbulence<TypeVector, TypeMatrix, TypeRef>::createFromDissipationRate(p_nu, p_eps, p_l);
    // And set inf and sup
    sKGen->inf = 2.0 * M_PI/l;
    sKGen->sup = 2.0 * M_PI/eta;
}

template<typename TypeVector, typename TypeMatrix, template<typename...> class TypeRef, template<typename...> class TypeContainer>
void Kinematic<TypeVector, TypeMatrix, TypeRef, TypeContainer>::createFromCharacteristicLengthsAndTime(const double& p_eta, const double& p_tauEta, const double& p_l) {
    Turbulence<TypeVector, TypeMatrix, TypeRef>::createFromCharacteristicLengthsAndTime(p_eta, p_tauEta, p_l);
    // And set inf and sup
    sKGen->inf = 2.0 * M_PI/l;
    sKGen->sup = 2.0 * M_PI/eta;
}

template<typename TypeVector, typename TypeMatrix, template<typename...> class TypeRef, template<typename...> class TypeContainer>
void Kinematic<TypeVector, TypeMatrix, TypeRef, TypeContainer>::init() {
    // Generate wave vectors
    *sK = sKGen->gen(*sK);
    // Ensure periodicity if needed
    if(not periodicBox.empty()) {
        TypeVector kBox;
        std::transform(std::begin(periodicBox), std::end(periodicBox), std::begin(kBox), [](const double& pxBox) {
            return 2 * M_PI / pxBox;
        });
        for(size_t n = 0; n < sK->size(); n++) {
            std::transform(std::begin((*sK)[n]), std::end((*sK)[n]), std::begin(kBox), std::begin((*sK)[n]), [](const double& kx, const double& kxBox) {
                return kxBox * std::round(kx/kxBox);
            });
        }
        // Resort just in case some k are misplaced
        std::sort(std::begin(*sK), std::end(*sK), [](const TypeVector& a, const TypeVector& b){
            return a.norm() < b.norm();
        });
    }
    // Generate amplitudes
    a = sAGen->gen(a);
    b = sAGen->gen(b);
    // Generate time pulsations
    omega = sOGen->gen(omega);
    // Creating data as tuple
    data = TypeContainer<std::tuple<TypeVector, TypeVector, TypeVector, double>> (*sNk);
    for(unsigned int n = 0; n < *sNk; n++) {
        data[n] = std::make_tuple((*sK)[n], a[n], b[n], omega[n]);
    }
}

template<typename TypeVector, typename TypeMatrix, template<typename...> class TypeRef, template<typename...> class TypeContainer>
TypeVector Kinematic<TypeVector, TypeMatrix, TypeRef, TypeContainer>::getVelocity(const TypeRef<const TypeVector>& x, const double& t) const {
    // Zero velocity creation
    TypeVector u; u.fill(0.0);
    // TODO : should be changed to reduce, c++17
    for(const auto& dn : data) {
        TypeVector kn; TypeVector an; TypeVector bn; double omegan;
        std::tie(kn, an, bn, omegan) = dn;
        const double y = kn.dot(x) + omegan*t;
        u += an * std::cos(y) + bn * std::sin(y);
    }
    return u;
}

template<typename TypeVector, typename TypeMatrix, template<typename...> class TypeRef, template<typename...> class TypeContainer>
TypeVector Kinematic<TypeVector, TypeMatrix, TypeRef, TypeContainer>::getAcceleration(const TypeRef<const TypeVector>& x, const double& t) const {
    // zero velocity creation
    TypeVector a; a.fill(0.0);
    // TODO : should be changed to reduce, c++17
    for(const auto& dn : data) {
        TypeVector kn; TypeVector an; TypeVector bn; double omegan;
        std::tie(kn, an, bn, omegan) = dn;
        const double y = kn.dot(x) + omegan*t;
        a += -omegan * std::sin(y) * an + omegan * std::cos(y) * bn;
    }
    return a;
}

template<typename TypeVector, typename TypeMatrix, template<typename...> class TypeRef, template<typename...> class TypeContainer>
TypeMatrix Kinematic<TypeVector, TypeMatrix, TypeRef, TypeContainer>::getVelocityGradients(const TypeRef<const TypeVector>& x, const double& t) const {
    TypeMatrix J;
    for(size_t i = 0; i < x.size(); i++) {
        // zero u creation
        TypeVector u; u.fill(0.0);
        // u accumulation
        for(const auto& dn : data) {
            TypeVector kn; TypeVector an; TypeVector bn; double omegan;
            std::tie(kn, an, bn, omegan) = dn;
            const double y = kn.dot(x) + omegan*t;
            u += -kn(i) * std::sin(y) * an + kn(i) * std::cos(y) * bn;
        }
        J.col(i) = u;
    }
    return J;
}

// Spectrum class

Spectrum::Spectrum() {
    
}

// SpectrumKolmogorov class

SpectrumKolmogorov::SpectrumKolmogorov() : Spectrum() {
    
}

const double SpectrumKolmogorov::operator()(const double& k) const {
    // Physical Parameters
    return (*sK0) * std::pow((*sEps), 2.0/3.0) * std::pow(k, -5.0/3.0);
}

// Generator class

template<typename TypeGenerated, template<typename...> class TypeContainer>
Generator<TypeGenerated, TypeContainer>::Generator() {

}

// GeneratorWaveVector class

template<typename TypeGenerated, template<typename...> class TypeContainer>
GeneratorWaveVector<TypeGenerated, TypeContainer>::GeneratorWaveVector() : Generator<TypeGenerated, TypeContainer>() {
    
}

// GeneratorWaveVectorLinear class

template<typename TypeGenerated, template<typename...> class TypeContainer>
const TypeContainer<TypeGenerated>& GeneratorWaveVectorLinear<TypeGenerated, TypeContainer>::gen(TypeContainer<TypeGenerated>& k) const {
    // Get Parameters from Pointers
    const std::size_t& nk = *sNk;

    //Resize vector
    k.resize(nk);
    
    // Generate wave vectors
    const double delta = (sup - inf)/(nk-1);
    double accum = inf;
    for(int n = 0; n < nk; n++) {
        k[n] = sRandomGenerator->template unitVector<TypeGenerated>() * accum;
        accum += delta;
    }
    return k;
}

// GeneratorWaveVectorGeometric class

template<typename TypeGenerated, template<typename...> class TypeContainer>
const TypeContainer<TypeGenerated>& GeneratorWaveVectorGeometric<TypeGenerated, TypeContainer>::gen(TypeContainer<TypeGenerated>& k) const {
    // Get Parameters from Pointers
    const std::size_t& nk = *sNk;
    
    //Resize vector
    k.resize(nk);
    
    // Generate wave vectors
    const double alpha = std::pow(sup/inf, 1/(nk-1));
    double accum = inf;
    for(int n = 0; n < nk; n++) {
        k[n] = sRandomGenerator->template unitVector<TypeGenerated>() * accum;
        accum *= alpha;
    }
    return k;
}

// GeneratorAmplitude class

template<typename TypeGenerated, template<typename...> class TypeContainer>
GeneratorAmplitude<TypeGenerated, TypeContainer>::GeneratorAmplitude() : Generator<TypeGenerated, TypeContainer>() {

}

template<typename TypeGenerated, template<typename...> class TypeContainer>
const TypeContainer<TypeGenerated>& GeneratorAmplitude<TypeGenerated, TypeContainer>::gen(TypeContainer<TypeGenerated>& a) const {
    // Get Parameters from Pointers
    const std::size_t& nk = *sNk;
    const TypeContainer<TypeGenerated>& k = *sK;
    const Spectrum& e = *sE;
    
    //Resize vector
    a.resize(nk);
    
    // Get width of each point of the spectrum
    TypeContainer<double> dk(nk);
    dk[0] = 0.5 * (k[1].norm() - k[0].norm());
    for(std::size_t n = 1; n < nk-1; n++) {
        dk[n] = 0.5 * (k[n+1].norm() - k[n-1].norm());
    }
    dk[nk-1] = 0.5 * (k[nk-1].norm() - k[nk-2].norm());
    
    // Generate wave amplitudes
    for(std::size_t n = 0; n < nk; n++) {
        TypeGenerated randVec = sRandomGenerator->template unitVector<TypeGenerated>();
        // Computing the direction of k[n]
        const double normKn = k[n].norm();
        TypeGenerated dirKn = k[n]/normKn;
        randVec -= randVec.dot(dirKn) * dirKn;
        randVec.normalize();
        // Setting the amplitude a_n
        a[n] = randVec * std::sqrt(2.0 * e(normKn) * dk[n]);
    }
    return a;
}

// GeneratorOmega class

template<typename TypeVector, template<typename...> class TypeContainer>
GeneratorOmega<TypeVector, TypeContainer>::GeneratorOmega() : Generator<double, TypeContainer>() {

}

// GeneratorOmegaLewis2000 class

template<typename TypeVector, template<typename...> class TypeContainer>
GeneratorOmegaLewis2000<TypeVector, TypeContainer>::GeneratorOmegaLewis2000() : GeneratorOmega<TypeVector, TypeContainer>() {

}

template<typename TypeVector, template<typename...> class TypeContainer>
const TypeContainer<double>& GeneratorOmegaLewis2000<TypeVector, TypeContainer>::gen(TypeContainer<double>& omega) const {
    // Get Parameters from Pointers
    const std::size_t& nk = *sNk;
    const TypeContainer<TypeVector>& k = *sK;
    const Spectrum& e = *sE;
    
    //Resize vector
    omega.resize(nk);
    // Generate wave pulsation
    for(std::size_t n = 0; n < nk; n++) {
        // Setting the amplitude T
        omega[n] = lambda * std::sqrt(std::pow(k[n].norm(), 3.0) * e(k[n].norm()));
    }
    return omega;
}

}

}

#endif
