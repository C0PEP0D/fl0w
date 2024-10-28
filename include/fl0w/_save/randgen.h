#ifndef FL0W_RANDGEN_H
#define FL0W_RANDGEN_H
#pragma once

#include <random>
#include <chrono>

namespace randgen {

template<typename TypeScalar>
class Generator {
    public:
        Generator();
        TypeScalar uniform();
        
        template<typename TypeVector>
        TypeVector unitVector();
        
        void setSeed(std::size_t seed);
        void randomizeSeed();
    public:
        std::random_device rd;
        std::mt19937_64 generator; //Mersenne Twister by Matsumoto and Nishimura, 2000
        std::uniform_real_distribution<TypeScalar> distrib;
};

template<typename TypeScalar>
Generator<TypeScalar>::Generator() : rd(), generator(rd()), distrib(0.0, std::nextafter(1.0, std::numeric_limits<TypeScalar>::max())) {
    randomizeSeed();
}

template<typename TypeScalar>
void Generator<TypeScalar>::setSeed(std::size_t val) {
    generator.seed(val);
}

template<typename TypeScalar>
void Generator<TypeScalar>::randomizeSeed() {
    setSeed(std::chrono::system_clock::now().time_since_epoch().count());
}

template<typename TypeScalar>
TypeScalar Generator<TypeScalar>::uniform() {
    return distrib(generator);
}

template<typename TypeScalar> 
template<typename TypeVector>
TypeVector Generator<TypeScalar>::unitVector() {
    TypeVector vec;
    for(std::size_t i = 0; i < vec.size(); i++) {
        vec(i) = 2.0 * uniform() - 1.0;
    }
    vec.normalize();
    return vec;
}

}

#endif
