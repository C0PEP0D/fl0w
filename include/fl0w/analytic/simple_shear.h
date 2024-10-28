#ifndef FL0W_SIMPLE_SHEAR_H
#define FL0W_SIMPLE_SHEAR_H
#pragma once

#include "fl0w/flow.h"
#include <limits>

namespace fl0w {

namespace analytic {

template<typename _tSpaceVector, typename _tSpaceMatrix, template<typename...> class _tView, double _ShearRate>
class FlowSimpleShear : public Flow<_tSpaceVector, _tSpaceMatrix, _tView> {
	public:
		using tBase = Flow<_tSpaceVector, _tSpaceMatrix, _tView>;
		using typename tBase::tSpaceVector;
		using typename tBase::tSpaceMatrix;
		template<typename... Args> using tView = typename tBase::tView<Args...>;
		constexpr static double ShearRate = _ShearRate;
    public:
        static tSpaceVector getVelocity(const double* pX, const double t) {
        	return {
        		ShearRate * pX[1],
        		0.0
        	};
        };
        static tSpaceMatrix getVelocityGradients(const double* pX, const double t) {
        	tSpaceMatrix velocityGradients = tSpaceMatrix::Zero();
       	    velocityGradients(0,1) = ShearRate;
       	    return velocityGradients;
        };
        static tSpaceVector getAcceleration(const double* pX, const double t) {
        	return tSpaceVector::Zero();
        };
};

}

}

#endif
