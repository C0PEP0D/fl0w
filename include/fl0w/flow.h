#ifndef FL0W_FLOW_H
#define FL0W_FLOW_H
#pragma once

// std includes
#include <cstddef> // size_t

namespace fl0w {

template<typename _tSpaceVector, typename _tSpaceMatrix, template<typename...> typename _tView>
class Flow {
	public:
		using tSpaceVector = _tSpaceVector;
		using tSpaceMatrix = _tSpaceMatrix;
		template<typename... Args> using tView = _tView<Args...>;
};

template<typename tFlow, double Time>
class FlowFrozen : public tFlow {
	public:
		using tBase = tFlow;
		using typename tBase::tSpaceVector;
		using typename tBase::tSpaceMatrix;
		template<typename... Args> using tView = typename tBase::tView<Args...>;
	public:
		static tSpaceVector getVelocity(const double* pX, const double t) {
			return tBase::getVelocity(pX, Time);
		}
		static tSpaceMatrix getVelocityGradients(const double* pX, const double t) {
			return tBase::getVelocityGradients(pX, time);
		}
		static tSpaceVector getAcceleration(const double* pX, const double t) {
			return tBase::getAcceleration(pX, time);
		}
};

template<typename tFlow, double TimeOrigin>
class FlowReverse : public tFlow {
	public:
		using tBase = tFlow;
		using typename tBase::tSpaceVector;
		using typename tBase::tSpaceMatrix;
		template<typename... Args> using tView = typename tBase::tView<Args...>;
	public:
		static tSpaceVector getVelocity(const double* pX, const double t) {
			return tBase::getVelocity(pX, TimeOrigin - t);
		}
		static tSpaceMatrix getVelocityGradients(const double* pX, const double t) {
			return tBase::getVelocityGradients(pX, TimeOrigin - t);
		}
		static tSpaceVector getAcceleration(const double* pX, const double t) {
			return tBase::getAcceleration(pX, TimeOrigin - t);
		}
};

}

#endif
