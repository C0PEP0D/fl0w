#ifndef FL0W_FLOW_H
#define FL0W_FLOW_H
#pragma once

// std includes
#include <cstddef> // size_t

namespace fl0w {

template<typename _TypeVector, typename _TypeMatrix, template<typename...> class _TypeRef>
class Flow {
	public:
		using TypeVector = _TypeVector;
		using TypeMatrix = _TypeMatrix;
		template<typename... Args>
		using TypeRef = _TypeRef<Args...>;
	public:
		Flow() {
			
		}
	public:
		virtual TypeVector getVelocity(const TypeRef<const TypeVector>& x, const double& t) const = 0;
		virtual TypeMatrix getVelocityGradients(const TypeRef<const TypeVector>& x, const double& t) const = 0;
		virtual TypeVector getAcceleration(const TypeRef<const TypeVector>& x, const double& t) const = 0;
};

template<typename TypeFlow>
class FlowFrozen : public TypeFlow {
	public:
		using Type = TypeFlow;
		using typename Type::TypeVector;
		using typename Type::TypeMatrix;
	public:
		FlowFrozen(const double& p_time) : TypeFlow(), time(p_time) {
			
		}
	public:
		TypeVector getVelocity(const typename Type::TypeRef<const TypeVector>& x, const double& t) const override {
			return Type::getVelocity(x, time);
		}
		TypeMatrix getVelocityGradients(const typename Type::TypeRef<const TypeVector>& x, const double& t) const override {
			return Type::getVelocityGradients(x, time);
		}
		TypeVector getAcceleration(const typename Type::TypeRef<const TypeVector>& x, const double& t) const override {
			return Type::getAcceleration(x, time);
		}
	public:
		double time;
};

template<typename TypeFlow>
class FlowReverse : public TypeFlow {
	public:
		using Type = TypeFlow;
		using typename Type::TypeVector;
		using typename Type::TypeMatrix;
	public:
		FlowReverse(const double& p_time) : TypeFlow(), time(p_time) {
			
		}
	public:
		TypeVector getVelocity(const typename Type::TypeRef<const TypeVector>& x, const double& t) const override {
			return -Type::getVelocity(x, time - t);
		}
		TypeMatrix getVelocityGradients(const typename Type::TypeRef<const TypeVector>& x, const double& t) const override {
			return -Type::getVelocityGradients(x, time - t);
		}
		TypeVector getAcceleration(const typename Type::TypeRef<const TypeVector>& x, const double& t) const override {
			return -Type::getAcceleration(x, time - t);
		}
	public:
		double time;
};

}

#endif
