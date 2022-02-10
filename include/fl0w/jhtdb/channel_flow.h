#ifndef FL0W_JHTDB_CHANNEL_FLOW_H
#define FL0W_JHTDB_CHANNEL_FLOW_H
#pragma once

#include "fl0w/jhtdb/jhtdb.h"

namespace fl0w {

namespace jhtdb {

template<typename TypeVector, typename TypeMatrix, template<typename...> class TypeRef>
class ChannelFlow : public JHTDB<TypeVector, TypeMatrix, TypeRef> {
	public:
		ChannelFlow() : JHTDB<TypeVector, TypeMatrix, TypeRef>::JHTDB("channel") {
			// physical parameter
			h = 1.0;
			// numerical parameters
			lx = 8.0 * M_PI * h;
			ly = 2.0 * h;
			lz = 3.0 * M_PI * h;
			nx = 2048;
			ny = 512;
			nz = 1536;
			dt = 0.0065;
			tMax = 25.9935;
			// physical parameters
			viscosity = 5e-5;
			meanPressureGradient = 0.0025;
			// statistics
			velocityBulk = 0.99994;
			velocityCenterline = 1.1312;
			velocityFriction = 4.9968e-2;
			viscousLengthScale = 1.0006e-3;
		}
	public:
		virtual void queryPointFromX(float* point, const TypeRef<const TypeVector>& x) const override {
			// x
			if(x[0] < 0) {
				point[0] = lx - std::fmod(-x[0], lx);
			} else {
				point[0] = std::fmod(x[0], lx);
			}
			// y
			if(x[1] + h < 0) {
				point[1] = ly - std::fmod(-(x[1] + h), ly) - h;
			} else {
				point[1] = std::fmod((x[1] + h), ly) - h;
			}
			// z
			if(x[2] < 0) {
				point[2] = lz - std::fmod(-x[2], lz);
			} else {
				point[2] = std::fmod(x[2], lz);
			}
		}
	public:
		// numerical parameters
		using JHTDB<TypeVector, TypeMatrix, TypeRef>::lx;
		using JHTDB<TypeVector, TypeMatrix, TypeRef>::ly;
		using JHTDB<TypeVector, TypeMatrix, TypeRef>::lz;
		using JHTDB<TypeVector, TypeMatrix, TypeRef>::nx;
		using JHTDB<TypeVector, TypeMatrix, TypeRef>::ny;
		using JHTDB<TypeVector, TypeMatrix, TypeRef>::nz;
		using JHTDB<TypeVector, TypeMatrix, TypeRef>::dt;
		using JHTDB<TypeVector, TypeMatrix, TypeRef>::tMax;
		// physical parameters
		using JHTDB<TypeVector, TypeMatrix, TypeRef>::viscosity;
		double h;
		double meanPressureGradient;
		// statistics
		double velocityBulk;
		double velocityCenterline;
		double velocityFriction;
		double viscousLengthScale;
};

}

}

#endif
