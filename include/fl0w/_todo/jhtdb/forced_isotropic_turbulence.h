#ifndef FL0W_JHTDB_FORCED_ISOTROPIC_TURBULENCE_H
#define FL0W_JHTDB_FORCED_ISOTROPIC_TURBULENCE_H
#pragma once

#include "fl0w/jhtdb/jhtdb.h"

namespace fl0w {

namespace jhtdb {

template<typename TypeVector, typename TypeMatrix, template<typename...> class TypeRef>
class ForcedIsotropicTurbulence : public JHTDB<TypeVector, TypeMatrix, TypeRef> {
	public:
		ForcedIsotropicTurbulence() : JHTDB<TypeVector, TypeMatrix, TypeRef>::JHTDB("isotropic1024coarse") {
			// numerical parameters
			l = 2.0 * M_PI;
			lx = l;
			ly = l;
			lz = l;
			n = 1024;
			nx = n;
			ny = n;
			nz = n;
			dt = 0.002;
			tMax = 10.056;
			// physical parameters
			viscosity = 0.000185;
			// statistics
			dissipation = 0.103;
			rmsVelocity = 0.686;
			taylorMicroScale = 0.113;
			taylorScaleReynolds = 418.0;
			kolmogorovTimeScale = 0.0424;
			kolmogorovLengthScale = 0.00280;
			integralLengthScale = 1.364;
			largeEddyTurnOverTime = 1.99;
		}
	public:
		// numerical parameters
		double l;
		using JHTDB<TypeVector, TypeMatrix, TypeRef>::lx;
		using JHTDB<TypeVector, TypeMatrix, TypeRef>::ly;
		using JHTDB<TypeVector, TypeMatrix, TypeRef>::lz;
		unsigned int n;
		using JHTDB<TypeVector, TypeMatrix, TypeRef>::nx;
		using JHTDB<TypeVector, TypeMatrix, TypeRef>::ny;
		using JHTDB<TypeVector, TypeMatrix, TypeRef>::nz;
		using JHTDB<TypeVector, TypeMatrix, TypeRef>::dt;
		using JHTDB<TypeVector, TypeMatrix, TypeRef>::tMax;
		// physical parameters
		using JHTDB<TypeVector, TypeMatrix, TypeRef>::viscosity;
		// statistics
		double dissipation;
		double rmsVelocity;
		double taylorMicroScale;
		double taylorScaleReynolds;
		double kolmogorovTimeScale;
		double kolmogorovLengthScale;
		double integralLengthScale;
		double largeEddyTurnOverTime;
};

}

}

#endif
