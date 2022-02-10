#ifndef FL0W_JHTDB_H
#define FL0W_JHTDB_H
#pragma once

#include "fl0w/flow.h"

// I did some really strange stuff here but it seems to work
#pragma push_macro("__cplusplus")
#undef __cplusplus
extern "C" {
	#include "turblib.h"
}
#pragma pop_macro("__cplusplus")

#include <thread>
#include <chrono>
#include <string>
#include <memory>
#include <unordered_map>

namespace fl0w {

namespace jhtdb {

template<typename TypeVector, typename TypeMatrix, template<typename...> class TypeRef>
class JHTDB : public Flow<TypeVector, TypeMatrix, TypeRef> {
	public:
		JHTDB(const std::string& dataset) : 
			sAuthtoken(std::make_shared<std::string>("edu.jhu.pha.turbulence.testing-201406")), 
			spatialInterp(Lag6), 
			spatialGradInterp(FD4Lag4), 
			temporalInterp(NoTInt), 
			attemptsNb(360), 
			sDataset(std::make_shared<std::string>(dataset)), 
			maxPointsNbPerQuery(4096)
		{
			// trublib
			::soapinit();
			::turblibSetExitOnError(0);
		}
		
		~JHTDB() {
			::soapdestroy();
		}
	public:
		TypeVector getVelocity(const TypeRef<const TypeVector>& x, const double& t) const override {
			std::string key = std::to_string(x[0]) + "_" + std::to_string(x[1]) + "_" + std::to_string(x[2]);
			if(preparedVelocities.count(key) == 1) {
				return preparedVelocities.at(key);
			} else {
				return queryVelocity(x, t);
			}
		}

		TypeMatrix getJacobian(const TypeRef<const TypeVector>& x, const double& t) const override {
			std::string key = std::to_string(x[0]) + "_" + std::to_string(x[1]) + "_" + std::to_string(x[2]);
			if(preparedJacobians.count(key) == 1) {
				return preparedJacobians.at(key);
			} else {
				return queryJacobian(x, t);
			}
		}

		TypeVector getAcceleration(const TypeRef<const TypeVector>& x, const double& t) const override {
			std::cout << "WARNING: JHTDB getAcceleration is not yet working" << std::endl;
			return TypeVector::Zero();
		}
	public:
		virtual void queryPointFromX(float* point, const TypeRef<const TypeVector>& x) const {
			// x
			if(x[0] < 0) {
				point[0] = lx - std::fmod(-x[0], lx);
			} else {
				point[0] = std::fmod(x[0], lx);
			}
			// y
			if(x[1] < 0) {
				point[1] = ly - std::fmod(-x[1], ly);
			} else {
				point[1] = std::fmod(x[1], ly);
			}
			// z
			if(x[2] < 0) {
				point[2] = lz - std::fmod(-x[2], lz);
			} else {
				point[2] = std::fmod(x[2], lz);
			}
		}
	public:
		TypeVector queryVelocity(const TypeRef<const TypeVector>& x, const double& t) const {
			// init
			float points[1][3];
			float result[1][3];
			// copy double vector x to float array points
			queryPointFromX(points[0], x);
			// get velocity into results
			int attempts = 0;
			while(::getVelocity(&((*sAuthtoken)[0]), &((*sDataset)[0]), t, Lag6, NoTInt, 1, points, result) != SOAP_OK) {
				if (attempts++ > attemptsNb) {
					std::cout << "ERROR: JHTDBFluid getVelocity Fatal Error: too many query failures." << std::endl;
					::exit(EXIT_FAILURE);
				} else {
					std::cout << "WARNING: JHTDBFluid getVelocity Error: " << ::turblibGetErrorString() << " . Trying again in 1s." << std::endl;
				}
				std::this_thread::sleep_for(std::chrono::seconds(1));
			}
			// copy float array results into vector u
			TypeVector u;
			for(unsigned int i = 0; i < 3; i++) {
				u[i] = result[0][i];
			}
			return u;
		}

		TypeMatrix queryJacobian(const TypeRef<const TypeVector>& x, const double& t) const {
			// init
			float points[1][3];
			float result[1][9];
			// copy double vector x to float array points
			queryPointFromX(points[0], x);
			// get gradient into results
			int attempts = 0;
			while(::getVelocityGradient(&((*sAuthtoken)[0]), &((*sDataset)[0]), t, FD4Lag4, NoTInt, 1, points, result) != SOAP_OK) {
				if (attempts++ > attemptsNb) {
					std::cout << "ERROR: JHTDBFluid getJacobian Fatal Error: too many query failures." << std::endl;
					::exit(EXIT_FAILURE);
				} else {
					std::cout << "WARNING: JHTDBFluid getJacobian Error: " << ::turblibGetErrorString() << " . Trying again in 1s." << std::endl;
				}
				std::this_thread::sleep_for(std::chrono::seconds(1));
			}
			// copy float array results into J
			TypeMatrix J;
			for(unsigned int i = 0; i < 3; i++) {
				for(unsigned int j = 0; j < 3; j++) {
					J(i, j) = result[0][i * 3 + j];
				}
			}
			return J;
		}
		
		std::vector<TypeVector> queryVelocities(const std::vector<TypeVector>& positions, const double& t) const {
			// init
			//float points[positions.size()][3];
			//float result[positions.size()][3];
			// c memory alloc
			float (*points)[3] = (float(*)[3]) std::malloc(3 * positions.size() * sizeof(float));
			float (*result)[3] = (float(*)[3]) std::malloc(3 * positions.size() * sizeof(float));
			// copy double vector x to float array points
			for(unsigned int i = 0; i < positions.size(); i++) {
				queryPointFromX(points[i], positions[i]);
			}
			// get velocity into results
			int startQueryIndex = 0;
			while(startQueryIndex < positions.size()) {
				// position pointers
				float (*points_tmp)[3] = (float(*)[3]) (((char*)points) + 3 * startQueryIndex * sizeof(float));
				float (*result_tmp)[3] = (float(*)[3]) (((char*)result) + 3 * startQueryIndex * sizeof(float));
				// query
				int attempts = 0;
				while(::getVelocity(&((*sAuthtoken)[0]), &((*sDataset)[0]), t, Lag6, NoTInt, std::min(positions.size() - startQueryIndex, maxPointsNbPerQuery), points_tmp, result_tmp) != SOAP_OK) {
					if (attempts++ > attemptsNb) {
						std::cout << "ERROR: JHTDBFluid getVelocity Fatal Error: too many query failures." << std::endl;
						::exit(EXIT_FAILURE);
					} else {
						std::cout << "WARNING: JHTDBFluid getVelocity Error: " << ::turblibGetErrorString() << " . Trying again in 1s." << std::endl;
					}
					std::this_thread::sleep_for(std::chrono::seconds(1));
				}
				startQueryIndex += maxPointsNbPerQuery;
			}
			// copy float array results into vector u
			std::vector<TypeVector> velocities(positions.size());
			for(unsigned int i = 0; i < velocities.size(); i++) {
				for(unsigned int j = 0; j < 3; j++) {
					velocities[i][j] = result[i][j];
				}
			}
			// c memory free
			std::free(points);
			std::free(result);
			// return
			return velocities;
		}

		std::vector<TypeMatrix> queryJacobians(const std::vector<TypeVector>& positions, const double& t) const {
			// init
			//float points[positions.size()][3];
			//float result[positions.size()][3];
			// c memory alloc
			float (*points)[3] = (float(*)[3]) std::malloc(3 * positions.size() * sizeof(float));
			float (*result)[9] = (float(*)[9]) std::malloc(9 * positions.size() * sizeof(float));
			// copy double vector x to float array points
			for(unsigned int i = 0; i < positions.size(); i++) {
				queryPointFromX(points[i], positions[i]);
			}
			// get velocity into results
			int startQueryIndex = 0;
			while(startQueryIndex < positions.size()) {
				// position pointers
				float (*points_tmp)[3] = (float(*)[3]) (((char*)points) + 3 * startQueryIndex * sizeof(float));
				float (*result_tmp)[9] = (float(*)[9]) (((char*)result) + 9 * startQueryIndex * sizeof(float));
				// query
				int attempts = 0;
				while(::getVelocityGradient(&((*sAuthtoken)[0]), &((*sDataset)[0]), t, FD4Lag4, NoTInt, std::min(positions.size() - startQueryIndex, maxPointsNbPerQuery), points_tmp, result_tmp) != SOAP_OK) {
					if (attempts++ > attemptsNb) {
						std::cout << "ERROR: JHTDBFluid getJacobian Fatal Error: too many query failures." << std::endl;
						::exit(EXIT_FAILURE);
					} else {
						std::cout << "WARNING: JHTDBFluid getJacobian Error: " << ::turblibGetErrorString() << " . Trying again in 1s." << std::endl;
					}
					std::this_thread::sleep_for(std::chrono::seconds(1));
				}
				startQueryIndex += maxPointsNbPerQuery;
			}
			// copy float array results into vector u
			std::vector<TypeMatrix> jacobians(positions.size());
			for(unsigned int i = 0; i < jacobians.size(); i++) {
				for(unsigned int j = 0; j < 3; j++) {
					for(unsigned int k = 0; k < 3; k++) {
						jacobians[i](j, k) = result[i][j * 3 + k];
					}
				}
			}
			// c memory free
			std::free(points);
			std::free(result);
			// return
			return jacobians;
		}
	public:
		void prepareVelocities(const std::vector<TypeVector>& positions, const double& t) {
			// get velocities
			std::vector<TypeVector> velocities = queryVelocities(positions, t);
			// set prepared velocities
			preparedVelocities.clear();
			for(unsigned int i = 0; i < positions.size(); i++) {
				preparedVelocities[std::to_string(positions[i][0]) + "_" + std::to_string(positions[i][1]) + "_" + std::to_string(positions[i][2])] = velocities[i];
			}
		}
		
		void prepareJacobians(const std::vector<TypeVector>& positions, const double& t) {
			// get jacobians
			std::vector<TypeMatrix> jacobians = queryJacobians(positions, t);
			// set prepared jacobians
			preparedJacobians.clear();
			for(unsigned int i = 0; i < positions.size(); i++) {
				preparedJacobians[std::to_string(positions[i][0]) + "_" + std::to_string(positions[i][1]) + "_" + std::to_string(positions[i][2])] = jacobians[i];
			}
		}
	public:
		// data
		std::unordered_map<std::string, TypeVector> preparedVelocities;
		std::unordered_map<std::string, TypeMatrix> preparedJacobians;
		// query
		std::shared_ptr<std::string> sAuthtoken;
		enum SpatialInterpolation spatialInterp;
		enum SpatialInterpolation spatialGradInterp;
		enum TemporalInterpolation temporalInterp;
		unsigned int attemptsNb;
		long unsigned int maxPointsNbPerQuery;
		// dataset
		std::shared_ptr<std::string> sDataset;
		// numerical parameters
        double lx;
        double ly;
        double lz;
        unsigned int nx;
        unsigned int ny;
        unsigned int nz;
        double dt;
        double tMax;
        // physical parameters
        double viscosity;
};

}

}

#endif
