#ifndef FL0W_JHTDB_H
#define FL0W_JHTDB_H
#pragma once

// #pragma push_macro("__cplusplus")
// #undef __cplusplus
// extern "C" {
// 	#include "turblib.h"
// }
// #pragma pop_macro("__cplusplus")

// thirdparty
#include "getjhtdb/getData.h"
// std
#include <string>
#include <mutex>

#include <chrono>
#include <thread>

namespace fl0w {

template<typename tVector, typename tMatrix, template<typename...> class tView>
struct Jhtdb {

	// query
	inline static std::string authToken = "edu.jhu.pha.turbulence.testing-201406";
	inline static std::string velocitySpatialInterpolation = "lag6";
	inline static std::string gradientsSpatialInterpolation = "fd4lag4";
	inline static std::string temporalInterpolation = "pchip";
	static constexpr unsigned int maxAttemptNumber = 256;
	static constexpr long unsigned int maxPointNumberPerQuery = 4096;

	// base queries

	static tVector queryVelocity(std::string& dataset, const float* pX, const double t) {
		// init
		float points[1][3];
		for(unsigned int i(0); i < 3; ++i) {
			points[0][i] = pX[i];
		}

		// query
		const std::string variable = "velocity";
		const std::string spatialOperator = "field";

		double **output = getData(authToken.data(), dataset.data(), t, velocitySpatialInterpolation.data(), temporalInterpolation.data(), variable.data(), spatialOperator.data(), points, 1);

		// copy float array output into vector
		tVector vOutput;
		for(unsigned int i(0); i < 3; ++i) {
			vOutput(i) = output[0][i];
		}

		// clean
		free(output);
		
		return vOutput;
	}

	static double** queryVelocityArray(std::string& dataset, float (*points)[3], const unsigned int n, const double t) {

		// query
		const std::string variable = "velocity";
		const std::string spatialOperator = "field";

		double **output = getData(authToken.data(), dataset.data(), t, velocitySpatialInterpolation.data(), temporalInterpolation.data(), variable.data(), spatialOperator.data(), points, n);
		return output;
	}

	static tMatrix queryVelocityGradients(std::string& dataset, const float* pX, const double t) {
		// init
		float points[1][3];
		for(unsigned int i(0); i < 3; ++i) {
			points[0][i] = pX[i];
		}
		
		// query
		const std::string variable = "velocity";
		const std::string spatialOperator = "gradient";
		
		double **output = getData(authToken.data(), dataset.data(), t, gradientsSpatialInterpolation.data(), temporalInterpolation.data(), variable.data(), spatialOperator.data(), points, 1);

		// copy float array output into matrix
		tMatrix mOutput;
		for(unsigned int i(0); i < 3; ++i) {
			for(unsigned int j(0); j < 3; ++j) {
				mOutput(i, j) = output[0][i * 3 + j];
			}
		}

		// clean
		free(output);
		
		return mOutput;
	}

	static double** queryVelocityGradientsArray(std::string& dataset, float (*points)[3], const unsigned int n, const double t) {
		// query
		const std::string variable = "velocity";
		const std::string spatialOperator = "gradient";

		double **output = getData(authToken.data(), dataset.data(), t, gradientsSpatialInterpolation.data(), temporalInterpolation.data(), variable.data(), spatialOperator.data(), points, n);
		return output;
	}

	struct Isotropic {
		// preparation
		inline static std::unordered_map<double, unsigned int> preparedVelocityPointIndex;
		inline static std::unordered_map<double, unsigned int> preparedVelocityGradientsPointIndex;
		inline static std::unordered_map<double, float (*)[3]> preparedVelocityPoints;
		inline static std::unordered_map<double, float (*)[3]> preparedVelocityGradientsPoints;
		// data
		inline static std::unordered_map<double, std::unordered_map<std::string, tVector>> preparedVelocity;
		inline static std::unordered_map<double, std::unordered_map<std::string, tMatrix>> preparedVelocityGradients;
		// dataset
		inline static std::string dataset = "isotropic1024coarse";
		// numerical parameters
		static constexpr double lx = 2.0 * M_PI;
		static constexpr double ly = 2.0 * M_PI;
		static constexpr double lz = 2.0 * M_PI;
		static constexpr unsigned int nx = 1024;
		static constexpr unsigned int ny = 1024;
		static constexpr unsigned int nz = 1024;
		inline static const double dt = 0.002;
		inline static const double tMax = 10.056;
		// physical parameters
		inline static const double viscosity = 0.000185;
		// statistics
		inline static const double dissipation = 0.103;
		inline static const double rmsVelocity = 0.686;
		inline static const double taylorMicroScale = 0.113;
		inline static const double taylorScaleReynolds = 418.0;
		inline static const double kolmogorovTimeScale = 0.0424;
		inline static const double kolmogorovLengthScale = 0.00280;
		inline static const double integralLengthScale = 1.364;
		inline static const double largeEddyTurnOverTime = 1.99;
		// mutex
		inline static std::mutex preparedVelocityMutex;
		inline static std::mutex updateVelocityMutex;
		
		inline static std::mutex preparedVelocityGradientsMutex;
		inline static std::mutex updateVelocityGradientsMutex;

		// get

		static tVector getVelocity(const double* pX, const double t) {
			// periodicity
			float pXPeriodic[3];
			xToXPeriodic(pX, pXPeriodic);
			const std::string xPeriodicKey = getKeyFromPoint(pXPeriodic);
			// query
			return preparedVelocity.at(t).at(xPeriodicKey);
		}
	
		static tMatrix getVelocityGradients(const double* pX, const double t) {
			// periodicity
			float pXPeriodic[3];
			xToXPeriodic(pX, pXPeriodic);
			const std::string xPeriodicKey = getKeyFromPoint(pXPeriodic);
			// query
			return preparedVelocityGradients.at(t).at(xPeriodicKey);
		}

		// prepare

		static void prepare(const double t) {
			updatePreparedVelocity(t);
			updatePreparedVelocityGradients(t);
		}

		static void prepareVelocity(const double* pX, const double t) {
			// periodicity
			float pXPeriodic[3];
			xToXPeriodic(pX, pXPeriodic);
			const std::string xPeriodicKey = getKeyFromPoint(pXPeriodic);
			if (preparedVelocity.count(t) == 0 || preparedVelocity.at(t).count(xPeriodicKey) == 0) {
				// lock
				preparedVelocityMutex.lock();
				// prepare
				if(preparedVelocityPoints.count(t) == 0) {
					preparedVelocityPointIndex[t] = 0;
					preparedVelocityPoints[t] = new float[maxPointNumberPerQuery][3]; // TODO: Free ?
				}
				preparedVelocityPoints[t][preparedVelocityPointIndex[t]][0] = pXPeriodic[0];
				preparedVelocityPoints[t][preparedVelocityPointIndex[t]][1] = pXPeriodic[1];
				preparedVelocityPoints[t][preparedVelocityPointIndex[t]][2] = pXPeriodic[2];
				preparedVelocityPointIndex[t] += 1;
				// init prepared
				preparedVelocity[t][xPeriodicKey] = tVector::Zero();
				// if prepared
				if(preparedVelocityPointIndex[t] >= maxPointNumberPerQuery){
					updatePreparedVelocity(t, true);
				}
				// unlock
				preparedVelocityMutex.unlock();
			}
		}

		static void updatePreparedVelocity(const double t, const bool isLocked = false) {
			updateVelocityMutex.lock();
			
			if(not isLocked) {
				preparedVelocityMutex.lock();
			}
			
			long unsigned int pointNumber = preparedVelocityPointIndex[t];
			float (*points)[3] = new float[pointNumber][3];
			for(unsigned int i = 0; i < pointNumber; ++i) {
				for(unsigned int j = 0; j < 3; ++j) {
					points[i][j] = preparedVelocityPoints[t][i][j];
				}
			}
			preparedVelocityPointIndex[t] = 0;

			if(not isLocked) {
				preparedVelocityMutex.unlock();
			}
			
			if (pointNumber > 0) {
				double** velocityArray = queryVelocityArray(dataset, points, pointNumber, t);
				for(unsigned int index = 0; index < pointNumber; ++index) {
					const std::string pointKey = getKeyFromPoint(points[index]);
					preparedVelocity[t][pointKey] = tView<tVector>(velocityArray[index]);
				}
				free(velocityArray);
			}
			free(points);

			updateVelocityMutex.unlock();
		}

		static void prepareVelocityGradients(const double* pX, const double t) {
			// periodicity
			float pXPeriodic[3];
			xToXPeriodic(pX, pXPeriodic);
			const std::string xPeriodicKey = getKeyFromPoint(pXPeriodic);
			if (preparedVelocityGradients.count(t) == 0 || preparedVelocityGradients.at(t).count(xPeriodicKey) == 0) {
				// lock
				preparedVelocityGradientsMutex.lock();
				// prepare
				if(preparedVelocityGradientsPoints.count(t) == 0) {
					preparedVelocityGradientsPointIndex[t] = 0;
					preparedVelocityGradientsPoints[t] = new float[maxPointNumberPerQuery][3]; // TODO: Free ?
				}
				preparedVelocityGradientsPoints[t][preparedVelocityGradientsPointIndex[t]][0] = pXPeriodic[0];
				preparedVelocityGradientsPoints[t][preparedVelocityGradientsPointIndex[t]][1] = pXPeriodic[1];
				preparedVelocityGradientsPoints[t][preparedVelocityGradientsPointIndex[t]][2] = pXPeriodic[2];
				preparedVelocityGradientsPointIndex[t] += 1;
				// init prepared
				preparedVelocity[t][xPeriodicKey] = tVector::Zero();
				// if prepared
				if(preparedVelocityGradientsPointIndex[t] >= maxPointNumberPerQuery){
					updatePreparedVelocityGradients(t);
				}
				// unlock
				preparedVelocityGradientsMutex.unlock();
			}
		}

		static void updatePreparedVelocityGradients(const double t) {
			preparedVelocityGradientsMutex.lock();
			preparedVelocityGradientsMutex.unlock();
			
			if(preparedVelocityGradientsPointIndex[t] > 0) {
				double** velocityGradientsArray = queryVelocityGradientsArray(dataset, preparedVelocityGradientsPoints[t], preparedVelocityGradientsPointIndex[t], t);
				for(unsigned int index = 0; index < preparedVelocityGradientsPointIndex[t]; ++index) {
					const std::string pointKey = getKeyFromPoint(preparedVelocityGradientsPoints[t][index]);
					preparedVelocityGradients[t][pointKey] = tView<tMatrix>(velocityGradientsArray[index]);
				}
				free(velocityGradientsArray);
				// restart preparation
				preparedVelocityGradientsPointIndex[t] = 0;
			}
		}

		// internal

		static void xToXPeriodic(const double* pX, float* pXPeriodic) {
			// x
			if(pX[0] < 0) {
				pXPeriodic[0] = lx - std::fmod(-pX[0], lx);
			} else {
				pXPeriodic[0] = std::fmod(pX[0], lx);
			}
			// y
			if(pX[1] < 0) {
				pXPeriodic[1] = ly - std::fmod(-pX[1], ly);
			} else {
				pXPeriodic[1] = std::fmod(pX[1], ly);
			}
			// z
			if(pX[2] < 0) {
				pXPeriodic[2] = lz - std::fmod(-pX[2], lz);
			} else {
				pXPeriodic[2] = std::fmod(pX[2], lz);
			}
		}
		
	};

	// iternal
	
	static std::string getKeyFromPoint(const float* pX) {
		return std::to_string(pX[0]) + "_" + std::to_string(pX[1]) + "_" + std::to_string(pX[2]);
	}

	static tVector getPointFromKey(const std::string& key) {
		// init
		tVector output;
		// util
		unsigned int start;
		unsigned int end;
		// x0
		start = 0;
		end = key.find('_');
		output[0] = float(key.substr(start, end));
		// x1
		start = end + 1;
		end = key.find(start, '_');
		output[1] = float(key.substr(start, end));
		// x2
		start = end + 1;
		end = key.find(start, '_');
		output[2] = float(key.substr(start, end));
		// output
		return output;
	}

	// batch queries

// 	// filtered queries
// 
// 		TypeVector queryVelocityFiltered(const TypeRef<const TypeVector>& x, const double& t, const double& width) const {
// 			// init
// 			float points[1][3];
// 			float result[1][3];
// 			// copy double vector x to float array points
// 			queryPointFromX(points[0], x);
// 			// get velocity into results
// 			int attempts = 0;
// 			while(::getBoxFilter(&((*sAuthtoken)[0]), &((*sDataset)[0]), &((*sField)[0]), t, width, 1, points, result) != SOAP_OK) {
// 				if (attempts++ > attemptsNb) {
// 					std::cout << "ERROR: JHTDBFluid getVelocity Fatal Error: too many query failures." << std::endl;
// 					::exit(EXIT_FAILURE);
// 				} else {
// 					std::cout << "WARNING: JHTDBFluid getVelocity Error: " << ::turblibGetErrorString() << " . Trying again in 1s." << std::endl;
// 				}
// 				std::this_thread::sleep_for(std::chrono::seconds(1));
// 			}
// 			// copy float array results into vector u
// 			TypeVector u;
// 			for(unsigned int i = 0; i < 3; i++) {
// 				u[i] = result[0][i];
// 			}
// 			return u;
// 		}
// 
// 		TypeMatrix queryVelocityGradientsFiltered(const TypeRef<const TypeVector>& x, const double& t, const double& width) const {
// 			// init
// 			float points[1][3];
// 			float result[1][9];
// 			// copy double vector x to float array points
// 			queryPointFromX(points[0], x);
// 			// get gradient into results
// 			int attempts = 0;
// 			while(::getBoxFilterGradient(&((*sAuthtoken)[0]), &((*sDataset)[0]), &((*sField)[0]), t, width, 0.25 * width, 1, points, result) != SOAP_OK) {
// 				if (attempts++ > attemptsNb) {
// 					std::cout << "ERROR: JHTDBFluid getVelocityGradients Fatal Error: too many query failures." << std::endl;
// 					::exit(EXIT_FAILURE);
// 				} else {
// 					std::cout << "WARNING: JHTDBFluid getVelocityGradients Error: " << ::turblibGetErrorString() << " . Trying again in 1s." << std::endl;
// 				}
// 				std::this_thread::sleep_for(std::chrono::seconds(1));
// 			}
// 			// copy float array results into J
// 			TypeMatrix J;
// 			for(unsigned int i = 0; i < 3; i++) {
// 				for(unsigned int j = 0; j < 3; j++) {
// 					J(i, j) = result[0][i * 3 + j];
// 				}
// 			}
// 			return J;
// 		}
		
// 		std::vector<TypeVector> queryVelocities(const std::vector<TypeVector>& positions, const double& t) const {
// 			// init
// 			//float points[positions.size()][3];
// 			//float result[positions.size()][3];
// 			// c memory alloc
// 			float (*points)[3] = (float(*)[3]) std::malloc(3 * positions.size() * sizeof(float));
// 			float (*result)[3] = (float(*)[3]) std::malloc(3 * positions.size() * sizeof(float));
// 			// copy double vector x to float array points
// 			for(unsigned int i = 0; i < positions.size(); i++) {
// 				queryPointFromX(points[i], positions[i]);
// 			}
// 			// get velocity into results
// 			int startQueryIndex = 0;
// 			while(startQueryIndex < positions.size()) {
// 				// position pointers
// 				float (*points_tmp)[3] = (float(*)[3]) (((char*)points) + 3 * startQueryIndex * sizeof(float));
// 				float (*result_tmp)[3] = (float(*)[3]) (((char*)result) + 3 * startQueryIndex * sizeof(float));
// 				// query
// 				int attempts = 0;
// 				while(::getVelocity(&((*sAuthtoken)[0]), &((*sDataset)[0]), t, Lag6, NoTInt, std::min(positions.size() - startQueryIndex, maxPointsNbPerQuery), points_tmp, result_tmp) != SOAP_OK) {
// 					if (attempts++ > attemptsNb) {
// 						std::cout << "ERROR: JHTDBFluid getVelocity Fatal Error: too many query failures." << std::endl;
// 						::exit(EXIT_FAILURE);
// 					} else {
// 						std::cout << "WARNING: JHTDBFluid getVelocity Error: " << ::turblibGetErrorString() << " . Trying again in 1s." << std::endl;
// 					}
// 					std::this_thread::sleep_for(std::chrono::seconds(1));
// 				}
// 				startQueryIndex += maxPointsNbPerQuery;
// 			}
// 			// copy float array results into vector u
// 			std::vector<TypeVector> velocities(positions.size());
// 			for(unsigned int i = 0; i < velocities.size(); i++) {
// 				for(unsigned int j = 0; j < 3; j++) {
// 					velocities[i][j] = result[i][j];
// 				}
// 			}
// 			// c memory free
// 			std::free(points);
// 			std::free(result);
// 			// return
// 			return velocities;
// 		}
// 
// 		std::vector<TypeMatrix> queryVelocityGradientss(const std::vector<TypeVector>& positions, const double& t) const {
// 			// init
// 			//float points[positions.size()][3];
// 			//float result[positions.size()][3];
// 			// c memory alloc
// 			float (*points)[3] = (float(*)[3]) std::malloc(3 * positions.size() * sizeof(float));
// 			float (*result)[9] = (float(*)[9]) std::malloc(9 * positions.size() * sizeof(float));
// 			// copy double vector x to float array points
// 			for(unsigned int i = 0; i < positions.size(); i++) {
// 				queryPointFromX(points[i], positions[i]);
// 			}
// 			// get velocity into results
// 			int startQueryIndex = 0;
// 			while(startQueryIndex < positions.size()) {
// 				// position pointers
// 				float (*points_tmp)[3] = (float(*)[3]) (((char*)points) + 3 * startQueryIndex * sizeof(float));
// 				float (*result_tmp)[9] = (float(*)[9]) (((char*)result) + 9 * startQueryIndex * sizeof(float));
// 				// query
// 				int attempts = 0;
// 				while(::getVelocityGradient(&((*sAuthtoken)[0]), &((*sDataset)[0]), t, FD4Lag4, NoTInt, std::min(positions.size() - startQueryIndex, maxPointsNbPerQuery), points_tmp, result_tmp) != SOAP_OK) {
// 					if (attempts++ > attemptsNb) {
// 						std::cout << "ERROR: JHTDBFluid getVelocityGradients Fatal Error: too many query failures." << std::endl;
// 						::exit(EXIT_FAILURE);
// 					} else {
// 						std::cout << "WARNING: JHTDBFluid getVelocityGradients Error: " << ::turblibGetErrorString() << " . Trying again in 1s." << std::endl;
// 					}
// 					std::this_thread::sleep_for(std::chrono::seconds(1));
// 				}
// 				startQueryIndex += maxPointsNbPerQuery;
// 			}
// 			// copy float array results into vector u
// 			std::vector<TypeMatrix> jacobians(positions.size());
// 			for(unsigned int i = 0; i < jacobians.size(); i++) {
// 				for(unsigned int j = 0; j < 3; j++) {
// 					for(unsigned int k = 0; k < 3; k++) {
// 						jacobians[i](j, k) = result[i][j * 3 + k];
// 					}
// 				}
// 			}
// 			// c memory free
// 			std::free(points);
// 			std::free(result);
// 			// return
// 			return jacobians;
// 		}
// 	public:
// 		void prepareVelocities(const std::vector<TypeVector>& positions, const double& t) {
// 			// get velocities
// 			std::vector<TypeVector> velocities = queryVelocities(positions, t);
// 			// set prepared velocities
// 			preparedVelocities.clear();
// 			for(unsigned int i = 0; i < positions.size(); i++) {
// 				preparedVelocities[std::to_string(positions[i][0]) + "_" + std::to_string(positions[i][1]) + "_" + std::to_string(positions[i][2])] = velocities[i];
// 			}
// 		}
// 		
// 		void prepareVelocityGradientss(const std::vector<TypeVector>& positions, const double& t) {
// 			// get jacobians
// 			std::vector<TypeMatrix> jacobians = queryVelocityGradientss(positions, t);
// 			// set prepared jacobians
// 			preparedVelocityGradientss.clear();
// 			for(unsigned int i = 0; i < positions.size(); i++) {
// 				preparedVelocityGradientss[std::to_string(positions[i][0]) + "_" + std::to_string(positions[i][1]) + "_" + std::to_string(positions[i][2])] = jacobians[i];
// 			}
// 		}
};

}

#endif
