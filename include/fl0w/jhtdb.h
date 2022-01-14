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

template<typename TypeVector, typename TypeMatrix, template<typename...> class TypeRef>
class JHTDB : public Flow<TypeVector, TypeMatrix, TypeRef> {
    public:
        JHTDB() : sAuthtoken(std::make_shared<std::string>("edu.jhu.pha.turbulence.testing-201406")), spatialInterp(Lag6), spatialGradInterp(FD4Lag4), temporalInterp(NoTInt), attemptsNb(360), sDataset(std::make_shared<std::string>("")), maxPointsNbPerQuery(4096) {
            selectDataset("isotropic1024coarse");
            ::soapinit();
            ::turblibSetExitOnError(0);
        }
        ~JHTDB() {
            ::soapdestroy();
        }

        void selectDataset(const std::string& dataset) {
            *sDataset = dataset;
            if(dataset == "isotropic1024coarse") {
                // simulation properies
                domainLength = 2.0 * M_PI;
                gridSize = 1024;
                viscosity = 0.000185;
                dt = 0.002; // saved data step
                tMax = 10.056;
                // turbulence statistics
                dissipation = 0.103;
                rmsVelocity = 0.686;
                taylorMicroScale = 0.113;
                taylorScaleReynolds = 418.0;
                kolmogorovTimeScale = 0.0424;
                kolmogorovLengthScale = 0.00280;
                integralLengthScale = 1.364;
                largeEddyTurnOverTime = 1.99;
            } else {
                std::cout << "WARNING: unknown JHTDB flow dataset" << std::endl;
            }
        }
        
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
        TypeVector queryVelocity(const TypeRef<const TypeVector>& x, const double& t) const {
            // init
            float points[1][3];
            float result[1][3];
            // copy double vector x to float array points
            for(unsigned int i = 0; i < 3; i++) {
                points[0][i] = std::fmod(std::abs(x[i]), domainLength);
                if (x[i] < 0) {
                    points[0][i] = domainLength - points[0][i];
                }
            }
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
            for(unsigned int i = 0; i < 3; i++) {
                points[0][i] = std::fmod(std::abs(x[i]), domainLength);
                if (x[i] < 0) {
                    points[0][i] = domainLength - points[0][i];
                }
            }
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
                for(unsigned int j = 0; j < 3; j++) {
                    points[i][j] = std::fmod(std::abs(positions[i][j]), domainLength);
                    if (positions[i][j] < 0) {
                        points[i][j] = domainLength - points[i][j];
                    }
                }
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
                for(unsigned int j = 0; j < 3; j++) {
                    points[i][j] = std::fmod(std::abs(positions[i][j]), domainLength);
                    if (positions[i][j] < 0) {
                        points[i][j] = domainLength - points[i][j];
                    }
                }
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
        // simulation properties
        double domainLength;
        unsigned int gridSize;
        double viscosity;
        double dt;
        double tMax;
        // turbulence statistics
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

#endif
