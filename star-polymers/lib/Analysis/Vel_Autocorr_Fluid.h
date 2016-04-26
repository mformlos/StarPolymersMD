/*
 * Vel_Autocorr_Fluid.h
 *
 *  Created on: Apr 17, 2016
 *      Author: maud
 */

#ifndef STAR_POLYMERS_LIB_ANALYSIS_VEL_AUTOCORR_FLUID_H_
#define STAR_POLYMERS_LIB_ANALYSIS_VEL_AUTOCORR_FLUID_H_

#include "Analysis.h"
#include <map>
#include <iterator>
#include <../eigen/Eigen/Dense>




class Vel_Autocorr_Fluid: public Analysis<Molecule> {
protected:
	unsigned Entries;
	unsigned NumberOfParticles;
	double DeltaT;
	unsigned ParticleCounter;
	unsigned EntryCounter;
	unsigned Counter;
	std::vector<std::vector<Vector3d>> Velocities;
	std::vector<std::vector<long double>> Autocorrelation;
	std::vector<long double> Autocorrelation_Average;


public:
	Vel_Autocorr_Fluid(unsigned N_p, unsigned N, double dt) :
		Entries {N},
		NumberOfParticles {N_p},
		DeltaT {dt},
		ParticleCounter { },
		EntryCounter { },
		Counter { },
		Velocities{N_p, std::vector<Vector3d>(N, Vector3d::Zero())},
		Autocorrelation(N_p, std::vector<long double>(N,0.)),
		Autocorrelation_Average(N, 0.)
		{}

	void operator() (const Particle& part) {

		Velocities[ParticleCounter].erase(Velocities[ParticleCounter].begin());
		Velocities[ParticleCounter].push_back(part.Velocity);

		ParticleCounter++;

		if (ParticleCounter == NumberOfParticles && EntryCounter < Entries) EntryCounter++;

 		//std::cout << "Entry Counter " << EntryCounter <<  std::endl;
		/*for (auto& el : Velocities) {
			for (auto& vel : el) {
				std::cout << vel.transpose() << " ";
			}
			std::cout << "\n";
		}*/
		if (ParticleCounter == NumberOfParticles && EntryCounter == Entries) {
			//std::cout << "calculating... " << std::endl;
			calculate();
		}
		if (ParticleCounter == NumberOfParticles) ParticleCounter = 0;
	}

	void calculate() {
		for (unsigned i = 0; i < NumberOfParticles; i++) {
			for (unsigned j = 0; j < Entries; j++) {
				Autocorrelation[i][j] += Velocities[i][0].dot(Velocities[i][j]);
				//std::cout << Autocorrelation[i][j] << " ";
			}
			//std::cout << "\n";
		}
		Counter++;
	}

	void calculate_average() {
		for (unsigned i = 0; i < Entries; i++) {
			Autocorrelation_Average[i] = 0.;
			for (unsigned j = 0; j < NumberOfParticles; j++) {
				Autocorrelation_Average[i] += Autocorrelation[j][i]/Counter;
			}
			Autocorrelation_Average[i] /= NumberOfParticles;
		}
	}

	double value() {
		double Diffusion{};
		for (auto& el : Autocorrelation_Average) {
			Diffusion += el*DeltaT/Counter;
		}
		Diffusion /= 3.;
		return Diffusion;
	}

	std::ostream& print_result(std::ostream& os) {
		calculate_average();
		os << "#Diffusion coefficient = " << value() << "\n";
		for (unsigned i = 0; i < Entries; i++) {
			os << i*DeltaT << " " << Autocorrelation_Average[i] << "\n";
		}
		return os;
	}


};




#endif /* STAR_POLYMERS_LIB_ANALYSIS_VEL_AUTOCORR_FLUID_H_ */
