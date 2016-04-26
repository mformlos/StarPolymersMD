/*
 * Fourier_Autocorr.h
 *
 *  Created on: Apr 18, 2016
 *      Author: maud
 */

#ifndef STAR_POLYMERS_LIB_ANALYSIS_FOURIER_AUTOCORR_H_
#define STAR_POLYMERS_LIB_ANALYSIS_FOURIER_AUTOCORR_H_

#include "Analysis.h"
#include <map>
#include <iterator>
#include <../eigen/Eigen/Dense>
#include <complex>

using namespace Eigen;

class Fourier_Autocorr: public Analysis<Molecule> {
protected:
	unsigned Entries;
	unsigned NumberOfParticles;
	double DeltaT;
	unsigned ParticleCounter;
	unsigned EntryCounter;
	unsigned Counter;
	std::vector<double> Nk;
    std::vector<Vector3cd> k;
    std::vector<Vector3cd> Velocity_Fourier;
	std::vector<std::vector<Vector3cd>> Velocities;
	std::vector<std::vector<long double>> Autocorrelation;


public:
	Fourier_Autocorr(unsigned N_p, unsigned N, double dt, double L) :
		Entries {N},
		NumberOfParticles {N_p},
		DeltaT {dt},
		ParticleCounter { },
		EntryCounter { },
		Counter { },
		Nk(3, 0.),
		k(3, Vector3cd::Zero()),
		Velocity_Fourier(3, Vector3cd::Zero()),
		Velocities(3, std::vector<Vector3cd>(N, Vector3cd::Zero())),
		Autocorrelation(3, std::vector<long double>(N,0.))
		{
			for (unsigned i = 0; i < 3; i++) {
				Nk[i] = 2*M_PI/L;
				k[i](i) = 1.;
				std::cout << k[i].transpose() << std::endl;
			}
		}

	void operator() (const Particle& part) {
		if (ParticleCounter == 0) {
			for (auto& el : Velocity_Fourier) el = Vector3cd::Zero();
		}
		std::complex<double> I(0.,1.);
		for(unsigned j = 0; j < 3; j++) {
			Velocity_Fourier[j] += part.Velocity*exp(I*Nk[j]*k[j].dot(part.Position));
		}

		ParticleCounter++;
		if (ParticleCounter == NumberOfParticles) {
			for (unsigned j = 0; j < 3; j++) {
				Velocity_Fourier[j] /= NumberOfParticles;
				Velocity_Fourier[j] -= k[j]*(k[j].dot(Velocity_Fourier[j]));
				Velocities[j].erase(Velocities[j].begin());
				Velocities[j].push_back(Velocity_Fourier[j]);
			}
		}

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
		for (unsigned i = 0; i < 3; i++) {
			for (unsigned j = 0; j < Entries; j++) {
				Autocorrelation[i][j] += (Velocities[i][0].dot((Velocities[i][j]))).real();
				//std::cout << Autocorrelation[i][j] << " ";
			}
			//std::cout << "\n";
		}
		Counter++;
	}



	/*double value() {
		double Diffusion{};
		for (auto& el : Autocorrelation) {
			Diffusion += el*DeltaT/Counter;
		}
		Diffusion /= 3.;
		return Diffusion;
	}*/

	std::ostream& print_result(std::ostream& os) {
		//os << "#Diffusion coefficient = " << value() << "\n";
		for (unsigned i = 0; i < Entries; i++) {
			os << i*DeltaT << " " << Autocorrelation[0][i]/Counter << " "  << Autocorrelation[1][i]/Counter << " " << Autocorrelation[2][i]/Counter << "\n";
		}
		return os;
	}

};



#endif /* STAR_POLYMERS_LIB_ANALYSIS_FOURIER_AUTOCORR_H_ */
