/*
 * Vel_Autocorr.h
 *
 *  Created on: Apr 12, 2016
 *      Author: maud
 */

#ifndef STAR_POLYMERS_LIB_ANALYSIS_VEL_AUTOCORR_H_
#define STAR_POLYMERS_LIB_ANALYSIS_VEL_AUTOCORR_H_



#include "Analysis.h"
#include <map>
#include <iterator>
#include <../eigen/Eigen/Dense>




class Vel_Autocorr: public Analysis<Molecule> {
protected:
	unsigned Entries;
	double DeltaT;
	unsigned EntryCounter;
	unsigned Counter;
	long double MeanSqVel;
	std::vector<Vector3d> Velocities;
	std::vector<long double> Autocorrelation;


public:
	Vel_Autocorr(unsigned N, double dt) :
		Entries {N},
		DeltaT {dt},
		EntryCounter{ },
		Counter { },
		MeanSqVel { },
		Velocities{},
		Autocorrelation(N, 0.) {}

	void operator() (const Molecule& mol) {
		Vector3d Velocity = mol.calculate_center_of_mass_velocity();
		if (EntryCounter == Entries) {
			Velocities.erase(Velocities.begin());
			Velocities.push_back(Velocity);
			for (unsigned i = 0; i < Entries; i++) {
				Autocorrelation[i] += Velocities[0].dot(Velocities[i]);
			}
			MeanSqVel += Velocity.dot(Velocity);
			Counter++;
		}
		else {
			Velocities.push_back(Velocity);
			Counter++;
		}
	}

	double value() {
		double Diffusion{};
		for (auto& el : Autocorrelation) {
			Diffusion += el*DeltaT/Counter;
		}
		Diffusion /= 3.;
		return Diffusion;
	}

	std::ostream& print_result(std::ostream& os) {
		os << "#Diffusion coefficient = " << value() << "\n";
		double MeanSqVel_result = MeanSqVel / Counter;
		for (unsigned i = 0; i < Entries; i++) {
			os << i*DeltaT << " " << Autocorrelation[i]/Counter << " " << MeanSqVel_result << "\n";
		}
		return os;
	}


};


#endif /* STAR_POLYMERS_LIB_ANALYSIS_VEL_AUTOCORR_H_ */
