/*
 * mpc_test.cpp
 *
 *  Created on: May 6, 2015
 *      Author: maud
 */

#include <iostream>
#include <fstream>
#include <forward_list>
#include <ctime>
#include <cmath>
#include <algorithm>
#include "Box.h"
#include "Thermostat_None.h"
#include "Lowe_Andersen.h"
#include "Nose_Hoover.h"
#include "Andersen.h"
#include "MPC.h"
#include "Analysis.h"
#include "Velocity_x.h"
#include "Vel_Autocorr_Fluid.h"
#include "Fourier_Autocorr.h"

bool is_number(const std::string &str)
{
	return str.find_first_not_of(".eE-0123456789") == std::string::npos;
}

double set_param(double def, char *array[], int length, int pos) {
	if (pos >= length) return def;
	if (is_number(array[pos])) return stod(array[pos]);
	return def;
}

struct nothing : std::binary_function<Particle, Particle, bool>
{
    bool operator()(const Particle& part1, const Particle& part2) const { return true; }
};

int main(int argc, char* argv[]) {

	int BoxX { }, BoxY { }, BoxZ { }, Steps_Equil { }, Steps_Total { }, Steps_Output { }, N { };
	double Temperature { }, Shear { }, StepSize { };
	bool AngularMomentumConservation { };
	stringstream ss_para { };


	//defaults für:  Particles, Temperature, BoxSize(x, y, z), Shear rate, stepsize, step_aufwärm, step_total, step_output, AMC
	double a_para[]{10, 0.5, 5, 5, 5, 0.5, 0.1, 1E3, 1E4, 100, 0};
	int a_para_size = sizeof(a_para) / sizeof(*a_para);
	int i_para{1};
	for (i_para = 1; i_para < min(a_para_size + 1, argc); ++i_para) {
		if (is_number(argv[i_para])) {
			a_para[i_para - 1] = stod(argv[i_para]);
		}
		else break;
	}

	N = a_para[0];
	Temperature = a_para[1];
	BoxX = a_para[2];
	BoxY = a_para[3];
	BoxZ = a_para[4];
	Shear = a_para[5];
	StepSize = a_para[6];
	Steps_Equil = (int)a_para[7];
	Steps_Total = (int)a_para[8];
	Steps_Output = (int)a_para[9];
	AngularMomentumConservation = (bool)a_para[10];

	Box box(BoxX, BoxY, BoxZ, Temperature, 1.0);



	MPC MPCroutine {box, Temperature, N, Shear, 1, AngularMomentumConservation};

	ss_para << "_N" << BoxX*BoxY*BoxZ*N;
	ss_para << "_Shear" << Shear;
	ss_para << "_h" << StepSize;
	if (AngularMomentumConservation) ss_para << "_AMC";

	ofstream fluid_file { };
	ofstream temperature_file { };
	ofstream autocorr_file { };
	ofstream fourier_file { };

	string fluid_file_name = "fluid" + ss_para.str() +".dat";
	string temperature_file_name = "temperature" + ss_para.str() + ".dat";
	string autocorr_file_name = "velocity" + ss_para.str() + ".dat";
	string fourier_file_name = "fourier" + ss_para.str() + ".dat";

	fluid_file.open(fluid_file_name, ios::out | ios::trunc);
	temperature_file.open(temperature_file_name, ios::out | ios::trunc);
	autocorr_file.open(autocorr_file_name, ios::out | ios::trunc);
	fourier_file.open(fourier_file_name, ios::out | ios::trunc);

	std::cout << "Temperature: " << Temperature <<  std::endl;
	std::cout << "MPC is turned ON with shear rate: " << Shear << std::endl;
	std::cout << "Box Size: " << BoxX << " " << BoxY << " " << BoxZ << std::endl;
	std::cout << "Step Size: " << StepSize << std::endl;
	std::cout << "Total Steps: " << Steps_Total << " Equilibration: " << Steps_Equil << " Output every: " << Steps_Output << std::endl;


	MPCroutine.initialize();

	VelocityX velocity_average {box, MPCroutine, 0.2};
	Vel_Autocorr_Fluid autocorr {MPCroutine.NumberOfParticles(), 50, 1.0};
	Fourier_Autocorr fourier {MPCroutine.NumberOfParticles(), 50, 1.0, (double)BoxX};
	clock_t begin = clock();

	for (int n = 0; n < Steps_Total; ++n) {
		MPCroutine.step(n, StepSize);
		if (n > Steps_Equil && !(n%Steps_Output)) {
			std::cout << n << " " << MPCroutine.calculateCurrentTemperature() << std::endl;
			/*for (auto& part : MPCroutine.Fluid) {
				std::cout << part.Position.transpose() << " ; " << part.Velocity.transpose() << std::endl;
			}*/
			MPCroutine(velocity_average);
			//std::cout << "autocorrelation start" << std::endl;
			//std::cout << "autocorrelation end" << std::endl;
			temperature_file << n << " " << MPCroutine.calculateCurrentTemperature() << std::endl;
			//MPCroutine.unitary(velocity_average);
		}
		if (n > Steps_Equil && !(n%10)){
			MPCroutine(autocorr);
			MPCroutine(fourier);
		}

	}
	velocity_average.print_result(fluid_file);
	autocorr.print_result(autocorr_file);
	fourier.print_result(fourier_file);

	clock_t end = clock();
	std::cout << "time: " << double(end-begin)/CLOCKS_PER_SEC << std::endl;
}



