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

	int BoxX { }, BoxY { }, BoxZ { }, Steps_Equil { }, Steps_Total { }, Steps_Output { };
	double Temperature { }, Shear { }, StepSize { };
	bool AngularMomentumConservation { };
	stringstream ss_para { };


	//defaults für: Temperature, BoxSize(x, y, z), Shear rate, stepsize, step_aufwärm, step_total, step_output
	double a_para[]{0.5, 5, 5, 5, 0.5, 0.1, 1E3, 1E4, 100, 0};
	int a_para_size = sizeof(a_para) / sizeof(*a_para);
	int i_para{1};
	for (i_para = 1; i_para < min(a_para_size + 1, argc); ++i_para) {
		if (is_number(argv[i_para])) {
			a_para[i_para - 1] = stod(argv[i_para]);
		}
		else break;
	}

	Temperature = a_para[0];
	BoxX = a_para[1];
	BoxY = a_para[2];
	BoxZ = a_para[3];
	Shear = a_para[4];
	StepSize = a_para[5];
	Steps_Equil = (int)a_para[6];
	Steps_Total = (int)a_para[7];
	Steps_Output = (int)a_para[8];
	AngularMomentumConservation = (bool)a_para[9];

	Box box(BoxX, BoxY, BoxZ, Temperature, 1.0);



	MPC MPCroutine {box, Temperature, 10, Shear, AngularMomentumConservation};

	ss_para << "_N" << BoxX*BoxY*BoxZ*10;
	ss_para << "_Shear" << Shear;
	ss_para << "_h" << StepSize;
	if (AngularMomentumConservation) ss_para << "_AMC";

	ofstream fluid_file { };

	string fluid_file_name = "fluid" + ss_para.str() +".dat";
	fluid_file.open(fluid_file_name, ios::out | ios::trunc);


	std::cout << "Temperature: " << Temperature <<  std::endl;
	std::cout << "MPC is turned ON with shear rate: " << Shear << std::endl;
	std::cout << "Box Size: " << BoxX << " " << BoxY << " " << BoxZ << std::endl;
	std::cout << "Step Size: " << StepSize << std::endl;
	std::cout << "Total Steps: " << Steps_Total << " Equilibration: " << Steps_Equil << " Output every: " << Steps_Output << std::endl;


	MPCroutine.initialize();

	VelocityX velocity_average {0.2};
	nothing no_function { };
	clock_t begin = clock();

	for (int n = 0; n < Steps_Total; ++n) {
		MPCroutine.step(0.01);
		if (n > Steps_Equil && !(n%Steps_Output)) {
			std::cout << n << " " << MPCroutine.calculateCurrentTemperature() << std::endl;
			/*for (auto& part : MPCroutine.Fluid) {
				std::cout << part.Position.transpose() << " ; " << part.Velocity.transpose() << std::endl;
			}*/
			//MPCroutine(velocity_average, no_function);
			MPCroutine(velocity_average);
			//MPCroutine.unitary(velocity_average);
		}

	}
	velocity_average.print_result(fluid_file);

	clock_t end = clock();
	std::cout << "time: " << double(end-begin)/CLOCKS_PER_SEC << std::endl;
}



