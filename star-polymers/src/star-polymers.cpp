//============================================================================
// Name        : star-polymers.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C, Ansi-style
//============================================================================

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

bool is_number(const std::string &str)
{
	return str.find_first_not_of(".eE-0123456789") == std::string::npos;
}

double set_param(double def, char *array[], int length, int pos) {
	if (pos >= length) return def;
	if (is_number(array[pos])) return stod(array[pos]);
	return def;
}

int main(int argc, char* argv[]) {

	int BoxX { }, BoxY { }, BoxZ { }, TypeA { }, TypeB { }, Steps_Equil { }, Steps_Total { }, Steps_Output { };
	double Temperature { }, Lambda { }, Shear { }, StepSize { };
	bool MPC_on {false};
	Thermostat *thermostat{};


	//defaults für: TypeA, TypeB, Lambda, Temperature, BoxSize(x, y, z), stepsize, step_aufwärm, step_total, step_output
	double a_para[]{5, 0, 1.0, 0.5, 10, 10, 10, 0.001, 1E3, 1E5, 1E3};
	int a_para_size = sizeof(a_para) / sizeof(*a_para);
	int i_para {1};
	for (i_para = 1; i_para < min(a_para_size + 1, argc); ++i_para) {
		if (is_number(argv[i_para])) a_para[i_para - 1] = stod(argv[i_para]);
		else break;
	}
	TypeA = (int)a_para[0];
	TypeB = (int)a_para[1];
	Lambda = a_para[2];
	Temperature = a_para[3];
	BoxX = a_para[4];
	BoxY = a_para[5];
	BoxZ = a_para[6];
	StepSize = a_para[7];
	Steps_Equil = (int)a_para[8];
	Steps_Total = (int)a_para[9];
	Steps_Output = (int)a_para[10];

	Box box(BoxX, BoxY, BoxZ, Temperature, Lambda);
	box.add_chain(TypeA, TypeB, 10., 1.01);

	while (i_para < argc - 1 && is_number(argv[i_para])) ++i_para;
	int i_MPC { i_para };
	++i_para;
	if (argc > 1 && i_MPC < argc) {
		if (strcmp(argv[i_MPC], "MPC") == 0) {
			MPC_on = true;
			thermostat = new Thermostat_None{ box, StepSize };
			Shear = set_param(0., argv, argc, i_MPC + 1);
		}
	}
	else thermostat = new Andersen{box, StepSize, Temperature, 1999};


	MPC MPCroutine{box, Temperature, Shear};

	std::cout << "Type A: " << TypeA << " Type B: " << TypeB << " Lambda: " << Lambda << std::endl;
	std::cout << "Temperature: " << Temperature <<  std::endl;
	std::cout << "Box Size: " << BoxX << " " << BoxY << " " << BoxZ << std::endl;
	std::cout << "Step Size: " << StepSize << std::endl;
	std::cout << "Total Steps: " << Steps_Total << "Equilibration: " << Steps_Equil << "Output every: " << Steps_Output << std::endl;
	if (MPC_on) {
		std::cout << "MPC is turned ON with shear rate: " << Shear << std::endl;
	}
	else std::cout << "MPC is turned OFF" << std::endl;

	ofstream temp_file;
	ofstream config_file;
	temp_file.open("temperature2.dat", ios::out | ios::trunc);
	//config_file.open("config3.dat", ios::out | ios::trunc);


	if (MPC_on) MPCroutine.initializeMPC();
	clock_t begin = clock();

	for (int n = 0; n < Steps_Total; n++) {

		thermostat -> propagate(false);
		if (n > Steps_Equil && !(n%Steps_Output)) {
			thermostat -> propagate(true);
			std::cout << n << " ";
			box.print_Epot(std::cout);
			box.print_Ekin(std::cout);
			if (MPC_on) std::cout << MPCroutine.calculateCurrentTemperature() << " ";
			else box.print_Temperature(std::cout);
			temp_file << n << " ";
			box.print_Epot(temp_file);
			box.print_Ekin(temp_file);
			box.print_Temperature(temp_file);
			box.print_radius_of_gyration(temp_file);
			temp_file << "\n";
			std::cout << '\n';
		}
		if (MPC_on && !(n%10)) {
			MPCroutine.MPCstep(0.01);
		}
	}

	clock_t end = clock();
	//box.print_molecules(std::cout);
	std::cout << "time: " << double(end-begin)/CLOCKS_PER_SEC << std::endl;
	delete thermostat;
}
