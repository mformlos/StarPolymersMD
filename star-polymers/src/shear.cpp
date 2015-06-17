/*
 * shear.cpp
 *
 *  Created on: Jun 16, 2015
 *      Author: maud
 */


#include <iostream>
#include <fstream>
#include <forward_list>
#include <ctime>
#include <cmath>
#include <algorithm>
#include <stdio.h>
#include "Box.h"
#include "Thermostat_None.h"
#include "Lowe_Andersen.h"
#include "Nose_Hoover.h"
#include "Andersen.h"
#include "Hydrodynamics_None.h"
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

	int BoxX { }, BoxY { }, BoxZ { }, TypeA { }, TypeB { }, Arms { }, Steps_Equil { }, Steps_Total { }, Steps_Output { };
	double Temperature { }, Lambda { }, Shear { }, StepSize { };

	stringstream ss_para { };



	//defaults für: TypeA, TypeB, Arms, Lambda, Temperature, BoxSize(x, y, z), stepsize, step_aufwärm, step_total, step_output
	double a_para[]{3, 0, 3, 1.0, 0.5, 10, 10, 50, 0.001, 1E2, 1E4, 1E2, 0.0};
	int a_para_size = sizeof(a_para) / sizeof(*a_para);
	int i_para { }, start_i_para { };
	if (argc > 1) {
		if (is_number(argv[1])) start_i_para = 1;
		else start_i_para = 2;
	}
	for (i_para = start_i_para; i_para < min(a_para_size + 2, argc); ++i_para) {
		if (is_number(argv[i_para])) {
			if ( i_para == 4 && strcmp(argv[1], "Chain") == 0) a_para[i_para- start_i_para +1] =stod(argv[i_para]);
			else a_para[i_para - start_i_para] = stod(argv[i_para]);
		}
		else break;
	}
	TypeA = (int)a_para[0];
	TypeB = (int)a_para[1];
	if (argc > 1 && strcmp(argv[1], "Chain") == 0) Arms = 0;
	else Arms = (int)a_para[2];
	Lambda = a_para[3];
	Temperature = a_para[4];
	BoxX = a_para[5];
	BoxY = a_para[6];
	BoxZ = a_para[7];
	StepSize = a_para[8];
	Steps_Equil = (int)a_para[9];
	Steps_Total = (int)a_para[10];
	Steps_Output = (int)a_para[11];
	Shear = a_para[12];

	Box box(BoxX, BoxY, BoxZ, Temperature, Lambda);

	while (i_para < argc - 1 && is_number(argv[i_para])) ++i_para;


	Andersen andersen{box, StepSize, Temperature, 1999};


	std::cout << "Type A: " << TypeA << " Type B: " << TypeB << " Arms: " << Arms << " Lambda: " << Lambda << std::endl;
	std::cout << "Temperature: " << Temperature <<  std::endl;
	std::cout << "Box Size: " << BoxX << " " << BoxY << " " << BoxZ << std::endl;
	std::cout << "Step Size: " << StepSize << std::endl;
	std::cout << "Total Steps: " << Steps_Total << " Equilibration: " << Steps_Equil << " Output every: " << Steps_Output << std::endl;
	std::cout << "MPC is turned ON with shear rate: " << Shear << std::endl;


	ss_para.precision(0);
	if (argc > 1 && strcmp(argv[1], "Chain") == 0) ss_para << "_Chain";
	else ss_para << "_Star";
	ss_para << "_A" << TypeA;
	ss_para << "_B" << TypeB;
	if (!(argc > 1 && strcmp(argv[1], "Chain") == 0)) ss_para << "_Arms" << Arms;
	ss_para << "_Lx" << BoxX;
	ss_para << "_Ly" << BoxY;
	ss_para << "_Lz" << BoxZ;
	ss_para << "_Lambda" << Lambda;
	ss_para << "_T" << Temperature;
	ss_para << "_run" << scientific << Steps_Total;
	ss_para << "_MPCON_Shear" << Shear;



	ofstream statistic_file { };
	//ofstream config_file { };
	//FILE* config_file { };

	string statistic_file_name = "../results/statistics"+ss_para.str()+".dat";
	string config_file_name = "../results/config"+ss_para.str()+".dat";
	statistic_file.open(statistic_file_name, ios::out | ios::trunc);
	//config_file = fopen(config_file_name.c_str(), "w");
	if (argc > 1 && strcmp(argv[1], "Chain") == 0) {
		box.add_chain(TypeA, TypeB, 10.);
		std::cout << "building a chain" << std::endl;
	}
	else	 {
		box.	add_star(TypeA, TypeB, Arms, 10.);
		std::cout << "building a star" << std::endl;
	}

	//box.print_molecules(std::cout);


	clock_t begin = clock();
	int n {  };
	for (n = 0; n < Steps_Equil; ++n) {

		andersen.propagate(false);

	}

	Thermostat_None integrator{box, StepSize};
	MPC mpc{box, Temperature, Shear};
	mpc.initialize();
	for (n = 0; n < Steps_Total;n++) {

		if (!(n%Steps_Output)) {
			integrator.propagate(true);
			std::list<unsigned> patches = box.calculate_patches();
			int number_of_patches {0};
			double av_patch_size {0.0};
			for (auto& element : patches) {
				if (element > 1) {
					number_of_patches++;
					av_patch_size += element;
				}
			}
			if (number_of_patches > 0) av_patch_size /= number_of_patches;
			Matrix3d gyr_tensor = box.calculate_gyration_tensor();
			//std::cout << '\n';
			//box.print_PDB(config_file, n);
			statistic_file << n << " ";
			box.print_Epot(statistic_file);
			box.print_Ekin(statistic_file);
			box.print_Temperature(statistic_file);
			box.print_radius_of_gyration(statistic_file);
			statistic_file << number_of_patches << " " << av_patch_size;
			for (int i = 0; i < gyr_tensor.size(); i++) {
				statistic_file << *(gyr_tensor.data()+i) << " ";
			}
			statistic_file << "\n";
			//std::cout << '\n';
		}
		else integrator.propagate(false);
		if (!(n%100)) {
			mpc.step(0.1);
		}
	}

	clock_t end = clock();
	//box.print_molecules(std::cout);
	std::cout << "time: " << double(end-begin)/CLOCKS_PER_SEC << std::endl;

}


