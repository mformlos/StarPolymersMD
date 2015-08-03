/*
 * continue_star_polymers.cpp
 *
 *  Created on: Jul 9, 2015
 *      Author: maud
 */

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

string find_parameter(string para_string, string para_name) {
	string para {"0.0"};
	std::string::size_type i = para_string.find(para_name);
	if (i == std::string::npos) return para;
	para = para_string.substr(i + para_name.length());
	i = para.find("_");
	if (i == std::string::npos) i = para.find(".");
	if (i != std::string::npos) para.erase(i);

	return para;
}

int main(int argc, char* argv[]) {

	int BoxX { }, BoxY { }, BoxZ { }, TypeA { }, TypeB { }, Arms { };
	long int Steps_Start { }, Steps_Total { }, Steps_Output { };
	double Temperature { }, Lambda { }, Shear { }, StepSize { };
	bool MPC_on {false};
	string s_para { };
	stringstream ss_para { }, ss_para_new { };
	Thermostat *thermostat{};
	Hydrodynamics *hydrodynamics{};

	s_para = std::string(argv[1]);


	std::string find_this = "end_config";
	std::string::size_type i = s_para.find(find_this);
	std::cout << i << std::endl;
	if (i != std::string::npos){
		s_para.erase(i, find_this.length());
		s_para.erase(0, i);
	}


	find_this = ".pdb";
	i = s_para.find(find_this);
	if (i != std::string::npos) s_para.erase(i, find_this.length());
	ss_para << s_para;

	std::cout << s_para << std::endl;
	std::cout << find_parameter(s_para, "Shear") << std::endl;
	BoxX = stoi(find_parameter(s_para,"Lx"));
	BoxY = stoi(find_parameter(s_para,"Ly"));
	BoxZ = stoi(find_parameter(s_para,"Lz"));
	TypeA = stoi(find_parameter(s_para,"A"));
	TypeB = stoi(find_parameter(s_para,"B"));
	Arms = stoi(find_parameter(s_para,"Arms"));
	Steps_Start = stol(find_parameter(s_para,"run"));
	Temperature = stod(find_parameter(s_para,"T"));
	Shear = stod(find_parameter(s_para, "Shear"));
	Lambda = stod(find_parameter(s_para,"Lambda"));

	StepSize = stod(argv[2]);
	Steps_Total = Steps_Start + stold(argv[3]);
	Steps_Output = stold(argv[4]);
	if (find_parameter(s_para, "MPC").compare("ON") == 0) {
		MPC_on = true;
	}

	std::cout << "blubb" << std::endl;
	int arg { };
	while (arg < argc) {
		if (strcmp(argv[arg], "MPC") == 0) {
			MPC_on = true;
			std::cout << arg << std::endl;
			if (arg +1 < argc) Shear = stod(argv[arg+1]);
			arg++;
			break;
		}
		arg++;
	}


	if (argc >= 8) {
		BoxX = stod(argv[5]);
		BoxY = stod(argv[6]);
		BoxZ = stod(argv[7]);
	}



	std::cout << "Type A: " << TypeA << " Type B: " << TypeB << " Arms: " << Arms << " Lambda: " << Lambda << std::endl;
	std::cout << "Temperature: " << Temperature <<  std::endl;
	std::cout << "Box Size: " << BoxX << " " << BoxY << " " << BoxZ << std::endl;
	std::cout << "Step Size: " << StepSize << std::endl;
	std::cout << "Start at: " << Steps_Start << " Stop at: " << Steps_Total << " Output every: " << Steps_Output << std::endl;
	if (MPC_on) {
		std::cout << "MPC is turned ON with shear rate: " << Shear << std::endl;
	}
	else std::cout << "MPC is turned OFF" << std::endl;



	ss_para_new.str(std::string());
	ss_para_new.precision(0);
	ss_para_new << "_Star";
	ss_para_new << "_A" << TypeA;
	ss_para_new << "_B" << TypeB;
	ss_para_new << "_Arms" << Arms;
	ss_para_new << "_Lx" << BoxX;
	ss_para_new << "_Ly" << BoxY;
	ss_para_new << "_Lz" << BoxZ;
	ss_para_new.precision(3);
	ss_para_new << "_Lambda" << Lambda;
	ss_para_new << "_T" << Temperature;
	ss_para_new << "_run" << scientific << Steps_Total;
	if (MPC_on) ss_para_new << "_MPCON_Shear" << Shear;
	else ss_para_new << "_MPCOFF";

	ofstream statistic_file { };
	FILE* config_file { };
	string oldname = "./results/statistics"+ss_para.str()+".dat";
	string newname = "./results/statistics"+ss_para_new.str()+".dat";
	std::cout << oldname << std::endl;
	std::cout << newname << std::endl;
	rename(oldname.c_str(), newname.c_str());
	string statistic_file_name = newname;
	oldname = "./results/config"+ss_para.str()+".pdb";
	newname = "./results/config"+ss_para_new.str()+".pdb";
	rename(oldname.c_str(), newname.c_str());
	string config_file_name = newname;
	statistic_file.open(statistic_file_name, ios::out | ios::app);
	config_file = fopen(config_file_name.c_str(), "a");

	Box box(BoxX, BoxY, BoxZ, Temperature, Lambda);


	if (argc > 1 && strcmp(argv[1], "Chain") == 0) {
		box.add_chain(TypeA, TypeB, 10.);
		std::cout << "building a chain" << std::endl;
	}
	else	 {
		box.add_star(std::string(argv[1]), TypeA, TypeB, Arms, 10.);
		std::cout << "building a star" << std::endl;
	}



	if (MPC_on) {
		thermostat = new Thermostat_None{ box, StepSize };
		hydrodynamics = new MPC{box, Temperature, Shear};

	}
	else {
		thermostat = new Andersen{box, StepSize, Temperature, 1999};
		hydrodynamics = new Hydrodynamics_None{box};
	}



	if (MPC_on) hydrodynamics -> initialize();
	clock_t begin = clock();

	std::cout << thermostat -> info() << std::endl;
	box.print_molecules(std::cout);

	for (long int n = Steps_Start + 1; n <= Steps_Total; ++n) {


		if (!(n%Steps_Output)) {
			thermostat -> propagate(true);
			std::cout << n << " ";
			box.print_Epot(std::cout);
			box.print_Ekin(std::cout);
			if (MPC_on) std::cout << hydrodynamics -> calculateCurrentTemperature() << " ";
			else box.print_Temperature(std::cout);
			std::list<unsigned> patches = box.calculate_patches();
			int number_of_patches {0};
			double av_patch_size {0.0};
			for (auto& element : patches) {
				if (element > 1) {
					number_of_patches++;
					av_patch_size += element;
				}
				//std::cout << element << ' ';
			}
			if (number_of_patches > 0) av_patch_size /= number_of_patches;
			Matrix3d gyr_tensor = box.calculate_gyration_tensor();
			std::cout << '\n';
			box.print_PDB(config_file, n);
			statistic_file << n << " ";
			box.print_Epot(statistic_file);
			box.print_Ekin(statistic_file);
			box.print_Temperature(statistic_file);
			box.print_radius_of_gyration(statistic_file);
			statistic_file << number_of_patches << " " << av_patch_size << " ";
			for (int i = 0; i < gyr_tensor.size(); i++) {
				statistic_file << *(gyr_tensor.data()+i) << " ";
			}
			statistic_file << "\n";
			//std::cout << '\n';
		}
		else thermostat -> propagate(false);
		if (MPC_on && !(n%100)) {
			hydrodynamics -> step(0.1);
		}
	}

	FILE* end_config_file { };
	string end_config_file_name = "./results/end_config"+ss_para_new.str()+".pdb";
	end_config_file = fopen(end_config_file_name.c_str(), "w");
	box.print_PDB(end_config_file, Steps_Total);


	clock_t end = clock();
	//box.print_molecules(std::cout);
	std::cout << "time: " << double(end-begin)/CLOCKS_PER_SEC << std::endl;

}



