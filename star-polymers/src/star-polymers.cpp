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
#include <csignal>
#include <tuple>
#include "Box.h"
#include "Thermostat_None.h"
#include "Lowe_Andersen.h"
#include "Nose_Hoover.h"
#include "Andersen.h"
#include "Hydrodynamics_None.h"
#include "MPC.h"
#include "Velocity_x.h"


int signal_caught { };
void signal_handler( int signum )
{
    signal_caught = signum;
    std::cout << "signal " << signum << " caught" << std::endl;
    // cleanup and close up stuff here
    // terminate program
}

void file_handling(ofstream& os, stringstream& ss_para, string type, string old_steps, string new_steps) {
	string oldname { };
	string newname { };
	oldname = "./results/"+ type+ss_para.str()+".dat";
	newname = oldname;
	newname.replace(newname.find(old_steps), old_steps.length(),new_steps);
	rename(oldname.c_str(), newname.c_str());
	os.close();
}
void file_handling(FILE* file, stringstream& ss_para, string type, string old_steps, string new_steps) {
	string oldname { };
	string newname { };
	oldname = "./results/"+ type+ss_para.str()+".dat";
	newname = oldname;
	newname.replace(newname.find(old_steps), old_steps.length(),new_steps);
	rename(oldname.c_str(), newname.c_str());
	fclose(file);
}


bool file_is_empty(std::ifstream& pFile)
{
    return pFile.peek() == std::ifstream::traits_type::eof();
}

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
    //if (i == std::string::npos) i = para.find(".");
	if (i != std::string::npos) para.erase(i);
	return para;
}

inline bool file_exists (const std::string& name) {
    return ( access( name.c_str(), F_OK ) != -1 );
}

int main(int argc, char* argv[]) {
	signal_caught = 0;
	signal(SIGINT, signal_handler);
	signal(SIGSEGV, signal_handler);
	int BoxX { }, BoxY { }, BoxZ { }, TypeA { }, TypeB { }, Arms { };
	long int Steps_Equil { }, Steps_Total { }, Steps_Output { }, Steps_Start { }; 
	long int Steps_pdb { }, Steps_fluid { };
	double Temperature { }, Lambda { }, Shear { }, StepSize { }, MPC_Step{ };
	bool MPC_on {false}, continue_run {false}, pdb_print {false}, fluid_print {false}, overwrite{true};
	string s_para { };
	stringstream ss_para { }, ss_para_old { };
	Thermostat *thermostat{};
	//Hydrodynamics *hydrodynamics{};


	//defaults für: TypeA, TypeB, Arms, Lambda, Temperature, BoxSize(x, y, z), stepsize, step_aufwärm, step_total, step_output
	long double a_para[]{3, 0, 3, 1.0, 0.5, 10, 10, 50, 0.01, 1E2, 1E10, 1E2};
	int a_para_size = sizeof(a_para) / sizeof(*a_para);
	int i_para { }, start_i_para { };

	
	std::string find_this = "end_config"; 
	std::string::size_type i { }; 
	if (argc > 1) {
		if (is_number(argv[1])) {
			start_i_para = 1;
			continue_run = false; 
		}
		else {
			start_i_para = 2;
			s_para = std::string(argv[1]);
			i = s_para.find(find_this); 
			if (!file_exists(s_para)) {
				std::cout << "input file does not exist";
				return 1;
			}
			if (i!= std::string::npos) {
				continue_run = true; 
				s_para.erase(i, find_this.length());
				s_para.erase(0, i);
			}
		}
	}
	for (i_para = start_i_para; i_para < min(a_para_size + 2, argc); ++i_para) {
		if (is_number(argv[i_para])) {
			if ( i_para == 4 && strcmp(argv[1], "Chain") == 0) a_para[i_para- start_i_para +1] =stold(argv[i_para]);
			else a_para[i_para - start_i_para] = stold(argv[i_para]);
		}
		else break;
	}
        if (!continue_run) {
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
		Steps_Equil = (long int)a_para[9];
		Steps_Total = (long int)a_para[10];
		Steps_Output = (long int)a_para[11];
		while (i_para < argc) {
			while (i_para < argc - 1 && is_number(argv[i_para])) ++i_para;
			if (argc > 1 && i_para < argc) {
				if (strcmp(argv[i_para], "MPC") == 0) {
					MPC_on = true;
					Shear = set_param(0., argv, argc, i_para + 1);
				}
				else if (strcmp(argv[i_para], "pdb") == 0) {
					pdb_print = true;
					if (i_para + 1 < argc && is_number(argv[i_para + 1])){
						i_para++;
						Steps_pdb = (long int)stold(argv[i_para]);
						std::cout << Steps_pdb <<std::endl;

					}
				}
				else if (strcmp(argv[i_para], "fluid") == 0) {
					fluid_print = true;
					if (i_para + 1 < argc && is_number(argv[i_para + 1])){
						i_para++;
						Steps_fluid = (long int)stold(argv[i_para]);
						std::cout << Steps_fluid <<std::endl;
					}
				}
			}
			i_para++;
		}
	}
	else {
		find_this = "./results/";
		i = s_para.find(find_this);
		if (i != std::string::npos){
			s_para.erase(i, find_this.length());
		}

		find_this = ".pdb";
		i = s_para.find(find_this);
		if (i != std::string::npos) s_para.erase(i, find_this.length());
		ss_para_old << s_para;

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

		StepSize = a_para[0];
		Steps_Total = (long int)a_para[1];
		Steps_Output = (long int)a_para[2];
		Steps_Equil = 0;

		if (i_para > 5) {
			BoxX = a_para[3];
			if (i_para > 6) {
				BoxY = a_para[4];
			 	if (i_para > 7) BoxZ = a_para[5]; 
			}
		}

		if (find_parameter(s_para, "MPC").compare("ON") == 0) {
			MPC_on = true;
		}
		while (i_para < argc) {
			if (strcmp(argv[i_para], "MPC") == 0) {
				if (!MPC_on) {
					Steps_Start = 0;
					overwrite = false;
				}
				MPC_on = true;
				if (i_para + 1 < argc && is_number(argv[i_para + 1])){
					i_para++;
					if (Shear != stod(argv[i_para])) {
						Steps_Start = 0;
						overwrite = false;
					}
					Shear = stod(argv[i_para]);
				}
			}
			else if (strcmp(argv[i_para], "pdb") == 0) {
				pdb_print = true;
				if (i_para + 1 < argc && is_number(argv[i_para + 1])) {
					i_para++;
					Steps_pdb = (long int)stold(argv[i_para]);
					std::cout << Steps_pdb << " ";
				}
			}
			else if (strcmp(argv[i_para], "fluid") == 0) {
				fluid_print = true;
				if (i_para + 1 < argc && is_number(argv[i_para + 1])) {
					i_para++;
					Steps_fluid = (long int)stold(argv[i_para]);
					std::cout << Steps_fluid <<std::endl;
				}
			}
			i_para++;
		}
	}


	ss_para.precision(0);
	if (argc > 1 && strcmp(argv[1], "Chain") == 0) ss_para << "_Chain";
	else ss_para << "_Star";
	ss_para << "_A" << TypeA;
	ss_para << "_B" << TypeB;
	if (!(argc > 1 && strcmp(argv[1], "Chain") == 0)) ss_para << "_Arms" << Arms;
	ss_para << "_Lx" << BoxX;
	ss_para << "_Ly" << BoxY;
	ss_para << "_Lz" << BoxZ;
	ss_para.precision(2);
	ss_para << fixed << "_Lambda" << Lambda;
	ss_para << fixed << "_T" << Temperature;
	ss_para << "_run" << scientific << Steps_Total;
	if (MPC_on) {
		ss_para.precision(2);
		ss_para << fixed << "_MPCON_Shear" << Shear;
	}
	else ss_para << "_MPCOFF";

	ofstream statistic_file { };
	FILE* config_file { };
	FILE* fluid_file { };
	ofstream fluid_profile { };
	string oldname = "./results/statistics"+ss_para_old.str()+".dat";
	string newname = "./results/statistics"+ss_para.str()+".dat";
	if(continue_run && overwrite) rename(oldname.c_str(), newname.c_str());
	string statistic_file_name = newname;
	ifstream statistic_file_check {};
	statistic_file_check.open(statistic_file_name, ios::in);
    bool add_stat_descriptor {false};
	if (file_is_empty(statistic_file_check)) add_stat_descriptor = true;
	statistic_file_check.close();
	statistic_file.open(statistic_file_name, ios::out | ios::app);
	if (add_stat_descriptor) statistic_file << "Step   Epot     Ekin     Temp   n_p s_p  R_gyr   G_xx      G_xy      G_xz       G_yx       G_yy      G_yz      G_zx     Gzy     Gzz \n";

	if(pdb_print) {
		oldname = "./results/config"+ss_para_old.str()+".pdb";
		newname = "./results/config"+ss_para.str()+".pdb";
		if (continue_run && overwrite) rename(oldname.c_str(), newname.c_str());
		string config_file_name = newname;
		config_file = fopen(config_file_name.c_str(), "a");
	}

	if(fluid_print) {
		oldname = "./results/fluid"+ss_para_old.str()+".dat";
		newname = "./results/fluid"+ss_para.str()+".dat";
		if (continue_run && overwrite) rename(oldname.c_str(), newname.c_str());
		string fluid_file_name = newname;
		fluid_file = fopen(fluid_file_name.c_str(), "a");
	}
	if(MPC_on) {
		oldname = "./results/fluid_profile"+ss_para_old.str()+".dat";
		newname = "./results/fluid_profile"+ss_para.str()+".dat";
		if (continue_run && overwrite) rename(oldname.c_str(), newname.c_str());
		string fluid_profile_name = newname;
		fluid_profile.open(fluid_profile_name, ios::out | ios::trunc);
	}



	FILE* end_config_file { };
	string end_config_file_name = "./results/end_config"+ss_para.str()+".pdb";
     
	ofstream output_file { };
	string output_file_name = "./results/output"+ss_para.str()+".dat";
	output_file.open(output_file_name, ios::out | ios::trunc);
	output_file << "Type A: " << TypeA << " Type B: " << TypeB << " Arms: " << Arms << " Lambda: " << Lambda << std::endl;
	output_file << "Temperature: " << Temperature <<  std::endl;
	output_file << "Box Size: " << BoxX << " " << BoxY << " " << BoxZ << std::endl;
	output_file << "Step Size: " << StepSize << std::endl;
	output_file << "Start at: " << Steps_Start << " Stop at: " << Steps_Total << " Output every: " << Steps_Output << std::endl;
	if (MPC_on) {
		output_file << "MPC is turned ON with shear rate: " << Shear << std::endl;
    }
	else output_file << "MPC is turned OFF" << std::endl;
	output_file << "pdb file generated: " << (pdb_print ? "yes" : "no") << std::endl;
	output_file << "fluid file generated: " << (fluid_print ? "yes" : "no") << std::endl;


	MPC_Step = 100.0*StepSize;
	Box box(BoxX, BoxY, BoxZ, Temperature, Lambda);
    VelocityX velocity_average_x{0.2};

	if (argc > 1 && strcmp(argv[1], "Chain") == 0) {
		box.add_chain(TypeA, TypeB, 10.);
		output_file << "building a chain" << std::endl;
    }
	else {
		if (continue_run) box.add_star(std::string(argv[1]), TypeA, TypeB, Arms, 10.);
		else box.add_star(TypeA, TypeB, Arms, 10.);
		output_file << "building a star" << std::endl;
	}

	MPC hydrodynamics{box, Temperature, Shear};

	if (MPC_on) {
		thermostat = new Thermostat_None{ box, StepSize };
		//hydrodynamics = new MPC{box, Temperature, Shear};
    }
	else {
		thermostat = new Andersen{box, StepSize, Temperature, 1999};
		//thermostat = new Lowe_Andersen{box, StepSize, Temperature, 20., 1.5};
		//hydrodynamics = new Hydrodynamics_None{box};
    }



	if (MPC_on) hydrodynamics.initialize(); //hydrodynamics -> initialize();
	clock_t begin = clock();

	long int n { }; 
	if (continue_run) n = Steps_Start + 1;
	else n = 0; 
	for ( ; n <= Steps_Total; ++n) {
		if (signal_caught) {
			string Step_total = find_parameter(ss_para.str(),"run");
			stringstream new_total { };
			string new_total_str { };
			new_total.precision(0);
			new_total << scientific << n;
			new_total_str = new_total.str();
			file_handling(statistic_file, ss_para, "statistics", Step_total, new_total_str);
			file_handling(output_file, ss_para, "output", Step_total, new_total_str);
            if (MPC_on) file_handling(fluid_profile, ss_para, "fluid_profile", Step_total, new_total_str);
            if (pdb_print) file_handling(config_file, ss_para, "config", Step_total, new_total_str);
            if (fluid_print) file_handling(fluid_file, ss_para, "fluid", Step_total, new_total_str);


			newname = "./results/end_config"+ss_para.str()+".pdb";
			newname.replace(newname.find(Step_total), Step_total.length(),new_total.str());
			end_config_file = fopen(newname.c_str(), "w"); 
			box.print_PDB_with_velocity(end_config_file, n);
			fclose(end_config_file);

			output_file << "signal " << signal_caught << "caught, data saved" << std::endl;

            if(continue_run && overwrite) remove(argv[1]);
			exit(signal_caught);
		}

		if (n > Steps_Equil && !(n%Steps_Output)) {
			thermostat -> propagate(true);
			std::tuple<double,double> patches= box.calculate_patches_new();
			std::tuple<double,Matrix3d> gyration = box.calculate_gyration_tensor();
			statistic_file << n << " ";
			box.print_Epot(statistic_file);
			box.print_Ekin_and_Temperature(statistic_file);
			statistic_file << std::get<0>(patches) << " " << std::get<1>(patches) << " ";
			statistic_file << std::get<0>(gyration) << " ";
			Matrix3d gyr_tensor {std::get<1>(gyration)};
			for (int i = 0; i < gyr_tensor.size(); i++) {
				statistic_file << *(gyr_tensor.data()+i) << " ";
			}
			statistic_file << "\n";
			statistic_file.flush();
			if(MPC_on) hydrodynamics(velocity_average_x);

			/*std::cout << n << " ";
			box.print_Temperature(std::cout);
			std::cout << std::endl;
			output_file << n << " ";
			box.print_Epot(output_file);
			box.print_Ekin(output_file);
			if (MPC_on) output_file << hydrodynamics -> calculateCurrentTemperature() << " ";
			else box.print_Temperature(output_file);*/
			/*std::list<unsigned> patches = box.calculate_patches();
			int number_of_patches {0};
			double av_patch_size {0.0};
			for (auto& element : patches) {
				if (element > 1) {
					number_of_patches++;
					av_patch_size += element;
				}
				//std::cout << element << ' ';
			}
			if (number_of_patches > 0) av_patch_size /= number_of_patches;*/

			//output_file << std::get<0>(patches) << " " << std::get<1>(patches) << " ";

			//output_file << '\n';
			//output_file.flush();
			if (pdb_print && (Steps_pdb > 0 ? !(n%Steps_pdb) : true)) {
				box.print_PDB(config_file, n);
				fflush(config_file);
			}
			if (fluid_print && (Steps_fluid > 0 ? !(n%Steps_fluid) : true)) {
				hydrodynamics.print_fluid_with_coordinates(fluid_file, n, (int)BoxZ/2 - 2, (int)BoxZ/2 + 2);
				//hydrodynamics -> print_fluid_with_coordinates(fluid_file, n, (int)BoxZ/2 - 2, (int)BoxZ/2 + 2);
				fflush(fluid_file);
			}

			//std::cout << '\n';

		}
		else thermostat -> propagate(false);
		if (MPC_on && !(n%100)) {
			hydrodynamics.step(MPC_Step); //hydrodynamics -> step(1.0);
		}
	}

    velocity_average_x.print_result(fluid_profile);
	end_config_file = fopen(end_config_file_name.c_str(), "w");
	box.print_PDB_with_velocity(end_config_file, Steps_Total);
    if(continue_run && overwrite) remove(argv[1]);


	clock_t end = clock();
	//box.print_molecules(std::cout);
	output_file << "time: " << double(end-begin)/CLOCKS_PER_SEC << std::endl;

}
