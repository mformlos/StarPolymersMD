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
#include "Andersen_stochastic.h"
#include "Hydrodynamics_None.h"
#include "MPC.h"
#include "Velocity_x.h"
#include "Vel_Autocorr.h"


int signal_caught { };
void signal_handler( int signum )
{
    signal_caught = signum;
    std::cout << "signal " << signum << " caught" << std::endl;
    // cleanup and close up stuff here
    // terminate program
}

void file_handling(ofstream& os, string ss_para, string type, string old_steps, string new_steps) {
	string oldname { };
	string newname { };
	oldname = "./results/"+ type+ss_para;
	newname = oldname;
	newname.replace(newname.find(old_steps), old_steps.length(),new_steps);
	rename(oldname.c_str(), newname.c_str());
	os.close();
}
void file_handling(FILE* file, string ss_para, string type, string old_steps, string new_steps) {
	string oldname { };
	string newname { };
	oldname = "./results/"+ type+ss_para;
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
	int BoxX { }, BoxY { }, BoxZ { }, N{ }, Mass_mono{ };
	long int Steps_Equil { }, Steps_Total { }, Steps_Output { }, Steps_Start { }; 
	long int Steps_pdb { }, Steps_fluid { };
	double Temperature { }, Shear { }, StepSize { }, MPC_Step{ };
	bool MPC_on {false}, continue_run {false}, pdb_print {false}, fluid_print {false}, overwrite{true}, set_zero{false}, fluid_profile_print {false};
	string s_para { };
	stringstream ss_para { }, ss_para_old { };
	Thermostat *thermostat{};
	//Hydrodynamics *hydrodynamics{};


	//defaults für:  N, Mass_mono, Temperature, BoxSize(x, y, z), stepsize, step_aufwärm, step_total, step_output
	long double a_para[]{20, 1, 1, 100, 100, 100, 0.001, 1E2, 1E5, 1E2};
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
		N = (int)a_para[0];
		Mass_mono = a_para[1];
		Temperature = a_para[2];
		BoxX = a_para[3];
		BoxY = a_para[4];
		BoxZ = a_para[5];
		StepSize = a_para[6];
		Steps_Equil = (long int)a_para[6];
		Steps_Total = (long int)a_para[8];
		Steps_Output = (long int)a_para[9];
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
				else if (strcmp(argv[i_para], "profile") == 0) {
					fluid_profile_print = true;
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
		Mass_mono = stoi(find_parameter(s_para,"Mass"));
		Steps_Start = stol(find_parameter(s_para,"run"));
		std::cout << "Steps Start: " << Steps_Start << std::endl;
		if (Steps_Start == 0) set_zero = true;
		Temperature = stod(find_parameter(s_para,"T"));
		Shear = stod(find_parameter(s_para, "Shear"));

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
			else if (strcmp(argv[i_para], "profile") == 0) {
				fluid_profile_print = true;
			}
			i_para++;
		}
	}


	ss_para.precision(0);
	ss_para << "_Gaussian";
	ss_para << "_N" << N;
	ss_para << "_Mass" << Mass_mono;
	ss_para << "_Lx" << BoxX;
	ss_para << "_Ly" << BoxY;
	ss_para << "_Lz" << BoxZ;
	ss_para.precision(2);
	ss_para << fixed << "_T" << Temperature;
	ss_para << "_run" << scientific << Steps_Total;
	if (MPC_on) {
		ss_para.precision(4);
		ss_para << fixed << "_MPCON_Shear" << Shear;
	}
	else ss_para << "_MPCOFF";

	ofstream statistic_file { };
	FILE* config_file { };
	FILE* fluid_file { };
	FILE* end_fluid_file { };
	ofstream fluid_profile { };
	ofstream autocorrelation_file { };
	string end_fluid_file_name {"./results/end_fluid"+ss_para.str()+".dat"};
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
	if (add_stat_descriptor) statistic_file << "Step   Epot     Ekin     Temp   n_p s_p  R_gyr   G_xx      G_xy      G_xz       G_yx       G_yy      G_yz      G_zx     Gzy     Gzz      wx      wy     wz\n";

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
	if(MPC_on && fluid_profile_print) {
		oldname = "./results/fluid_profile"+ss_para_old.str()+".dat";
		newname = "./results/fluid_profile"+ss_para.str()+".dat";
		if (continue_run && overwrite) rename(oldname.c_str(), newname.c_str());
		string fluid_profile_name = newname;
		fluid_profile.open(fluid_profile_name, ios::out | ios::trunc);
	}

	if(MPC_on) {
		oldname = "./results/vel_autocorr"+ss_para_old.str()+".dat";
		newname = "./results/vel_autocorr"+ss_para.str()+".dat";
		if (continue_run && overwrite) rename(oldname.c_str(), newname.c_str());
		string autocorrelation_name = newname;
		autocorrelation_file.open(autocorrelation_name, ios::out | ios::trunc);
	}

	FILE* end_config_file { };
	string end_config_file_name = "./results/end_config"+ss_para.str()+".pdb";
     
	ofstream output_file { };
	string output_file_name = "./results/output"+ss_para.str()+".dat";
	output_file.open(output_file_name, ios::out | ios::trunc);
	output_file << "Monomers: " << N <<  " Mass: " << Mass_mono << std::endl;
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


	MPC_Step = 0.1;
	Box box(BoxX, BoxY, BoxZ, Temperature, 1.);

	MPC hydrodynamics{box, Temperature, Mass_mono, Shear, 100, false};
	VelocityX velocity_average_x{box, hydrodynamics,0.2};
	Vel_Autocorr autocorrelation{100,MPC_Step};

	if (continue_run) {
			box.add_gaussian(std::string(argv[1]), N, double(Mass_mono), set_zero);
			box.print_molecules(output_file);
	}
	else {
		box.add_gaussian(N, double(Mass_mono));
		output_file << "building a gaussian chain" << std::endl;
		box.print_molecules(std::cout);
	}

	std::cout << '\n';
	box.print_center_of_mass(std::cout);


    if (continue_run && MPC_on && fluid_profile_print) {
      	velocity_average_x.initialize("./results/fluid_profile"+ss_para_old.str()+".dat");
      	std::cout << "initialized" << std::endl;
    }

	if (MPC_on) {
		thermostat = new Thermostat_None{ box, StepSize, true };
		//hydrodynamics = new MPC{box, Temperature, Shear};
    }
	else {
		//thermostat = new Andersen{box, StepSize, Temperature, 1999};
		thermostat = new AndersenStochastic{box, StepSize, Temperature, 0.1, true};
		//thermostat = new Lowe_Andersen{box, StepSize, Temperature, 20., 1.5};
		//hydrodynamics = new Hydrodynamics_None{box};
    }


	if (MPC_on) {
		string continue_fluid {"./results/end_fluid"+ss_para_old.str()+".dat"};
		if (file_exists(continue_fluid)) {
			hydrodynamics.initialize(continue_fluid);
			std::cout << "yeah" << std::endl;
		}
		else hydrodynamics.initialize(); //hydrodynamics -> initialize();
	}
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
			file_handling(statistic_file, ss_para.str()+".dat", "statistics", Step_total, new_total_str);
			file_handling(output_file, ss_para.str()+".dat", "output", Step_total, new_total_str);
            if (MPC_on && fluid_profile_print) {
            	std::cout << "printing fluid profile..." << std::endl;
    		    velocity_average_x.print_result(fluid_profile);
            	file_handling(fluid_profile, ss_para.str()+".dat", "fluid_profile", Step_total, new_total_str);
            }
            if (pdb_print) {
            	std::cout << "printing pdb..." << std::endl;
            	file_handling(config_file, ss_para.str()+".pdb", "config", Step_total, new_total_str);
            }
            if (fluid_print) file_handling(fluid_file, ss_para.str()+".dat", "fluid", Step_total, new_total_str);
            file_handling(autocorrelation_file, ss_para.str()+".dat", "vel_autocorr", Step_total, new_total_str);
            autocorrelation.print_result(autocorrelation_file);

			newname = "./results/end_config"+ss_para.str()+".pdb";
			newname.replace(newname.find(Step_total), Step_total.length(),new_total.str());
			end_config_file = fopen(newname.c_str(), "w"); 
			box.print_PDB_with_velocity(end_config_file, n);
			fclose(end_config_file);
			if (MPC_on) {
				newname =  "./results/end_fluid"+ss_para.str()+".dat";
				newname.replace(newname.find(Step_total), Step_total.length(),new_total.str());
				end_fluid_file = fopen(newname.c_str(), "w");
				hydrodynamics.print_fluid_complete(end_fluid_file);
				fclose(end_fluid_file);
				string overwrite_fluid {"./results/end_fluid"+ss_para_old.str()+".dat"};
				if (continue_run && overwrite) remove(overwrite_fluid.c_str());
			}


			output_file << "signal " << signal_caught << "caught, data saved" << std::endl;

            if(continue_run && overwrite) remove(argv[1]);
			exit(signal_caught);
		}

		if (n > Steps_Equil && !(n%Steps_Output)) {
			thermostat -> propagate_gaussian(true);
			box.center_of_mass_reference_frame();
			std::tuple<double,Matrix3d> gyration = box.calculate_gyration_tensor();
			statistic_file << n << " ";
			box.print_Epot(statistic_file);
			box.print_Ekin_and_Temperature(statistic_file);
			statistic_file << std::get<0>(gyration) << "\n";

			if(MPC_on && fluid_profile_print) hydrodynamics(velocity_average_x);

			std::cout << n << " ";
			box.print_Epot(std::cout);
			output_file << n << " ";
			box.print_center_of_mass(output_file);
			output_file << '\n';
			output_file.flush();
			box.print_Temperature(std::cout);
			std::cout << " " << std::get<0>(gyration) << '\n';

			if (pdb_print && (Steps_pdb > 0 ? !(n%Steps_pdb) : true)) {
				box.print_PDB_with_velocity(config_file, n);
				fflush(config_file);
			}
			if (fluid_print && (Steps_fluid > 0 ? !(n%Steps_fluid) : true)) {
				hydrodynamics.print_fluid(fluid_file, n, (int)BoxZ/2, (int)BoxZ/2);
				//hydrodynamics.print_fluid_with_coordinates(fluid_file, n, (int)BoxZ/2 - 2, (int)BoxZ/2 + 2);
				//hydrodynamics -> print_fluid_with_coordinates(fluid_file, n, (int)BoxZ/2 - 2, (int)BoxZ/2 + 2);
				fflush(fluid_file);
			}

			//std::cout << '\n';

		}
		else thermostat -> propagate_gaussian(false);

		if (MPC_on) {
			hydrodynamics.step(n, MPC_Step); //hydrodynamics -> step(1.0);
			if (n > Steps_Equil) box(autocorrelation);
		}
	}

    if (MPC_on && fluid_profile_print) velocity_average_x.print_result(fluid_profile);
	end_config_file = fopen(end_config_file_name.c_str(), "w");
	box.print_PDB_with_velocity(end_config_file, Steps_Total);
    if(continue_run && overwrite) remove(argv[1]);

    if (MPC_on) {
    	autocorrelation.print_result(autocorrelation_file);
    	end_fluid_file = fopen(end_fluid_file_name.c_str(), "w");
    	hydrodynamics.print_fluid_complete(end_fluid_file);
    	fclose(end_fluid_file);
    	string overwrite_fluid {"./results/end_fluid"+ss_para_old.str()+".dat"};
    	if (continue_run && overwrite) remove(overwrite_fluid.c_str());
    }

	clock_t end = clock();
	//box.print_molecules(std::cout);
	output_file << "time: " << double(end-begin)/CLOCKS_PER_SEC << std::endl;

}
