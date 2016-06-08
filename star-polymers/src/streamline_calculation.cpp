/*
 * Streamlines.cpp
 *
 *  Created on: Jun 3, 2016
 *      Author: formanek
 */

#include <iostream>
#include <fstream>
#include <forward_list>
#include <ctime>
#include <cmath>
#include <algorithm>
#include <stdio.h>
#include <csignal>
#include <tuple>
#include "Streamlines.h"
#include <../eigen/Eigen/Dense>


int main(int argc, char* argv[]) {
	double d_sep{}, Lx{}, Ly{};
	int n_directions{};
	//bool finished{false};
	std::string filename{};


	d_sep = std::stod(argv[1]);
	Lx = std::stod(argv[2]);
	Ly = std::stod(argv[3]);
	n_directions = std::stoi(argv[4]);
	filename = argv[5];


    ESADStreamlines StreamlineRoutine(d_sep, Lx, Ly, n_directions);
    StreamlineRoutine.initialize(filename);
    StreamlineRoutine.print_vel();


	Vector2d seedpoint {3,2};
	StreamlineRoutine.make_streamline(seedpoint);
	//StreamlineRoutine.print_CellList();
	StreamlineRoutine.iterate();
	/*streamline_queue.push_back(calculate_streamline(seedpoint));
	//Streamline* current = streamline_queue.front();
	std::list<Streamline>::iterator current;
	current = streamline_queue.begin();
	while (!finished) {
		for(auto& p : current) {
			for (auto& dir : directions) {
				seedpoint.x = p.x + dir.x;
				seedpoint.y = p.y + dir.y;
				bool valid = check_valid(seedpoint);
				if (valid) streamline_queue.push_back(calculate_streamline(seedpoint));
			}
		}
		if (current == streamline_queue.end()) finished = true;
		else current++;
	}*/
    return 0;
}
