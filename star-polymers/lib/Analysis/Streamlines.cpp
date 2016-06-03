/*
 * Streamlines.cpp
 *
 *  Created on: Jun 3, 2016
 *      Author: formanek
 */
#include "Streamlines.h"

ESADStreamlines::ESADStreamlines(double d, double L_x, double L_y, unsigned n_directions) :
d_sep{d},
Lx{L_x},
Ly{L_y}
{
	CellSize(0) = (int)ceil(Lx/d_sep);
	CellSize(1) = (int)ceil(Ly/d_sep);
	CellList = std::vector<std::vector<std::forward_list<Vector2d*>>>(CellSize(0), std::vector<std::forward_list<Vector2d*>>(CellSize(1), std::forward_list<Vector2d*>()));
	Directions = std::vector<Vector2d>(n_directions, Vector2d());
	double d_phi = 2.*M_PI/n_directions;
	for (unsigned i = 0; i < n_directions; i++)
	{
		Directions[i](0) = d_sep*cos(d_phi*i);
		Directions[i](1) = d_sep*sin(d_phi*i);
	}

}

void ESADStreamlines::initialize(string filename) {
	std::ifstream filestream(filename);
	string line;
	while(getline(filestream, line)) {
		std::vector<Vector2d> row {};
		std::stringstream stream{line};
		Vector2d vel{};
		while(stream >> vel(0) >> vel(1)) {
			row.push_back(vel);
		}
		Velocities.push_back(row);
	}
}

double ESADStreamlines::distance(Vector2d first, Vector2d second) {
	return (first-second).norm();
}

bool ESADStreamlines::check_valid(Vector2d seedpoint) {
	int cell_x {(int)(seedpoint(0)/d_sep)};
	int cell_y {(int)(seedpoint(1)/d_sep)};
	for (int i = cell_x - 1; i < cell_x + 2; i++) {
		if (i < 0 || i >= CellSize(0)) continue;
		for (int j = cell_y - 1; j < cell_y + 2; j++) {
			if (j < 0 || j >= CellSize(1)) continue;
			for (auto& point : CellList[i][j]) {
				if (distance(seedpoint, point) < d_sep) return false;
			}
		}
	}
	return true;
}

void ESADStreamlines::make_streamline(Vector2d seedpoint) {

}


