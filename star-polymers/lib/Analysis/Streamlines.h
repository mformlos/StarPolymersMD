/*
 * Streamlines.h
 *
 *  Created on: Jun 3, 2016
 *      Author: formanek
 */

#ifndef STREAMLINES_H_
#define STREAMLINES_H_

#include <iostream>
#include <array>
#include <vector>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <sstream>
#include <list>
#include <forward_list>
#include <../eigen/Eigen/Dense>
using namespace Eigen;


class ESADStreamlines {
public:
	struct Streamline {
		std::list<Vector2d> points;
	};
	double d_sep;
	double d_test;
	double Lx;
	double Ly;
	std::vector<std::vector<Vector2d>> Velocities;
	std::list<Streamline> Streamline_queue;
	std::vector<std::vector<std::forward_list<Vector2d>>> CellList;
	std::vector<Vector2d> Directions;
	Vector2i CellSize;


	ESADStreamlines(double d, double L_x, double L_y, unsigned n_directions);

	void initialize(std::string filename);
	Vector2d velocity(Vector2d position);
	bool check_valid_seed(Vector2d& seedpoint);
	bool check_valid_line(Vector2d position);
	void make_streamline(Vector2d seedpoint);

	void iterate();

	void print_vel();
	void print_streamline(Streamline& line);
	void print_CellList();
};

#endif /* STREAMLINES_H_ */
