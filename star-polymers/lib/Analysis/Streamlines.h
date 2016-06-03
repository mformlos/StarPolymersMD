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
	std::vector<std::vector<std::forward_list<Vector2d*>>> CellList;
	std::vector<Vector2d> Directions;
	Vector2i CellSize;


	ESADStreamlines(double d, unsigned n_directions);

	void initialize(string filename);
	double distance(Vector2d first, Vector2d second);
	Streamline calculate_streamline(Vector2d seedpoint);
	bool check_valid(Vector2d seedpoint);
	void make_streamline(Vector2d seedpoint);


};

#endif /* STREAMLINES_H_ */
