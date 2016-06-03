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

struct Point {
	double x;
	double y;
}
struct Streamline {
	std::list<Point> points;
}

int main(int argc, char* argv[]) {
	double d_sep{}, Lx{}, Ly{}, d_phi{};
	int cellsize_x{}, cellsize_y{};
	int n_directions{}, n_cells{};
	bool finished{false};
	std::list<Streamline> streamline_queue {};
	std::vector<std::vector<Point*>> CellList{};
	std::vector<Points> directions {};

	d_sep = (double)argv[1];
	n_directions = (int)argv[2];

	d_phi = 2.*M_PI/n_directions;
	for (unsigned i = 0; i < n_directions; i++) {
		directions.push_back(Point{d_sep*cos(d_phi*i), d_sep*sin(d_phi*i)});
	}

	cellsize_x = (int)(Lx/d_sep);
	cellsize_y = (int)(Ly/d_sep);
	n_cells = cellsize_x*cellsize_y;
	CellList = std::vector<std::vector<Point*>>(n_cells+1, std::vector<Point*>());



	Point seedpoint {0.,0.};
	streamline_queue.push_back(calculate_streamline(seedpoint));
	Streamline* current = streamline_queue.front();
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
	}
}
