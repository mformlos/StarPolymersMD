/*
 * Streamlines.cpp
 *
 *  Created on: Jun 3, 2016
 *      Author: formanek
 */
#include "Streamlines.h"

ESADStreamlines::ESADStreamlines(double d, double L_x, double L_y, unsigned n_directions) :
d_sep{d},
d_test{0.5*d},
Lx{L_x},
Ly{L_y}
{
	CellSize(0) = (int)ceil(Lx/d_sep);
	CellSize(1) = (int)ceil(Ly/d_sep);
	CellList = std::vector<std::vector<std::forward_list<Vector2d>>>(CellSize(0), std::vector<std::forward_list<Vector2d>>(CellSize(1), std::forward_list<Vector2d>()));
	Directions = std::vector<Vector2d>(n_directions, Vector2d());
	double d_phi = 2.*M_PI/n_directions;
	for (unsigned i = 0; i < n_directions; i++)
	{
		Directions[i](0) = d_sep*cos(d_phi*i);
		Directions[i](1) = d_sep*sin(d_phi*i);
	}

}

void ESADStreamlines::initialize(std::string filename) {
	std::ifstream filestream(filename);
	std::string line;
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



bool ESADStreamlines::check_valid_seed(Vector2d& seedpoint) {
	if (seedpoint(0) > Lx || seedpoint(0) < 0.0) return false;
	if (seedpoint(1) > Ly || seedpoint(1) < 0.0) return false;
	int cell_x {(int)(seedpoint(0)/d_sep)};
	int cell_y {(int)(seedpoint(1)/d_sep)};
	for (int i = cell_x - 1; i < cell_x + 2; i++) {
		if (i < 0 || i >= CellSize(0)) continue;
		for (int j = cell_y - 1; j < cell_y + 2; j++) {
			if (j < 0 || j >= CellSize(1)) continue;
			for (auto &point : CellList[i][j]) {
				if ((seedpoint - point).norm() < d_sep) return false;
			}
		}
	}
	return true;
}

bool ESADStreamlines::check_valid_line(Vector2d point) {
	std::cout << "check start" << std::endl;
    //print_CellList();
	if (point(0) > Lx || point(0) < 0.0) return false;
	if (point(1) > Ly || point(1) < 0.0) return false;
	int cell_x {(int)(point(0)/d_sep)};
	int cell_y {(int)(point(1)/d_sep)};
	for (int i = cell_x - 1; i < cell_x + 2; i++) {
		if (i < 0 || i >= CellSize(0)) continue;
		for (int j = cell_y - 1; j < cell_y + 2; j++) {
			if (j < 0 || j >= CellSize(1)) continue;
			for (auto& other : CellList[i][j]) {
				double distance {(point-other).norm()};
				std::cout << "point: " << point.transpose() << " other: " << other.transpose() << " distance: " << distance << std::endl;
				if (distance < d_test) return false;
			}
		}
	}
	std::cout << "check stop" << std::endl;

	return true;
}

Vector2d ESADStreamlines::velocity(Vector2d position) {
	int x1 {(int)position(0)};
	int x2 {x1+1};
	int y1 {(int)position(1)};
	int y2 {y1+1};
	Vector2d vel{};
	vel = ((double)x2-position(0))*(Velocities[x1][y1]*((double)y2-position(1)) + Velocities[x1][y2]*(position(1)-(double)y1));
	vel += (position(0)-(double)x1)*(Velocities[x2][y1]*((double)y2-position(1)) + Velocities[x2][y2]*(position(1)-(double)y1));
	return vel;
}


void ESADStreamlines::make_streamline(Vector2d seedpoint) {
	Streamline line{};
	bool finished{false};
	Vector2d x {};
	Vector2d v {};
    std::cout << "starting streamline" << std::endl;
    x = seedpoint;
    line.points.push_back(seedpoint);
	int cell_x {(int)(x(0)/d_sep)};
	int cell_y {(int)(x(1)/d_sep)};
	CellList[cell_x][cell_y].push_front(seedpoint);
    /*for (int i = 0; i < 10; i++) {
    	v = velocity(x);
    	x += velocity(x+d_sep*0.5*v);
		std::cout << v.transpose() << " " << x.transpose() << std::endl;

    }*/
	while (!finished) {
		v = velocity(x);
		Vector2d new_x {};
		v = velocity(x+d_sep*0.5*v);
		v /= v.norm();
		new_x = x + d_sep*v;
		x = new_x;
		std::cout << "x: " << new_x.transpose() << "v: " << velocity(new_x).transpose() << std::endl;
		if (v.squaredNorm() <= 0.000001) finished = true;
		if (check_valid_line(new_x)) {
			line.points.push_back(new_x);
			cell_x = (int)(new_x(0)/d_sep);
			cell_y = (int)(new_x(1)/d_sep);
			CellList[cell_x][cell_y].push_front(new_x);
		}
		else finished = true;
	}
	if (line.points.size() >= 2) Streamline_queue.push_back(line);
	std::cout << "finishing streamline" << std::endl;
	print_streamline(line);
}

void ESADStreamlines::iterate() {
	std::list<Streamline>::iterator current {Streamline_queue.begin()};
	unsigned counter {};
	bool done{false};
	Vector2d seedpoint{};
	while (!done) {
		current = Streamline_queue.begin();
		std::advance(current,counter);
		for(auto& p : current -> points) {
			for (auto& dir : Directions) {
				seedpoint = p + dir;
				bool valid = check_valid_seed(seedpoint);
				if (valid) {
					make_streamline(seedpoint);
				}
			}
		}
		current = Streamline_queue.begin();
		std::advance(current,counter);
		if (current == Streamline_queue.end()) done = true;
		else {
			counter++;
		}
	}
}

void ESADStreamlines::print_vel() {
	for (auto& row : Velocities) {
		for (auto& vel : row) {
			std::cout << vel.transpose() << " ";
		}
		std::cout << "\n";
	}
}

void ESADStreamlines::print_streamline(Streamline& line){
	for (auto& point : line.points) {
		std::cout << point.transpose() << std::endl;
	}
}

void ESADStreamlines::print_CellList() {
	std::cout << "printing celllist..." << std::endl;
	for (auto& row : CellList) {
		for (auto& list : row) {
			for (auto& el : list) {
				std::cout << el.transpose() << std::endl;
			}
		}
	}
}
