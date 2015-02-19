#include "Box.h"

Box::Box(double Lx, double Ly, double Lz, double temperature) :
	SystemTime { },
	Temperature { temperature } {
		Size[0] = Lx;
		Size[1] = Ly;
		Size[2] = Lz;
		std::array<int, 3> cellSize { int(Size[0]/1.5), int(Size[1]/1.5), int(Size[2]/1.5)};
		CellList = std::vector<std::vector<std::vector<std::forward_list<Particle*>>>>(cellSize[0], std::vector<std::vector<std::forward_list<Particle*>>>(cellSize[1], std::vector<std::forward_list<Particle*>>(cellSize[1], std::forward_list<Particle*>())));
}

void Box::add_chain(unsigned N, double mass, double bondLength) {
	Molecule chain{N, mass};
	chain.initialize_straight_chain(bondLength);
	Molecules.push_back(chain);
}

std::ostream& Box::print_molecules(std::ostream& os) const {
	for (auto& mol : Molecules) {
		os << mol << std::endl;
	}
	return os;
}


