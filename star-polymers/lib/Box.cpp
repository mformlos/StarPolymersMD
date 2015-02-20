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



inline MatVec& Box::wrap(MatVec& pos) {
	return pos -=  floor(pos/Size) % Size;
}

inline MatVec Box::wrap(MatVec&& pos) {
	return pos -=  floor(pos/Size) % Size;
}

inline void Box::wrap(Particle& part) {
	part.Position = wrap(part.Position);
}

inline void Box::wrap(Molecule& mol) {
	for (auto& mono : mol.Monomers) wrap(mono);
}

void Box::wrap() {
	for (auto& mol : Molecules) wrap(mol);
}

MatVec Box::relative_position(Particle& one, Particle& two) {
	MatVec result { };
	result = two.Position - one.Position;
	result -= round(result/Size) % Size;
	return result;
}

void Box::add_chain(unsigned N, double mass, double bondLength) {
	Molecules.push_back(Molecule{N, mass});
	Molecules.back().initialize_straight_chain(bondLength);
	wrap(Molecules.back());
	calculate_forces();
}

std::ostream& Box::print_molecules(std::ostream& os) const {
	for (auto& mol : Molecules) {
		os << mol << std::endl;
	}
	return os;
}

std::ostream& Box::print_Epot(std::ostream& os) const {
	double PotentialEnergy { };
	for (auto& mol : Molecules) {
		PotentialEnergy += mol.Epot;
	}
	os << PotentialEnergy << std::endl;
	return os;
}

void Box::calculate_forces() {
	double radius2 { };
	double force_abs { };
	MatVec distance { };
	MatVec force { };
	for (auto& mol : Molecules) {
		mol.Epot = 0.0;
		for (auto& mono : mol.Monomers) mono.Force *= 0.0;
		for (unsigned i = 0; i < mol.Monomers.size(); i++) {
			for (unsigned j = i + 1; j < mol.Monomers.size(); j++) {
				distance = relative_position(mol[i], mol[j]);
				radius2 = distance*distance;
				if (mol[i].AmphiType == 1 && mol[j].AmphiType == 1) { // BB Type
					mol.Epot += TypeBB_Potential(radius2, 1.0);
					force_abs = TypeBB_Force(radius2, 1.0);
					force = distance*force_abs;
					mol[i].Force += force;
					mol[j].Force -= force;
				}
				else { // AA Type
					mol.Epot += TypeAA_Potential(radius2);
					force_abs = TypeAA_Force(radius2);
					force = distance*force_abs;
					mol[i].Force += force;
					mol[j].Force -= force;
				}
				for (auto& neighbor : mol[i].Neighbors) { //Fene bonds
					if (neighbor == &(mol.Monomers[j])) {
						mol.Epot += Fene_Potential(radius2);
						force_abs = Fene_Force(radius2);
						force = distance*force_abs;
						mol[i].Force += force;
						mol[j].Force -= force;
					}
				}
			}
		}
	}
}


