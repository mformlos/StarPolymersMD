#include "Box.h"

Box::Box(double Lx, double Ly, double Lz, double temperature, double lambda) :
	SystemTime { },
	Temperature { temperature },
	Lambda { lambda },
	Cutoff { 1.5 },
	VerletRadius { 2.0 },
	VerletRadius2 { 4.0 },
	NumberOfMonomers { } {
		BoxSize[0] = Lx;
		BoxSize[1] = Ly;
		BoxSize[2] = Lz;
		CellSize[0] = int(BoxSize[0]/VerletRadius);
		CellSize[1] = int(BoxSize[1]/VerletRadius);
		CellSize[2] = int(BoxSize[2]/VerletRadius);
		CellSideLength[0] = BoxSize[0]/(double)CellSize[0];
		CellSideLength[1] = BoxSize[1]/(double)CellSize[1];
		CellSideLength[2] = BoxSize[2]/(double)CellSize[2];
		CellList = std::vector<std::vector<std::vector<std::forward_list<MDParticle*>>>>(CellSize[0], std::vector<std::vector<std::forward_list<MDParticle*>>>(CellSize[1], std::vector<std::forward_list<MDParticle*>>(CellSize[2], std::forward_list<MDParticle*>())));
}


void Box::add_chain(unsigned A, unsigned B, double mass, double bondLength) {
	Molecules.push_back(Molecule{A+B, mass});
	Molecules.back().initialize_straight_chain(A, B, Temperature, bondLength);
	wrap(Molecules.back());
	NumberOfMonomers += (A+B);
}

void Box::add_star(unsigned A, unsigned B, unsigned Arms, double Mass, double Bond, double AnchorBond) {
	Molecules.push_back(Molecule{(A+B)*Arms+1, Mass});
	Molecules.back().initialize_open_star(A, B, Arms, Temperature, Bond, AnchorBond);
	wrap(Molecules.back());
	NumberOfMonomers += (A+B)*Arms + 1;
}

inline Vector3d& Box::wrap(Vector3d& pos) {
	for (unsigned i = 0; i < 3; i++) {
		pos(i) -= floor(pos(i)/BoxSize[i]) * BoxSize[i];
	}
	return pos;
}

inline Vector3d Box::wrap(Vector3d&& pos) {
	for (unsigned i = 0; i < 3; ++i) {
		pos(i) -= floor(pos(i)/BoxSize[i]) * BoxSize[i];
	}
	return pos;
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

Vector3d Box::relative_position(Particle& one, Particle& two) {
	Vector3d result {two.Position - one.Position};
	for (unsigned i = 0; i < 3; ++i) {
		result(i) -= round(result(i)/BoxSize[i]) * BoxSize[i];
	}
	return result;
}

void Box::update_VerletLists() {
	std::array<int,3> CellNumber { };
	double radius2 { };
	Vector3d distance { };
	//clear all Lists;
	for (auto& sheet : CellList) {
		for (auto& row : sheet) {
			for (auto& list : row) {
				list.clear();
			}
		}
	}
	for (auto& mol : Molecules) {
		for (auto& mono : mol.Monomers) {
			mono.clear_VerletList();
		}
	}

	//sort into CellLists
	for (auto& mol: Molecules) {
		for (auto& mono : mol.Monomers) {
			for (int i = 0; i < 3; ++i) {
				CellNumber[i] = (int)(mono.Position(i)/CellSideLength[i]);
			}
			mono.VerletPosition = mono.Position;
			CellList[CellNumber[0]][CellNumber[1]][CellNumber[2]].push_front(&mono);
		}
	}

	//make VerletLists
	int p { }, q { }, r { };
	for (auto& mol : Molecules) {
		for (auto& mono : mol.Monomers) {
			for (int i = 0; i < 3; ++i) {
				CellNumber[i] = (int)(mono.Position(i)/CellSideLength[i]);
			}
			for (int j = CellNumber[0]-1; j < CellNumber[0]+2; ++j) {
				for (int k = CellNumber[1]-1; k < CellNumber[1]+2; ++k) {
					for (int l = CellNumber[2]-1; l < CellNumber[2]+2; ++l) {

						p = my_modulus(j, CellSize[0]);
						q = my_modulus(k, CellSize[1]);
						r = my_modulus(l, CellSize[2]);

						for (auto& other : CellList[p][q][r]) {
							if (other == &mono) continue;
							distance = relative_position(mono, *other);
							radius2 = distance.dot(distance);
							if (radius2 <= VerletRadius2) {
								mono.VerletList.push_front(other);
							}
						}
					}
				}
			}
		}
	}
}

void Box::check_VerletLists() {
	Vector3d displacement {Vector3d::Zero()};
	for (auto& mol : Molecules) {
		for (auto& mono : mol.Monomers) {
			displacement = mono.Position - mono.VerletPosition;
			for (unsigned i = 0; i < 3; ++i) {
				displacement(i) -= round(displacement(i)/BoxSize[i]) * BoxSize[i];
			}
			if (displacement.norm() > (VerletRadius - Cutoff)*0.5) {
				update_VerletLists();
				return;
			}
		}
	}
	return;
}

void Box::calculate_forces(bool calc_epot) {
	double radius2 { };
	double force_abs { };
	Vector3d distance { };
	Vector3d force { };
	for (auto& mol : Molecules) {
		if (calc_epot) mol.Epot = 0.0;
		for (auto& mono : mol.Monomers) mono.Force *= 0.0;
		for (unsigned i = 0; i < mol.Monomers.size(); i++) {
			for (unsigned j = i + 1; j < mol.Monomers.size(); j++) {
				distance = relative_position(mol[i], mol[j]);
				radius2 = distance.dot(distance);
				if (mol[i].AmphiType == 1 && mol[j].AmphiType == 1) { // BB Type
					if (calc_epot) mol.Epot += TypeBB_Potential(radius2, Lambda);
					force_abs = TypeBB_Force(radius2, Lambda);
					force = distance*force_abs;
					mol[i].Force -= force;
					mol[j].Force += force;
				}
				else { // AA Type
					if (calc_epot) mol.Epot += TypeAA_Potential(radius2);
					force_abs = TypeAA_Force(radius2);
					if (force_abs > 1e5) std::cout << "AA! " << i <<  " " << j << " " << radius2 <<  std::cout;
					force = distance*force_abs;
					mol[i].Force -= force;
					mol[j].Force += force;
				}
			}
			for (auto& neighbor : mol[i].Neighbors) { //Fene bonds
				distance = relative_position(mol[i], *neighbor);
				radius2 = distance.dot(distance);
				if (mol[i].Anchor) {
					if (calc_epot) mol.Epot += Fene_Anchor_Potential(radius2);
					force_abs = Fene_Anchor_Force(radius2);
				}
				else {
					if (calc_epot) mol.Epot += Fene_Potential(radius2);
					force_abs = Fene_Force(radius2);
				}
				force = distance*force_abs;
				mol[i].Force -= force;
				neighbor -> Force += force;

			}
		}
	}
}

void Box::calculate_forces_verlet(bool calc_epot) {
	//double radius2 { };
	double force_abs { };
	//Vector3d distance { };
	Vector3d force { };
	for (auto& mol : Molecules) {
		mol.Epot = 0.0;
		for (auto& mono : mol.Monomers) mono.Force *= 0.0;
	}
	for (auto& mol : Molecules) {
		for (auto& mono : mol.Monomers) {

			for (auto& other : mono.VerletList) {
				Vector3d distance {relative_position(mono, *other)};
				double radius2 {distance.dot(distance)};
				//distance = relative_position(mono, *other);
				//radius2 = distance.dot(distance);
				if (mono.AmphiType == 1 && other -> AmphiType == 1) { // BB Type
					if (calc_epot) mol.Epot += 0.5*TypeBB_Potential(radius2, Lambda);
					force_abs = TypeBB_Force(radius2, Lambda);
					force = distance*force_abs;
					mono.Force -= force;
				}
				else { // AA Type
					if (calc_epot) mol.Epot += 0.5*TypeAA_Potential(radius2);
					force_abs = TypeAA_Force(radius2);
					if (force_abs > 1e5) std::cout << "AA! " << mono.Position.transpose() <<  " " << other -> Position.transpose() << " " << radius2 <<  std::cout;

					force = distance*force_abs;
					mono.Force -= force;
				}
			}
			for (auto& neighbor : mono.Neighbors) { //Fene bonds
				Vector3d distance {relative_position(mono, *neighbor)};
				double radius2 {distance.dot(distance)};
				//distance = relative_position(mono, *neighbor);
				//radius2 = distance.dot(distance);
				if (mono.Anchor) {
					if (calc_epot) mol.Epot += Fene_Anchor_Potential(radius2);
					force_abs = Fene_Anchor_Force(radius2);
				}
				else {
					if (calc_epot) mol.Epot += Fene_Potential(radius2);
					force_abs = Fene_Force(radius2);
				}
				force = distance*force_abs;
				mono.Force -= force;
				neighbor -> Force += force;
			}
		}
	}
}

double Box::calculate_ekin() {
	double KineticEnergy { };
	for (auto& mol : Molecules) {
		KineticEnergy += mol.calculate_Ekin();
	}
	return KineticEnergy;
}

double Box::calculate_epot(MDParticle& part1, MDParticle& part2) {
	double epot { };
	Vector3d distance = relative_position(part1, part2);
	double radius2 = distance.dot(distance);
	if (part1.AmphiType == 1 && part2.AmphiType == 1) epot += TypeBB_Potential(radius2, Lambda);
	else epot += TypeAA_Potential(radius2);
	for (auto& neighbor : part1.Neighbors) {
		if (neighbor == &part2) {
			if (part1.Anchor || part2.Anchor) epot += Fene_Anchor_Potential(radius2);
			else epot += Fene_Potential(radius2);
		}
	}
	return epot;
}

double Box::calculate_radius_of_gyration() {
	double r_gyr { };
	Vector3d distance { };
	for (auto& mol : Molecules) {
		double r_gyr_mol { };
		for (unsigned i = 0; i < mol.Monomers.size(); ++i) {
			for (unsigned j = i+1; j < mol.Monomers.size(); ++j) {
				Vector3d distance {relative_position(mol[i], mol[j])};
				r_gyr_mol += distance.dot(distance);
			}
		}
		r_gyr_mol /= (mol.Monomers.size()*mol.Monomers.size());
		r_gyr_mol = sqrt(r_gyr_mol);
		r_gyr += r_gyr_mol;
	}
	r_gyr /= Molecules.size();
	return r_gyr;
}

Matrix3d Box::calculate_gyration_tensor() {
	Matrix3d gyr_tensor {Matrix3d::Zero()};
	for (auto& mol: Molecules) {
		Matrix3d gyr_tensor_mol {Matrix3d::Zero()};
		Vector3d shift_anchor_to_center {Vector3d::Zero()};
		Vector3d center_of_mass {Vector3d::Zero()};
		for (int i = 0; i < 3; i++) {
			shift_anchor_to_center(i) = BoxSize[i]*0.5 - mol.Monomers[0].Position(i);
		}
		for (auto& mono : mol.Monomers) {
			center_of_mass += (mono.Position + shift_anchor_to_center);
		}
		center_of_mass /= mol.NumberOfMonomers;
		for (auto& mono : mol.Monomers) {
			Vector3d shifted_position = mono.Position +shift_anchor_to_center - center_of_mass;
			for (int alpha = 0; alpha < 3; alpha++) {
				for (int beta = 0; beta < 3; beta++) {
					gyr_tensor_mol(alpha, beta) += shifted_position(alpha)*shifted_position(beta);
				}
			}
		}
		gyr_tensor_mol /= mol.NumberOfMonomers;
		gyr_tensor += gyr_tensor_mol;
	}
	gyr_tensor /= Molecules.size();
	return gyr_tensor;
}

std::list<unsigned> Box::calculate_clusters() {
	std::list<unsigned> cluster_sizes { };
	for (auto& mol : Molecules) {
		std::list<MDParticle*> copy_Particles { };
		for(auto& mono : mol.Monomers) {
			copy_Particles.push_back(&mono);
		}
		std::list<MDParticle*> Particles_in_cluster { };
		while (!copy_Particles.empty()) {
			MDParticle* part1 { };
			if (Particles_in_cluster.empty()) {
				part1 = &*copy_Particles.back();
				copy_Particles.pop_back();

				cluster_sizes.push_back(1);
			}
			else {
				part1 = &*Particles_in_cluster.back();
				Particles_in_cluster.pop_back();
			}
			for (auto& part2 : part1 -> VerletList) {
				auto iter = std::find(copy_Particles.begin(), copy_Particles.end(), part2);
 				if (iter != copy_Particles.end()){
					if (calculate_epot(*part1, *part2) < 0.0) {
						Particles_in_cluster.splice(Particles_in_cluster.end(), copy_Particles, iter);
						cluster_sizes.back()++;
					}
				}
			}
		}
	}
	cluster_sizes.sort();
	return cluster_sizes;
}

std::list<unsigned> Box::calculate_patches() {
	std::list<unsigned> patches { };
	for (auto& mol : Molecules) {
		std::list<unsigned> copy_particles { };
		std::list<unsigned> copy_arms { };
		for (unsigned i = 1; i < mol.NumberOfMonomers; i++) copy_particles.push_back(i);
		for (unsigned i = 0; i < mol.Arms; i++) copy_arms.push_back(i);

		std::list<unsigned> arms_in_patch { };
		unsigned f1 {};
		while (!copy_arms.empty()){
			if (arms_in_patch.empty()) {
				f1 = copy_arms.front();
				copy_arms.pop_front();
				patches.push_back(1);
			}
			else {
				f1 = arms_in_patch.front();
				arms_in_patch.pop_front();
			}
			for (unsigned i = 1; i <= (mol.AType+mol.BType); i++) {
				unsigned part1 { f1*(mol.AType+mol.BType) + i };
				copy_particles.remove(part1);
				std::list<unsigned>::iterator part_it = copy_particles.begin();
				while (part_it != copy_particles.end()){
					unsigned part2 {*part_it};
					if (part2 < part1 || part2 > (mol.AType+mol.BType)*(f1+1)) {
						if (calculate_epot(mol.Monomers[part1], mol.Monomers[part2]) < 0.0) {
							unsigned f2 {unsigned((part2-1)/(mol.AType+mol.BType))};
							int reg {(part2-1)%(mol.AType+mol.BType)};
							arms_in_patch.push_back(f2);
							copy_arms.remove(f2);
							patches.back()++;
							std::list<unsigned>::iterator erase_begin {part_it};
							std::advance(erase_begin, -reg);
							std::list<unsigned>::iterator erase_end {erase_begin};
							std::advance(erase_end, mol.AType+mol.BType);
							part_it = copy_particles.erase(erase_begin, erase_end);

						}
						else part_it++;
					}
					else part_it++;
				}
			}
		}
	}
	patches.sort();
	return patches;
}

unsigned Box::numberOfMonomers() {
	unsigned N { };
	for(auto& mol: Molecules) {
		N += mol.Monomers.size();
	}
	return N;
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
	PotentialEnergy /= NumberOfMonomers;
	os << PotentialEnergy<< " ";
	return os;
}

std::ostream& Box::print_Ekin(std::ostream& os) {
	double KineticEnergy { };
	for (auto& mol : Molecules) {
		KineticEnergy += mol.calculate_Ekin();
	}
	KineticEnergy /= NumberOfMonomers;
	os << KineticEnergy << " ";
	return os;
}

std::ostream& Box::print_Temperature(std::ostream& os) {
	double Temperature { };
	for (auto& mol : Molecules) {
		Temperature += mol.calculate_Ekin()/mol.NumberOfMonomers;
	}
	Temperature *= 2./3.;
	os << Temperature << " ";
	return os;
}

std::ostream& Box::print_radius_of_gyration(std::ostream& os) {
	os << calculate_radius_of_gyration() << " ";
	return os;
}


template<class UnitaryFunc>
UnitaryFunc Box::unitary(UnitaryFunc&& func) const {
	return for_each(Molecules.cbegin(), Molecules.cend(), func);
}

template<class UnaryFunc, class BinaryFunc>
void Box::operator() (UnaryFunc& ufunc, BinaryFunc& bfunc) const {
	auto first = Molecules.cbegin(), last = Molecules.cend();
	for(; first != last; ++first) {
		*first(ufunc, bfunc);
	}
}

