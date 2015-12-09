#include "Box.h"
#include "MPC.h"

Box::Box(int Lx, int Ly, int Lz, double temperature, double lambda) :
	SystemTime { },
	Temperature { temperature },
	Lambda { lambda },
	Cutoff { 1.5 },
	VerletRadius { 2.0 },
	VerletRadius2 { 4.0 },
	NumberOfMonomers { },
	COM_Pos {Vector3d::Zero()},
	COM_Vel {Vector3d::Zero()} {
		if (Lx%2) Lx++;
		if (Ly%2) Ly++;
		if (Lz%2) Lz++;
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
	//wrap(Molecules.back());
	NumberOfMonomers += (A+B);
	Molecules_com_reference_frame.push_back(Molecule{Molecules.back()});
}

void Box::add_star(unsigned A, unsigned B, unsigned Arms, double Mass, double Bond, double AnchorBond) {
	Molecules.push_back(Molecule{(A+B)*Arms+1, Mass});
	Molecules.back().initialize_open_star(A, B, Arms, Temperature, Bond, AnchorBond);
	set_center_of_mass_to_zero(Molecules.back());
	NumberOfMonomers += (A+B)*Arms + 1;
	Molecules_com_reference_frame.push_back(Molecule{Molecules.back()});
	center_of_mass_reference_frame();
}

void Box::add_star(string filename, unsigned A, unsigned B, unsigned Arms, double Mass, bool set_zero) {
	Molecules.push_back(Molecule{(A+B)*Arms+1, Mass});
	Molecules.back().star_from_file(filename, A, B, Arms, set_zero);
	NumberOfMonomers += (A+B)*Arms + 1;
	Molecules_com_reference_frame.push_back(Molecule{Molecules.back()});
	center_of_mass_reference_frame();
}


void Box::resize(double Lx, double Ly, double Lz) {
	BoxSize[0] = Lx;
	BoxSize[1] = Ly;
	BoxSize[2] = Lz;
}

void Box::set_center_of_mass_to_zero(Molecule& mol) {
	Vector3d com {mol.calculate_center_of_mass()};
	for (auto& mono : mol.Monomers) {
		mono.Position -= com;
	}
}

void Box::center_of_mass_reference_frame() {
	for (unsigned n = 0; n < Molecules.size(); n++) {
		Molecule centered_molecule {Molecules[n]};
		Vector3d center_of_mass {Vector3d::Zero()};
		Vector3d center_of_mass_velocity {Vector3d::Zero()};
		for (auto& mono : centered_molecule.Monomers) {
			center_of_mass += mono.Position;
			center_of_mass_velocity += mono.Velocity;
		}
		center_of_mass /= centered_molecule.NumberOfMonomers;
		center_of_mass_velocity /= centered_molecule.NumberOfMonomers;
		for(auto& mono : centered_molecule.Monomers) {
			mono.Position -= center_of_mass;
			mono.Velocity -= center_of_mass_velocity;
		}
		if (n ==0) {
			COM_Pos = center_of_mass; //- shift_anchor_to_center;
			COM_Vel = center_of_mass_velocity;
		}
		Molecules_com_reference_frame[n] = centered_molecule;
	}
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
	center_of_mass_reference_frame();

	for (auto& mol: Molecules) {
		for (auto& mono : mol.Monomers) {
			for (int i = 0; i < 3; ++i) {
				CellNumber[i] = (int)((mono.Position(i)-COM_Pos(i)+BoxSize[i]*0.5)/CellSideLength[i]);
			}
			mono.VerletPosition = mono.Position;
			if (CellNumber[0] < 0 || CellNumber[0] >= CellSize[0] || CellNumber[1] < 0 || CellNumber[1] >= CellSize[1] || CellNumber[2] < 0 || CellNumber[2] >= CellSize[2]) {
				std::cout << mono.Position.transpose() << " " << COM_Pos.transpose() << std::endl;
				FILE* dump_file { };
				dump_file = fopen("core_dump.pdb", "a");
				print_PDB_with_velocity(dump_file,0);
			}
			CellList[CellNumber[0]][CellNumber[1]][CellNumber[2]].push_front(&mono);
		}
	}

	//make VerletLists
	int p { }, q { }, r { };
	for (auto& mol : Molecules) {
		for (auto& mono : mol.Monomers) {
			for (int i = 0; i < 3; ++i) {
				CellNumber[i] = (int)((mono.Position(i)-COM_Pos(i)+BoxSize[i]*0.5)/CellSideLength[i]);
			}
			for (int j = CellNumber[0]-1; j < CellNumber[0]+2; ++j) {
				for (int k = CellNumber[1]-1; k < CellNumber[1]+2; ++k) {
					for (int l = CellNumber[2]-1; l < CellNumber[2]+2; ++l) {

						p = my_modulus(j, CellSize[0]);
						q = my_modulus(k, CellSize[1]);
						r = my_modulus(l, CellSize[2]);

						for (auto& other : CellList[p][q][r]) {
							if (other == &mono) continue;
							distance = other -> Position - mono.Position;
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
				distance = mol[j].Position - mol[i].Position;
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
					force = distance*force_abs;
					mol[i].Force -= force;
					mol[j].Force += force;
				}
			}
			for (auto& neighbor : mol[i].Neighbors) { //Fene bonds
				distance = neighbor -> Position - mol[i].Position;
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
	double force_abs { };
	Vector3d force { };
	for (auto& mol : Molecules) {
		mol.Epot = 0.0;
		for (auto& mono : mol.Monomers) mono.Force *= 0.0;
	}
	for (auto& mol : Molecules) {
		for (auto& mono : mol.Monomers) {

			for (auto& other : mono.VerletList) {
				Vector3d distance {other -> Position - mono.Position};
				double radius2 {distance.dot(distance)};
				if (mono.AmphiType == 1 && other -> AmphiType == 1) { // BB Type
					if (calc_epot) mol.Epot += 0.5*TypeBB_Potential(radius2, Lambda);
					force_abs = TypeBB_Force(radius2, Lambda);
					force = distance*force_abs;
					mono.Force -= force;
				}
				else { // AA Type
					if (calc_epot) mol.Epot += 0.5*TypeAA_Potential(radius2);
					force_abs = TypeAA_Force(radius2);
					force = distance*force_abs;
					mono.Force -= force;
				}
			}
			for (auto& neighbor : mono.Neighbors) { //Fene bonds
				Vector3d distance {neighbor -> Position - mono.Position};//{relative_position(mono, *neighbor)};
				double radius2 {distance.dot(distance)};
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
	Vector3d distance {part2.Position - part1.Position};
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
	for (auto& mol : Molecules_com_reference_frame) {
		double r_gyr_mol { };
		for (unsigned i = 0; i < mol.Monomers.size(); ++i) {
			for (unsigned j = i+1; j < mol.Monomers.size(); ++j) {
				Vector3d distance {mol[i].Position - mol[j].Position};
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

std::tuple<double, Matrix3d> Box::calculate_gyration_tensor() {
	Matrix3d gyr_tensor {Matrix3d::Zero()};
	double r_gyr { };
	for (auto& mol: Molecules_com_reference_frame) {
		double r_gyr_mol { };
		Matrix3d gyr_tensor_mol {Matrix3d::Zero()};
		for (auto& mono : mol.Monomers) {
			for (int alpha = 0; alpha < 3; alpha++) {
				for (int beta = 0; beta < 3; beta++) {
					gyr_tensor_mol(alpha, beta) += mono.Position(alpha)*mono.Position(beta);
				}
			}
		}
		gyr_tensor_mol /= mol.NumberOfMonomers;
		gyr_tensor += gyr_tensor_mol;

		r_gyr_mol = gyr_tensor_mol(0,0) + gyr_tensor_mol(1,1) + gyr_tensor_mol(2,2);
		r_gyr_mol = sqrt(r_gyr_mol);
		r_gyr += r_gyr_mol;
	}
	r_gyr /= Molecules.size();
	gyr_tensor /= Molecules.size();
	return std::make_tuple(r_gyr, gyr_tensor);
}

Vector3d Box::calculate_rotation_frequency() {
	Vector3d omega {Vector3d::Zero()};
	for (auto& mol : Molecules_com_reference_frame) {
		Vector3d omega_mol {Vector3d::Zero()};
		Matrix3d inertia_tensor {Matrix3d::Zero()};
		Vector3d angular_momentum {Vector3d::Zero()};
		for (auto& mono : mol.Monomers) {
			double rsqr {mono.Position.squaredNorm()};
			for (unsigned i = 0; i < 3; i++) {
				for (unsigned j = i; j < 3; j++) {
					inertia_tensor(i,j) -= mono.Position(i)*mono.Position(j);
					if (i == j) inertia_tensor(i,j) += rsqr;
				}
			}
			angular_momentum += mono.Position.cross(mono.Velocity);
			inertia_tensor(1,0) = inertia_tensor(0,1);
			inertia_tensor(2,1) = inertia_tensor(1,2);
			inertia_tensor(2,0) = inertia_tensor(0,2);
		}
		omega_mol = inertia_tensor.ldlt().solve(angular_momentum);
		omega += omega_mol;
	}
	omega /= Molecules.size();
	return omega;
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

std::tuple<double, double> Box::calculate_patches_new() {
	double patch_number { };
	double patch_size { };
	for (auto& mol : Molecules) {
		int patch_number_mol { };
		double patch_size_mol { };
		std::vector<std::vector<bool>> arm_bond(mol.Arms, std::vector<bool>(mol.Arms, false));
		for (unsigned i = 0; i < mol.Arms; i++) {    // for every arm
			for (unsigned k = 0; k < mol.Arms; k++) {   // for every Arm j != i
                if (arm_bond[i][k] || i==k) continue;
				for (unsigned j = mol.AType + 1; j <= mol.AType+mol.BType; j++) {   //for every attractive monomer in Arm i
					unsigned j_n { i*(mol.AType+mol.BType)+j };
					for (unsigned l = mol.AType +1; l <= mol.AType+mol.BType; l++) { // ever attractive monomer in Arm j
						unsigned l_n {k*(mol.AType+mol.BType)+l};
						if (calculate_epot(mol.Monomers[j_n], mol.Monomers[l_n]) < 0.0) { //stop loop if one attractive energy is found
							arm_bond[i][k] = true;
							break;
						}
					}
				}
			}
		}// arm_bond matrix filled

		std::vector<std::vector<unsigned>> clusters(mol.Arms, std::vector<unsigned>());
		for (unsigned i = 0; i < mol.Arms; i++) clusters[i].push_back(i); //each arm is one cluster


		for (unsigned i = 0; i < mol.Arms; i++) { //loop trough cluster i
			std::vector<unsigned>::iterator iter {};
			if (clusters[i].empty()) continue;
			iter = clusters[i].begin();
			int n {0};
			do {
				unsigned arm_first{};
				if (clusters[i].empty()) continue;
				arm_first = *iter; // first arm in cluster
				for (unsigned j = i+1; j < mol.Arms; j++) { //loop trough all clusters j > i
					 for (auto& arm_second : clusters[j]) { //loop through all arms in cluster j

						 if (arm_bond[arm_first][arm_second]) {
							 clusters[i].insert(clusters[i].end(), clusters[j].begin(), clusters[j].end()); //move all to cluster i
							 clusters[j].clear(); //delete all in j
							 if (n > 0) iter = clusters[i].begin()+n;
							 else iter = clusters[i].begin();
						 }
					 }

				 }
				 n++;
			     if (iter != clusters[i].end()) iter++;
			}while(iter != clusters[i].end());

		}
		for (auto& cluster : clusters) {
			if (cluster.size()>1) {
				patch_number_mol++;
				for (auto& arm : cluster) {
					patch_size_mol += 1.0;
				}
			}
		}
		if (patch_number_mol > 1) patch_size_mol /= patch_number_mol;
		patch_number += patch_number_mol;
		patch_size += patch_size_mol;
	}
	patch_number /= Molecules.size();
	patch_size /= Molecules.size();
	return std::make_tuple(patch_number, patch_size);
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

void Box::print_PDB(FILE* pdb, int step) {
	int mol_count {0};
	fprintf(pdb, "MODEL     %d \n", step);
	//for (auto& mol : Molecules_com_reference_frame) {
	for (auto& mol : Molecules) {
		mol_count++;
		for (unsigned i = 0; i < mol.NumberOfMonomers; i++) {
			Vector3d pos_print {mol.Monomers[i].Position};
			if (mol.Monomers[i].AmphiType) fprintf(pdb, "ATOM %6d  C   GLY    %2d     %7.3f %7.3f %7.3f \n", i+1, mol_count, pos_print(0), pos_print(1), pos_print(2));
			else fprintf(pdb, "ATOM %6d  O   GLY    %2d     %7.3f %7.3f %7.3f \n", i+1, mol_count, pos_print(0), pos_print(1), pos_print(2));
		}
		fprintf(pdb, "TER \n");
		unsigned count{1};
		unsigned arm {0};
		while (arm < mol.Arms) {
			fprintf(pdb, "CONECT %4i ", 1);

			while(count%5) {
				fprintf(pdb, "%4i ", arm*(mol.AType+mol.BType)+2);
				arm++;
				count++;
				if (arm >= mol.Arms) break;
			}
			count = 1;
			fprintf(pdb, "\n");
		}
		for (unsigned i = 0; i < mol.Arms; i++) {
			fprintf(pdb, "CONECT %4u %4u \n", i*(mol.AType+mol.BType) + 2, 1); // , i*(mol.AType+mol.BType) + 3);
			for(unsigned j = 1; j < (mol.AType+mol.BType); j++) {
				unsigned number{i*(mol.AType+mol.BType) + j};
				fprintf(pdb, "CONECT %4u %4u %4u \n", number +1 , number, number+2);
			}
			fprintf(pdb, "CONECT %4u %4u \n", (i+1)*(mol.AType+mol.BType) +1 , (i+1)*(mol.AType+mol.BType));
		}
	}
	fprintf(pdb, "ENDMDL \n");

}

void Box::print_PDB_with_velocity(FILE* pdb, int step) {
	int mol_count {0};
	fprintf(pdb, "MODEL     %d \n", step);
	//for (auto& mol : Molecules_com_reference_frame) {
	for (auto& mol : Molecules) {
    	mol_count++;
		for (unsigned i = 0; i < mol.NumberOfMonomers; i++) {
            Vector3d pos_print {mol.Monomers[i].Position};
			fprintf(pdb, "ATOM %6d  C   GLY    %2d     %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f \n", i, mol_count, pos_print(0), pos_print(1), pos_print(2), mol.Monomers[i].Velocity(0), mol.Monomers[i].Velocity(1), mol.Monomers[i].Velocity(2));
		}
		fprintf(pdb, "TER \n");
		unsigned count{1};
		unsigned arm {0};
		while (arm < mol.Arms) {
			fprintf(pdb, "CONECT %4i ", 1);

			while(count%5) {
				fprintf(pdb, "%4i ", arm*(mol.AType+mol.BType)+2);
				arm++;
				count++;
				if (arm >= mol.Arms) break;
			}
			count = 1;
			fprintf(pdb, "\n");
		}
		for (unsigned i = 0; i < mol.Arms; i++) {
			fprintf(pdb, "CONECT %4u %4u \n", i*(mol.AType+mol.BType) + 2, 1);
			for(unsigned j = 1; j < (mol.AType+mol.BType); j++) {
				unsigned number{i*(mol.AType+mol.BType) + j};
				fprintf(pdb, "CONECT %4u %4u %4u \n", number +1 , number, number+2);
			}
			fprintf(pdb, "CONECT %4u %4u \n", (i+1)*(mol.AType+mol.BType) +1 , (i+1)*(mol.AType+mol.BType));
		}
	}
	fprintf(pdb, "ENDMDL \n");

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

std::ostream& Box::print_Ekin_and_Temperature(std::ostream& os) {
	double Temperature { };
	double Ekin { };
	unsigned NumberOfMonomers { };
	for (auto& mol : Molecules) {
		Ekin += mol.calculate_Ekin();
		NumberOfMonomers += mol.NumberOfMonomers;
	}
	Temperature = Ekin*2./(3.*NumberOfMonomers);
	os << Ekin << " " << Temperature << " ";
	return os;
}


std::ostream& Box::print_radius_of_gyration(std::ostream& os) {
	os << calculate_radius_of_gyration() << " ";
	return os;
}

std::ostream& Box::print_center_of_mass(std::ostream& os) {
	for (auto& mol : Molecules) {
		Vector3d com {mol.calculate_center_of_mass()};
		os << com.transpose() << " " << mol.calculate_center_of_mass_velocity().transpose() << '\n'; // << (com-BoxCenter).norm() << " ";
	}
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

