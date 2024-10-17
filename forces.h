#ifndef FORCES_H_
#define FORCES_H_

#include "structs.h"

/**
 * @brief Calculate total forces on particles, including both bonded and non-bonded forces.
 * 
 * This function calculates the forces acting on all particles in the system by 
 * combining the contributions from both non-bonded interactions (e.g., Lennard-Jones potential)
 * and bonded interactions (e.g., bond-stretch, angle-bend, dihedral-torsion).
 * 
 * @param p_parameters used members: num_part, L
 * @param p_nbrlist used members: num_nbrs, nbr
 * @param p_vectors used members: r, f
 * @return double potential energy from all interactions
 */
double calculate_forces(struct Parameters *p_parameters, struct Nbrlist *p_nbrlist, struct Vectors *p_vectors);

/**
 * @brief Calculate non-bonded forces (e.g., Lennard-Jones interactions) between particles.
 * 
 * Non-bonded forces are computed based on the pairwise interactions between particles that
 * are not directly bonded. These forces are usually calculated using the Lennard-Jones potential,
 * with a cutoff distance applied to limit the interaction range. The Lorentz-Berthelot mixing rules
 * may be implemented for cross-species interactions in the case of multi-component systems.
 * 
 * @param p_parameters used members: num_part
 * @param p_nbrlist used members: num_nbrs, nbr
 * @param[out] p_vectors used members: f
 * @return double potential energy from non-bonded interactions
 */
double calculate_forces_nb(struct Parameters *p_parameters, struct Nbrlist *p_nbrlist, struct Vectors *p_vectors);

/**
 * @brief Calculate bond-stretch forces for directly bonded particles.
 * 
 * This function calculates the bond-stretch forces between particles that are connected by 
 * covalent bonds. A harmonic potential is typically used to model these interactions,
 * where the bond force is proportional to the deviation from the equilibrium bond length.
 * 
 * @param p_parameters used members: num_part, L
 * @param p_vectors used members: r, f
 * @return double potential energy from bond-stretch interactions
 */
double calculate_forces_bond(struct Parameters *p_parameters, struct Vectors *p_vectors);

/**
 * @brief Calculate angle-bend forces between sets of three bonded particles (angles).
 * 
 * This function computes the forces associated with the bending of angles between 
 * bonded triplets of particles (i.e., particles i-j-k). These forces are typically modeled using 
 * a harmonic potential, where the force depends on the deviation of the bond angle from its 
 * equilibrium value.
 * 
 * @param p_parameters used members: num_part, L
 * @param p_vectors used members: r, f
 * @return double potential energy from angle-bend interactions
 */
double calculate_forces_angle(struct Parameters *p_parameters, struct Vectors *p_vectors);

/**
 * @brief Calculate dihedral-torsion forces between sets of four bonded particles (dihedrals).
 * 
 * Dihedral forces arise from the torsional interaction between four connected particles 
 * (i-j-k-l). These forces are modeled using a potential that depends on the dihedral angle 
 * formed between the bonds. Dihedral potentials can be periodic and are critical for simulating
 * complex molecules like ethane.
 * 
 * @param p_parameters used members: num_part, L
 * @param p_vectors used members: r, f
 * @return double potential energy from dihedral-torsion interactions
 */
double calculate_forces_dihedral(struct Parameters *p_parameters, struct Vectors *p_vectors);

double generateZeta();

double calculate_pressure(struct Parameters *p_parameters, struct Vectors *p_vectors);

void calculate_density_profile(FILE *density_file, struct Parameters *p_parameters, struct Vectors *p_vectors, size_t step);

#endif /* FORCES_H_ */
