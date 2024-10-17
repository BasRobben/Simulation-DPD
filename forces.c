#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "forces.h"
#include "constants.h"
#include "structs.h"
#include "nbrlist.h"

// This function calculates all forces acting on the particles (bonded and non-bonded).
// It initializes the forces array, then calculates bond-stretch, angle-bend, dihedral-torsion,
// and non-bonded forces. The total potential energy is returned.
double calculate_forces(struct Parameters *p_parameters, struct Nbrlist *p_nbrlist, struct Vectors *p_vectors)
{
    struct Vec3D *f = p_vectors->f;
    size_t num_part = p_parameters->num_part;

    // Initialize the forces to zero for all particles
    for (size_t i = 0; i < num_part; i++)
        f[i] = (struct Vec3D){0.0, 0.0, 0.0};

    // Calculate the forces and accumulate the potential energy from each type of interaction
    double Epot = calculate_forces_bond(p_parameters, p_vectors);
    Epot += calculate_forces_angle(p_parameters, p_vectors);
    Epot += calculate_forces_dihedral(p_parameters, p_vectors);
    Epot += calculate_forces_nb(p_parameters, p_nbrlist, p_vectors);

    return Epot;
}

// This function calculates bond-stretch forces based on the current positions of the bonded particles.
// It applies the minimum image convention to calculate the distance between bonded pairs and then
// computes the force and potential energy due to the bond interaction.
double calculate_forces_bond(struct Parameters *p_parameters, struct Vectors *p_vectors){
    double Epot = 0.0;
    struct Bond *bonds = p_vectors->bonds;
    size_t num_bonds = p_vectors->num_bonds;
    struct Vec3D *f = p_vectors->f;
    struct Vec3D *r = p_vectors->r;
    struct Vec3D L = p_parameters->L;
    struct Vec3D rij;
    struct Vec3D fi = {0.0, 0.0, 0.0};
    double kT = p_parameters->kT;
    double r_cut = p_parameters->r_cut;
    double C = 2 * kT / (r_cut * r_cut);


    // Loop through each bond and calculate the forces
    for (size_t q = 0; q < num_bonds; ++q)
    {
        size_t i = bonds[q].i;
        size_t j = bonds[q].j;

        // Apply the minimum image convention for calculating distances
        rij.x = r[i].x - r[j].x;
        rij.x = rij.x - L.x * floor(rij.x / L.x + 0.5);
        rij.y = r[i].y - r[j].y;
        rij.y = rij.y - L.y * floor(rij.y / L.y + 0.5);
        rij.z = r[i].z - r[j].z;
        rij.z = rij.z - L.z * floor(rij.z / L.z + 0.5);

        // Calculate the spring force
        fi.x = -C * rij.x;
        fi.y = -C * rij.y;
        fi.z = -C * rij.z;

        // Update the forces on particles i and j
        f[i].x += fi.x;
        f[i].y += fi.y;
        f[i].z += fi.z;
        f[j].x -= fi.x;
        f[j].y -= fi.y;
        f[j].z -= fi.z;

        // Calculate the potential energy due to the bond interaction
        Epot += 0.5 * C * (rij.x * rij.x + rij.y * rij.y + rij.z * rij.z);
    }

    return Epot;  // Return the potential energy due to bond-stretch interactions
}

// This function calculates angle-bend forces based on the current positions of the angle-defined particles.
// It uses the minimum image convention and computes forces due to angle interactions.
double calculate_forces_angle(struct Parameters *p_parameters, struct Vectors *p_vectors){
    double Epot = 0.0;
    struct Angle *angles = p_vectors->angles;
    size_t num_angles = p_vectors->num_angles;
    struct Vec3D *f = p_vectors->f;
    struct Vec3D *r = p_vectors->r;
    struct Vec3D L = p_parameters->L;
    struct Vec3D rij, rkj;
    struct Vec3D fi = {0.0, 0.0, 0.0}, fk = {0.0, 0.0, 0.0};

    // // Loop through each angle and calculate the forces
    // for (size_t q = 0; q < num_angles; ++q)
    // {
    //     size_t i = angles[q].i;
    //     size_t j = angles[q].j;
    //     size_t k = angles[q].k;

    //     // Apply the minimum image convention for calculating distances
    //     rij.x = r[i].x - r[j].x;
    //     rij.x = rij.x - L.x * floor(rij.x / L.x + 0.5);
    //     rij.y = r[i].y - r[j].y;
    //     rij.y = rij.y - L.y * floor(rij.y / L.y + 0.5);
    //     rij.z = r[i].z - r[j].z;
    //     rij.z = rij.z - L.z * floor(rij.z / L.z + 0.5);

    //     rkj.x = r[k].x - r[j].x;
    //     rkj.x = rkj.x - L.x * floor(rkj.x / L.x + 0.5);
    //     rkj.y = r[k].y - r[j].y;
    //     rkj.y = rkj.y - L.y * floor(rkj.y / L.y + 0.5);
    //     rkj.z = r[k].z - r[j].z;
    //     rkj.z = rkj.z - L.z * floor(rkj.z / L.z + 0.5);

    //     // TODO: Provide the angle force calculation and assign forces to particles i, j, and k

    //     f[i].x += fi.x;
    //     f[i].y += fi.y;
    //     f[i].z += fi.z;
    //     f[j].x -= (fi.x + fk.x);
    //     f[j].y -= (fi.y + fk.y);
    //     f[j].z -= (fi.z + fk.z);
    //     f[k].x += fk.x;
    //     f[k].y += fk.y;
    //     f[k].z += fk.z;
    // }

    return Epot;  // Return the potential energy due to angle-bend interactions
}

// This function calculates dihedral-torsion forces based on the positions of four connected particles.
// It uses the minimum image convention and computes the forces resulting from the dihedral-torsion interaction.
double calculate_forces_dihedral(struct Parameters *p_parameters, struct Vectors *p_vectors){
    double Epot = 0.0;
    struct Dihedral *dihedrals = p_vectors->dihedrals;
    size_t num_dihedrals = p_vectors->num_dihedrals;
    struct Vec3D *f = p_vectors->f;
    struct Vec3D *r = p_vectors->r;
    struct Vec3D L = p_parameters->L;
    struct Vec3D rij, rkj, rkl;
    struct Vec3D fi = {0.0, 0.0, 0.0}, fk = {0.0, 0.0, 0.0}, fl = {0.0, 0.0, 0.0};

    // // Loop through each dihedral and calculate the forces
    // for (size_t q = 0; q < num_dihedrals; ++q)
    // {
    //     size_t i = dihedrals[q].i;
    //     size_t j = dihedrals[q].j;
    //     size_t k = dihedrals[q].k;
    //     size_t l = dihedrals[q].l;

    //     // Apply the minimum image convention for calculating distances
    //     rij.x = r[i].x - r[j].x;
    //     rij.x = rij.x - L.x * floor(rij.x / L.x + 0.5);
    //     rij.y = r[i].y - r[j].y;
    //     rij.y = rij.y - L.y * floor(rij.y / L.y + 0.5);
    //     rij.z = r[i].z - r[j].z;
    //     rij.z = rij.z - L.z * floor(rij.z / L.z + 0.5);

    //     rkj.x = r[k].x - r[j].x;
    //     rkj.x = rkj.x - L.x * floor(rkj.x / L.x + 0.5);
    //     rkj.y = r[k].y - r[j].y;
    //     rkj.y = rkj.y - L.y * floor(rkj.y / L.y + 0.5);
    //     rkj.z = r[k].z - r[j].z;
    //     rkj.z = rkj.z - L.z * floor(rkj.z / L.z + 0.5);

    //     rkl.x = r[l].x - r[k].x;
    //     rkl.x = rkl.x - L.x * floor(rkl.x / L.x + 0.5);
    //     rkl.y = r[l].y - r[k].y;
    //     rkl.y = rkl.y - L.y * floor(rkl.y / L.y + 0.5);
    //     rkl.z = r[l].z - r[k].z;
    //     rkl.z = rkl.z - L.z * floor(rkl.z / L.z + 0.5);

    //     // TODO: Provide the dihedral-torsion force calculation and assign forces to particles i, j, k, and l

    // }

    return Epot;  // Return the potential energy due to dihedral-torsion interactions
}

// This function calculates non-bonded forces between particles using the neighbor list.
// The total potential energy and forces are calculated using the three DPD forces.
double calculate_forces_nb(struct Parameters *p_parameters, struct Nbrlist *p_nbrlist, struct Vectors *p_vectors){
    struct DeltaR rij;
    struct Pair *nbr = p_nbrlist->nbr;
    const size_t num_nbrs = p_nbrlist->num_nbrs;
    struct Vec3D *f = p_vectors->f;
    size_t num_part = p_parameters->num_part;
    struct Vec3D *v = p_vectors->v;

    double gamma = p_parameters->gamma;
    double sigma = p_parameters->sigma;
    double r_cut = p_parameters->r_cut;
    double dt = p_parameters->dt;

    double aij = 0.0;
    double Epot = 0.0;
    struct Vec3D df = {0.0, 0.0, 0.0};
    double fc, fd, fr, ftot, wD, wR, r, zeta;

    // Loop through the neighbor list and calculate the forces for each particle pair
    for (size_t k = 0; k < num_nbrs; k++) {
        rij = nbr[k].rij;
        size_t i = nbr[k].i;
        size_t j = nbr[k].j;
        r = sqrt(rij.sq);
        zeta = generateZeta();

        int type_i = p_vectors->type[i];
        int type_j = p_vectors->type[j];

        // Compute forces if the distance is smaller than the cutoff distance
        if (r < r_cut) {

            if (type_i == 0 && type_j == 0) { // A-A interaction
                aij = p_parameters->a_AA;
            } else if (type_i == 1 && type_j == 1) { // B-B interaction
                aij = p_parameters->a_BB;
            } else { // A-B interaction
                aij = p_parameters->a_AB;
            }

            // Weight functions
            wR = 1 - r;
            wD = pow(wR, 2);

            // // Compute conservative force FC
            fc = aij * (1 - r) / r;

            // Compute dissipative force FD
            fd = -gamma * wD * (rij.x * (v[i].x - v[j].x) + rij.y * (v[i].y - v[j].y) + rij.z * (v[i].z - v[j].z)) / rij.sq;

            // Compute random force FR
            fr = sigma * wR * zeta / (r * sqrt(dt));

            // Calculate total force
            ftot = (fc + fd + fr);
            df.x = ftot * rij.x;
            df.y = ftot * rij.y;
            df.z = ftot * rij.z;

            // Add to total force on particles i and j
            f[i].x += df.x;
            f[i].y += df.y;
            f[i].z += df.z;
            f[j].x -= df.x;
            f[j].y -= df.y;
            f[j].z -= df.z;

            // Calculate potential energy
            Epot += 0.5 * aij * pow(1 - r, 2);
        }
    }

    return Epot;  // Return the potential energy due to non-bonded interactions
}

// Generate zeta value for DPD forces between -sqrt(3) and sqrt(3)
double generateZeta() {
    double u = (double)rand() / RAND_MAX;
    double zeta = (2 * sqrt(3) * u) - sqrt(3);
    return zeta;
}

// Calculate the pressure of the system based on the virial theorem.
// double calculate_pressure(struct Parameters *p_parameters, struct Vectors *p_vectors) {
//     struct Vec3D *v = p_vectors->v;
//     struct Vec3D *r = p_vectors->r;
//     struct Vec3D *f = p_vectors->f;
//     size_t num_part = p_parameters->num_part;
//     double volume = p_parameters->L.x * p_parameters->L.y * p_parameters->L.z;

//     double virial_sum = 0.0;
//     double rho = (double)num_part / volume; // Number density

//     // Compute virial sum
//     for (size_t i = 0; i < num_part; i++) {
//         for (size_t j = i + 1; j < num_part; j++) {
//             // Distance vector between particles i and j
//             double rx = r[i].x - r[j].x;
//             double ry = r[i].y - r[j].y;
//             double rz = r[i].z - r[j].z;

//             // Force vector between particles i and j
//             double fx = f[i].x - f[j].x;
//             double fy = f[i].y - f[j].y;
//             double fz = f[i].z - f[j].z;

//             // Virial term contribution
//             virial_sum += (fx * rx + fy * ry + fz * rz);
//         }
//     }

//     // Average virial term
//     double num_pairs = (num_part * (num_part - 1)) / 2.0;
//     double virial_term = virial_sum / (3.0 * volume * num_pairs); // Average over pairs

//     // Total pressure
//     double P = rho + virial_term;
//     printf("Pressure: %f\n", P);
//     return P;
// }

void calculate_density_profile(FILE *density_file, struct Parameters *p_parameters, struct Vectors *p_vectors, size_t step) {
    struct Vec3D *r = p_vectors->r;
    size_t num_part = p_parameters->num_part;
    struct Vec3D L = p_parameters->L;
    double bin_width = 0.1;
    size_t num_bins = (size_t)(L.z / bin_width);

    // Allocate memory for density profiles
    double *density_profile_A = (double *)calloc(num_bins, sizeof(double));
    double *density_profile_B = (double *)calloc(num_bins, sizeof(double));
    double *density_profile_total = (double *)calloc(num_bins, sizeof(double));

    double bin_volume = L.x * L.y * bin_width;  // Volume of a slab along z-axis

    // Calculate density profiles
    for (size_t i = 0; i < num_part; i++) {
        size_t bin = (size_t)(r[i].z / bin_width);

        // Ensure particles stay within valid bin range
        if (bin >= num_bins) {
            continue; // Skip particles out of the bin range
        }

        // Count particles in each bin
        density_profile_total[bin] += 1.0;
        if (p_vectors->type[i] == 0) {
            density_profile_A[bin] += 1.0;
        } else if (p_vectors->type[i] == 1) {
            density_profile_B[bin] += 1.0;
        }
    }

    // Normalize counts by the bin volume to get densities
    for (size_t i = 0; i < num_bins; i++) {
        density_profile_A[i] /= bin_volume;
        density_profile_B[i] /= bin_volume;
        density_profile_total[i] /= bin_volume;
    }

    // Output total density profile
    
    for (size_t i = 0; i < num_bins; i++) {
        fprintf(density_file, "%lu\t\t%f\t\t%f\t\t%f\t\t%f\n", (long unsigned)step, (i + 0.5) * bin_width, density_profile_A[i], density_profile_B[i], density_profile_total[i]);
    }

    // Free allocated memory
    free(density_profile_A);
    free(density_profile_B);
    free(density_profile_total);
}