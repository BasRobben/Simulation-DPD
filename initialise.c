#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "constants.h"
#include "structs.h"
#include "random.h"
#include "initialise.h"

// This function initializes the particle types.
// Particle type 0 is assigned to all particles for now, but this could be modified
// later to initialize different types (e.g., methane and ethane in a binary mixture).
void initialise_types(struct Parameters *p_parameters, struct Vectors *p_vectors){
    size_t num_part = p_parameters->num_part;               // Number of particles
    double ratio_AB = p_parameters->ratio_AB;               // Ratio A/B
    double num_A = num_part * ratio_AB;                     // Number of particles of type A
    for (size_t i = 0; i < num_part; i++) {
        if (i < num_A) {
            p_vectors->type[i] = 0;       // Particle type 0 corresponding to A
        }
        else {
            p_vectors->type[i] = 1;       // Particle type 1 corresponding to B
        }
    }
}

// This function initializes the bond connectivity between particles.
// This will be important for handling bonded interactions in the simulation.
void initialise_bond_connectivity(struct Parameters *p_parameters, struct Vectors *p_vectors){

    size_t chain_length = p_parameters->chain_length;       // Number of monomers per chain
    size_t num_part = p_parameters->num_part;               // Number of particles

    // Check if the number of particles is a multiple of chain length
    if (num_part % chain_length != 0){
        fprintf(stderr, "Error: Number of particles is not a multiple of chain length\n");
        exit(EXIT_FAILURE);
    }

    size_t num_chains = num_part / chain_length;            // Number of chains
    size_t num_bonds = num_chains * (chain_length - 1);     // Total number of bonds
    
    // Allocate memory for the bonds
    struct Bond *bonds = (struct Bond *)malloc(num_bonds * sizeof(struct Bond));
    if (bonds == NULL){
        fprintf(stderr, "Error: Memory allocation for bonds failed\n");
        exit(EXIT_FAILURE);
    }

    // Initialize bonds for each chain
    size_t bond_index = 0;  // Index to keep track of bonds in the bonds array
    for (size_t chain = 0; chain < num_chains; chain++){
        size_t start_particle = chain * chain_length;  // The starting particle of the current chain
        for (size_t i = 0; i < chain_length - 1; i++){  // Each chain has (chain_length - 1) bonds
            bonds[bond_index].i = start_particle + i;       // Set the i-th particle in the bond
            bonds[bond_index].j = start_particle + i + 1;   // Set the (i+1)-th particle in the bond
            bond_index++;  // Move to the next bond
        }
    }

    // Update the vectors structure
    p_vectors->num_bonds = num_bonds;
    p_vectors->bonds = bonds;
}


// This function initializes the molecular structure, including bonds, angles, and dihedrals,
// and computes 1-2, 1-3, and 1-4 bonded interactions for the neighbor list used in force calculations.
void initialise_structure(struct Parameters *p_parameters, struct Vectors *p_vectors, struct Nbrlist *p_nbrlist)
{
    initialise_bond_connectivity(p_parameters, p_vectors); // Initialize bonds
    
    struct Bond *bonds = p_vectors->bonds;
    size_t num_bonds = p_vectors->num_bonds;
    size_t num_part = p_parameters->num_part;
    size_t *cnt = (size_t *)calloc(num_part + 1, sizeof(size_t));
    size_t cnt_tot = 0;
    for (size_t i = 0; i < num_bonds; ++i)
    {
        ++cnt[bonds[i].i + 1];
        ++cnt[bonds[i].j + 1];
    }
    size_t *head12 = (size_t *)malloc((num_part + 1) * sizeof(size_t));
    head12[0] = 0;
    for (size_t i = 1; i <= num_part; ++i)
    {
        head12[i] = cnt[i] + head12[i - 1];
        cnt[i] = head12[i];
    }
    size_t *pairs12 = (size_t *)malloc(cnt[num_part] * sizeof(size_t));
    for (size_t i = 0; i < num_bonds; ++i)
    {
        pairs12[cnt[bonds[i].i]++] = bonds[i].j;
        pairs12[cnt[bonds[i].j]++] = bonds[i].i;
    }
    size_t num_angles = 0;
    for (size_t i = 1; i <= num_part; ++i)
    {
        size_t num_pairs = head12[i] - head12[i - 1];
        num_angles += (num_pairs * (num_pairs - 1)) / 2;
    }
    struct Angle *angles = (struct Angle *)malloc(num_angles * sizeof(struct Angle));
    size_t m = 0;
    for (size_t i = 0; i < num_part; ++i)
        for (size_t j = head12[i]; j <= head12[i + 1]; ++j)
            for (size_t k = j + 1; k < head12[i + 1]; ++k)
            {
                angles[m].i = pairs12[j];
                angles[m].j = i;
                angles[m].k = pairs12[k];
                ++m;
            }
    size_t num_dihdr = 0;
    for (size_t i = 0; i < num_bonds; ++i)
    {
        num_dihdr += (head12[bonds[i].i + 1] - head12[bonds[i].i] - 1) * (head12[bonds[i].j + 1] - head12[bonds[i].j] - 1);
    }
    struct Dihedral *dihedrals = (struct Dihedral *) malloc(num_dihdr * sizeof(struct Dihedral));
    size_t k = 0;
    for (size_t i = 0; i < num_bonds; ++i)
    {
        struct Dihedral dihdr;
        dihdr.j = bonds[i].i;
        dihdr.k = bonds[i].j;
        for (size_t j = head12[dihdr.j]; j < head12[dihdr.j + 1]; ++j)
        {
            if (pairs12[j] == dihdr.k)
                continue;
            dihdr.i = pairs12[j];
            for (size_t j = head12[dihdr.k]; j < head12[dihdr.k + 1]; ++j)
            {
                if (pairs12[j] != dihdr.j)
                {
                    dihdr.l = pairs12[j];
                    dihedrals[k++] = dihdr;
                }
            }
        }
    }
    if (p_parameters->exclude_13_nb)
    {
        for (size_t i = 0; i <= num_part; ++i)
            cnt[i] = 0;
        for (size_t i = 0; i < num_angles; ++i)
        {
            ++cnt[angles[i].i + 1];
            ++cnt[angles[i].k + 1];
        }
        size_t *head13 = (size_t *)malloc((num_part + 1) * sizeof(size_t));
        head13[0] = 0;
        for (size_t i = 1; i <= num_part; ++i)
        {
            head13[i] = cnt[i] + head13[i - 1];
            cnt[i] = head13[i];
        }
        size_t *pairs13 = (size_t *)malloc(cnt[num_part] * sizeof(size_t));
        for (size_t i = 0; i < num_angles; ++i)
        {
            pairs13[cnt[angles[i].i]++] = angles[i].k;
            pairs13[cnt[angles[i].k]++] = angles[i].i;
        }
        p_nbrlist->head13 = head13;
        p_nbrlist->pairs13 = pairs13;
    }

    for (size_t i = 0; i <= num_part; ++i)
        cnt[i] = 0;
    for (size_t i = 0; i < num_dihdr; ++i)
    {
        ++cnt[dihedrals[i].i + 1];
        ++cnt[dihedrals[i].l + 1];
    }
    size_t *head14 = (size_t *)malloc((num_part + 1) * sizeof(size_t));
    head14[0] = 0;
    for (size_t i = 1; i <= num_part; ++i)
    {
        head14[i] = cnt[i] + head14[i - 1];
        cnt[i] = head14[i];
    }
    size_t *pairs14 = (size_t *)malloc(cnt[num_part] * sizeof(size_t));
    for (size_t i = 0; i < num_dihdr; ++i)
    {
        pairs14[cnt[dihedrals[i].i]++] = dihedrals[i].l;
        pairs14[cnt[dihedrals[i].l]++] = dihedrals[i].i;
    }
    p_nbrlist->head14 = head14;
    p_nbrlist->pairs14 = pairs14;
    p_vectors->num_angles = num_angles;
    p_vectors->angles = angles;
    p_vectors->num_dihedrals = num_dihdr;
    p_vectors->dihedrals = dihedrals;
    if (p_parameters->exclude_12_nb)
    {
        p_nbrlist->head12 = head12;
        p_nbrlist->pairs12 = pairs12;
    }
    else
    {
        free(head12);
        free(pairs12);
    }
    free(cnt);
}


// This function initializes the simulation by calling subroutines to initialize
// particle types, positions, velocities, and bond connectivity.
void initialise(struct Parameters *p_parameters, struct Vectors *p_vectors, struct Nbrlist *p_nbrlist, size_t *p_step, double *p_time)
{
    initialise_types(p_parameters, p_vectors);  // Initialize particle types
    initialise_structure(p_parameters, p_vectors, p_nbrlist);  // Initialize structure (bonds, angles, dihedrals)
    srand(13);  // Seed random number generator
    initialise_positions(p_parameters, p_vectors);  // Initialize particle positions
    initialise_velocities(p_parameters, p_vectors);  // Initialize particle velocities
    *p_step = 0;   // Initialize step to zero
    *p_time = 0.0; // Initialize time to zero
    return;
}


// This function initializes particle positions on a cubic lattice.
// Particles are placed in a grid with spacing based on the number of particles and the box dimensions.
void initialise_positions(struct Parameters *p_parameters, struct Vectors *p_vectors){
    struct Vec3D dr;  // Displacement vector for positioning particles
    struct Index3D n; // Number of grid cells along each axis
    double dl;        // Lattice spacing
    int ipart = 0;    // Particle index

    // Calculate lattice spacing based on particle number and box dimensions
    dl = pow(p_parameters->L.x * p_parameters->L.y * p_parameters->L.z / ((double)p_parameters->num_part), 1.0 / 3.0);
    n.i = (int)ceil(p_parameters->L.x / dl);
    n.j = (int)ceil(p_parameters->L.y / dl);
    n.k = (int)ceil(p_parameters->L.z / dl);
    dr.x = p_parameters->L.x / (double)n.i;
    dr.y = p_parameters->L.y / (double)n.j;
    dr.z = p_parameters->L.z / (double)n.k;

    // Initialize positions of particles
    for (size_t k = 0; k < n.k; ++k) {
        for (size_t j = 0; j < n.j; ++j) {
            for (size_t i = 0; i < n.i; ++i, ++ipart) {
                if (ipart >= p_parameters->num_part) {
                    break;
                }
                p_vectors->r[ipart].x = (i + 0.5) * dr.x; // x-position
                p_vectors->r[ipart].y = (j + 0.5) * dr.y; // y-position

                // Assign z-position based on particle type
                if (p_vectors->type[ipart] == 0) {
                    // Particles A on the left side of the box (z = 0 to L.z / 2)
                    p_vectors->r[ipart].z = (k + 0.5) * dr.z / 2; // z-position for type A in lower half
                } else if (p_vectors->type[ipart] == 1) {
                    // Particles B on the right side of the box (z = L.z / 2 to L.z)
                    p_vectors->r[ipart].z = p_parameters->L.z / 2 + (k + 0.5) * dr.z / 2; // z-position for type B in upper half
                }
            }
        }
    }
}


// This function initializes the velocities of particles based on the Maxwell-Boltzmann distribution.
// The total momentum is also removed to ensure zero total momentum (important for stability).
void initialise_velocities(struct Parameters *p_parameters, struct Vectors *p_vectors)
{
    double sqrtktm = sqrt(p_parameters->kT / p_parameters->mass);
    struct Vec3D sumv = {0.0, 0.0, 0.0};  // Total velocity (to remove later)

    // Assign random velocities to each particle
    for (size_t i = 0; i < p_parameters->num_part; i++)
    {
        p_vectors->v[i].x = sqrtktm * gauss();
        p_vectors->v[i].y = sqrtktm * gauss();
        p_vectors->v[i].z = sqrtktm * gauss();
        sumv.x += p_vectors->v[i].x;
        sumv.y += p_vectors->v[i].y;
        sumv.z += p_vectors->v[i].z;
    }

    // Remove the average velocity to ensure zero total momentum
    sumv.x /= ((double)(p_parameters->num_part));
    sumv.y /= ((double)(p_parameters->num_part));
    sumv.z /= ((double)(p_parameters->num_part));
    for (size_t i = 0; i < p_parameters->num_part; i++)
    {
        p_vectors->v[i].x -= sumv.x;
        p_vectors->v[i].y -= sumv.y;
        p_vectors->v[i].z -= sumv.z;
    }
}