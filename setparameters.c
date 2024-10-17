#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "constants.h"
#include "structs.h"

// Set the parameters of this simulation
void set_parameters(struct Parameters *p_parameters)
{
  // Units of the system
  p_parameters->kT = 1.0;                                   // Unit of energy (kT)
  p_parameters->mass = 1.0;                                 // Unit of mass is particle mass (m)
  p_parameters->r_cut = 1.0;                                // Unit of length is cutoff distance (r_cut)

  // Simulation parameters from paper
  p_parameters->num_part = 3840;                            // Number of particles
  p_parameters->num_dt_steps = 100000;                      // Number of time steps
  p_parameters->L = (struct Vec3D){8, 8, 20};               // Box sizes in 3 direction
  p_parameters->dt = 0.04;                                  // Integration time step
  p_parameters->gamma = 4.5;                                // DPD gamma
  p_parameters->sigma = 3.0;                                // DPD sigma

  // Binary Mixture
  p_parameters->a_AA = 25.0;                                // DPD maximum repulsion AA
  p_parameters->a_BB = 25.0;                                // DPD maximum repulsion BB
  p_parameters->a_AB = 37.0;                                // DPD maximum repulsion AB
  p_parameters->ratio_AB = 0.5;                             // Ratio of A and B particles

  // Number of monomers per chain
  p_parameters->chain_length = 4;                           // Number of monomers per chain

  p_parameters->exclude_12_nb = 0;                          // 1-2 connected atoms exluded from non-bonded interactions 
  p_parameters->exclude_13_nb = 0;                          // 1-3 connected atoms exluded from non-bonded interactions    
  p_parameters->r_shell = 0.2;                              // Shell thickness for neighbor list
  p_parameters->num_dt_pdb = 100;                           // Number of time steps in between pdb outputs
  strcpy(p_parameters->filename_pdb, "trajectories");       // Filename (without extension) for pdb file
  p_parameters->rescale_output = 1;                         // Factor used to rescale output lengthscale (Most visualisation programs identify bonds based on distances of order 1)
  p_parameters->load_restart = 0;                           // If equal 1 restart file is loaded
  strcpy(p_parameters->restart_in_filename, "restart.dat"); // Filename for loaded restart file
  p_parameters->num_dt_restart = 1000;                      // Number of time steps between saves
  strcpy(p_parameters->restart_out_filename, "restart.dat");// Filename for saved restart file


  if (p_parameters->r_cut > p_parameters->L.x / 2.0)
    fprintf(stderr, "Warning! r_cut > Lx/2");
  if (p_parameters->r_cut > p_parameters->L.y / 2.0)
    fprintf(stderr, "Warning! r_cut > Ly/2");
  if (p_parameters->r_cut > p_parameters->L.z / 2.0)
    fprintf(stderr, "Warning! r_cut > Lz/2");
}
