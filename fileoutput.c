#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "constants.h"
#include "memory.h"
#include "structs.h"

// Write the particle positions to a pdf file
// The filename (without extension) is given by p_parameters->filename_pdb.
// If reset = 1 the data is written to the file deleting data it possibly contained.
// If reset = 0 the data is appended.
void record_trajectories_pdb(int reset, struct Parameters *p_parameters, struct Vectors *p_vectors, double time)
{
  FILE *fp_traj;
  char filename[1024];
  double rs = p_parameters->rescale_output;

  snprintf(filename, 1024, "%s%s", p_parameters->filename_pdb, ".pdb");
  if (reset == 1)
  {
    fp_traj = fopen(filename, "w");
  }
  else
  {
    fp_traj = fopen(filename, "a");
  }

  fprintf(fp_traj, "MODEL\n");
  fprintf(fp_traj, "REMARK TIME = %f\n", time);
  fprintf(fp_traj, "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %-10s%-3s\n", rs*p_parameters->L.x, rs*p_parameters->L.y, rs*p_parameters->L.z, 90.0, 90.0, 90.0, "P 1", "1");
      for (size_t i = 0; i < p_parameters->num_part; i++)
    {
        // Adjust the atom type and color based on the particle type
        char atom_type[3];
        char color[8];
        
        if (p_vectors->type[i] == 0)
        {
            snprintf(atom_type, sizeof(atom_type), "A"); // Use "A" for type A
        }
        else if (p_vectors->type[i] == 1)
        {
            snprintf(atom_type, sizeof(atom_type), "B"); // Use "B" for type B
        }

        // Include color as a comment or a separate field if the format supports it
        fprintf(fp_traj, "HETATM%5u  %s  UNK A   1    %8.3f%8.3f%8.3f  1.00  0.00           C ; %s\n", 
                (unsigned int)i % 100000, atom_type, 
                rs*p_vectors->r[i].x, rs*p_vectors->r[i].y, rs*p_vectors->r[i].z, color);
    }

  fprintf(fp_traj, "ENDMDL\n");

  fclose(fp_traj);
}

// Write the particle positions to a xyz file
// The filename (without extension) is given by p_parameters->filename_xyz.
// If reset = 1 the data is written to the file deleting data it possibly contained.
// If reset = 0 the data is appended.
void record_trajectories_xyz(int reset, struct Parameters *p_parameters, struct Vectors *p_vectors, double time)
{
  FILE *fp_traj;
  char filename[1024];
  double rs = p_parameters->rescale_output;

  snprintf(filename, 1024, "%s%s", p_parameters->filename_xyz, ".xyz");
  if (reset == 1)
  {
    fp_traj = fopen(filename, "w");
  }
  else
  {
    fp_traj = fopen(filename, "a");
  }

  fprintf(fp_traj, "%lu\n", p_parameters->num_part);
  fprintf(fp_traj, "time = %f\n", time);
  struct Vec3D *r = p_vectors->r;
  for (size_t i = 0; i < p_parameters->num_part; i++)
  {
    fprintf(fp_traj, "  C        %10.5f %10.5f %10.5f\n", rs*r[i].x, rs*r[i].y, rs*r[i].z);
  }

  fclose(fp_traj);
}

// save arrays in vectors to binary file
void save_restart(struct Parameters *p_parameters, struct Vectors *p_vectors)
{
  FILE* p_file = fopen( p_parameters->restart_out_filename, "wb");
  size_t num_part = p_parameters->num_part;
  size_t sz = num_part*sizeof(struct Vec3D);

  fwrite(&num_part, sizeof(size_t), 1, p_file);
  fwrite(p_vectors->r, sz, 1, p_file);
  fwrite(p_vectors->v ,sz, 1, p_file);
  fwrite(p_vectors->f, sz, 1, p_file);
  fclose(p_file);
}

// load arrays in vectors from binary file
void load_restart(struct Parameters *p_parameters, struct Vectors *p_vectors)
{
  FILE* p_file = fopen( p_parameters->restart_in_filename, "rb" );
  size_t num_part;
  fread(&num_part, sizeof(size_t), 1, p_file);
  size_t sz = num_part*sizeof(struct Vec3D);
  alloc_vectors(p_vectors,num_part);
  p_parameters->num_part = num_part;
  fread(p_vectors->r, sz, 1, p_file);
  fread(p_vectors->v, sz, 1, p_file);
  fread(p_vectors->f, sz, 1, p_file);
  fclose(p_file);
}

void collect_velocity_data(struct Parameters *p_parameters, struct Vectors *p_vectors, size_t step, int is_first_call) {
  size_t num_part = p_parameters->num_part;

  FILE *velocity_file;
  if (is_first_call) {
      velocity_file = fopen("velocity_data.txt", "w"); // Overwrite file on first call
      fprintf(velocity_file, "Step\tVelocity\n");
  } else {
      velocity_file = fopen("velocity_data.txt", "a"); // Append to file on subsequent calls
  }

  for (size_t i = 0; i < num_part; ++i) {
    double vx = p_vectors->v[i].x;
    double vy = p_vectors->v[i].y;
    double vz = p_vectors->v[i].z;
    double velocity = sqrt(vx * vx + vy * vy + vz * vz);
    fprintf(velocity_file, "%d\t%f\n", (long unsigned)step, velocity);
  }

  fclose(velocity_file);
}

