#ifndef MEMORY_H_
#define MEMORY_H_

#include "structs.h"

/**
 * @brief Allocate the arrays in 'vectors' needed to store information of all particles
 * 
 * @param p_vectors struct with all vectors that need to be allocated
 * @param num_part number of particles
 */
void alloc_vectors(struct Vectors *p_vectors, size_t num_part);

/**
 * @brief Free all allocated vectors
 * 
 * @param p_vectors 
 */
void free_vectors(struct Vectors *p_vectors);

/**
 * @brief Allocate all variables needed in the MD simulation
 * 
 * @param p_parameters Contains information to be used to determine sizes
 * @param p_vectors Contains vectors to be allocated
 * @param p_nbrlist neighbor list with arrays that need initialization
 */
void alloc_memory(struct Parameters *p_parameters, struct Vectors *p_vectors, struct Nbrlist *p_nbrlist);

/**
 * @brief free all allocated memory
 * 
 * @param p_vectors all vectors used in the program
 * @param p_nbrlist struct that stores neighbor lists
 */
void free_memory(struct Vectors *p_vectors, struct Nbrlist *p_nbrlist);

#endif /* MEMORY_H_ */