#ifndef NBRLIST_H_
#define NBRLIST_H_

#include <stddef.h>

/**
 * @brief Allocate arrays needed to store the cell-linked-list data
 *
 * @param p_parameters
 * @param p_cellist
 */
void alloc_celllist(struct Parameters *p_parameters, struct Celllist *p_celllist);

/** 
 * @brief Free arrays used for the cell-linked-list
 *
 * @param p_cellist
 */
void free_celllist(struct Celllist *p_cellist);

/**
 * @brief Build the cell-linked-list
 * 
 * @param p_parameters 
 * @param p_vectors 
 * @param p_celllist 
 */
void build_celllist(struct Parameters *p_parameters, struct Vectors * p_vectors, struct Celllist *p_celllist);

/**
 * @brief Allocate arrays needed to store the neighbor list 
 * 
 * @param p_parameters 
 * @param p_nbrlist 
 */
void alloc_nbrlist(struct Parameters *p_parameters, struct Nbrlist *p_nbrlist);

/**
 * @brief Free arrays use to store the neighbor list
 * 
 * @param p_nbrlist 
 */
void free_nbrlist(struct Nbrlist *p_nbrlist);


/**
 * @brief is_connected_12 returns 1 if i and j are 1-2 connected
 * 
 * @param i, particle index
 * @param j, particle index
 * @param p_nbrlist pointer to neighbor list, members used head12, pairs12
 * @return int Returns 1 if particles are connected via a 1-2 (bonded) interaction, excluding them from non-bonded force calculations when p_parameters->exclude_12_nb == 1
 */
int is_connected_12(size_t i, size_t j, struct Nbrlist *p_nbrlist);

/**
 * @brief is_connected_13 returns 1 if i and j are 1-3 connected
 * 
 * @param i, particle index
 * @param j, particle index
 * @param p_nbrlist pointer to neighbor list, members used head13, pairs13
 * @return int Returns 1 when particle i and j are 1-3 connected, , excluding them from non-bonded force calculations when p_parameters->exclude_13_nb == 1
 */
int is_connected_13(size_t i, size_t j, struct Nbrlist *p_nbrlist);

/**
 * @brief is_connected_14 returns 1 if i and j are 1-4 connected
 * 
 * @param i, particle index
 * @param j, particle index
 * @param p_nbrlist pointer to neighbor list, members used head13, pairs13
 * @return int Returns 1 when particle i and j are 1-4 connected, 0 otherwise
 */
int is_connected_14(size_t i, size_t j, struct Nbrlist *p_nbrlist);

/**
 * @brief Build the neighbor list
 * 
 * @param p_parameters used members: rcut, rshell
 * @param p_vectors used members: r
 * @param p_nbrlist pointer to neighbor list
 */
void build_nbrlist(struct Parameters *p_parameters, struct Vectors *p_vectors, struct Nbrlist *p_nbrlist);

/**
 * @brief Update the neighbor list
 * Checks if the neigbor lists needs to be rebuild be compairing the maximum displacement with rcut+rshell. 
 * If so build_nbrlist is called. If not, it updates positions and squared-distances of all pairs
 * 
 * @param p_parameters 
 * @param p_vectors 
 * @param p_nbrlist 
 * @return int Returns 1 if nbrlist is rebuild and 0 if it is only updated.
 */
int update_nbrlist(struct Parameters *p_parameters, struct Vectors *p_vectors, struct Nbrlist *p_nbrlist);

#endif /* NBRLIST */