#ifndef RANDOM_H_
#define RANDOM_H_

/**
 * @brief Draw a uniformly (pseudo) random number that is uniformly distributed between 0 and 1.
 * 
 * @return double random number in (0,1)
 */
double generate_uniform_random(void);

/**
 * @brief Generate Gaussian distributed random number with mean 0 and standard deviations 1
 * 
 * @return double Normal distributed random number.
 */
double gauss(void);

#endif /* RANDOM_H_ */
