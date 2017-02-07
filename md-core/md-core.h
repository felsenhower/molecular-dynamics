/******************************************************************************/
/*                                                                            */
/*                             Molecular dynamics                             */
/*                             ==================                             */
/*                                                                            */
/* Copyright (c) 2017 Ruben Felgenhauer, Leonhard Reichenbach                 */
/*                                                                            */
/* This file is released under the MIT license.                               */
/*                                                                            */
/******************************************************************************/

#ifndef MDCORE_H
#define MDCORE_H

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>
#include <stdbool.h>
#include <pthread.h>
#include <inttypes.h>

#define ASSERT(COND, MSG)\
  if(!(COND)){\
    fprintf(stderr,"Error: \"%s\" in %s, line:%d\n", (MSG),__FILE__, __LINE__);\
    exit(EXIT_FAILURE);\
  }

typedef struct {
  double x;
  double y;
} Coordinates;

/**
  Struct for a molecule, describing its position and velocity
*/
typedef struct {
  Coordinates curr_pos; /**< The current position */
  Coordinates prev_pos; /**< The previous position */
  Coordinates vel;      /**< The (current) velocity */
} Molecule;

typedef struct {
  double target_temperature;
  uint64_t steps;
  uint64_t current_step;
  double alpha_tilde;
  double steps_to_the_alpha_tilde;
} HeatBath;

typedef struct {
  Molecule *molecules;         /**< The array of molecules */
  uint64_t **index_mappings;   /**< The 2d-array of index mappings */
  Coordinates **forces_buffer; /**< The buffer to store the forces inside */
  double kinetic_energy;       /**< The total kinetic energy of the system */
  double potential_energy;     /**< The total potential energy of the system */
  double energy;               /**< The energy of the system = E_kin + E_pot */
  double temperature;          /**< The temperature of the system */
  HeatBath heat_bath;
  bool heat_bath_enabled;
} SimData;

typedef struct {
  uint64_t number_of_threads; /**< The number of threads to use */
  uint64_t N;                 /**< The number of molecules */
  double L;                   /**< The width of the square lattice */
  double timestep;            /**< The timestep */
  double max_deviation;       /**< The maximum deviation of the molecules' positions */
  double max_v0;              /**< The maximum initial velocity of the molecules */
  bool periodic_boundaries;   /**< Whether to use periodic boundaries */
  double dist_threshold;      /**< The max distance to still calculate forces */
  uint64_t iterations;        /**< unused */
} SimParams;

typedef struct {
  Coordinates coordinates;
  double length;
} Vector;



/**
 * Seeds the drand48, erand48, jrand48, lrand48, mrand48, nrand48 methods
 * for pseudo-random number generation with micro-time.
 * Only works on POSIX systems.
 */
void seed();

/**
 * Initializes all simulation data (molecules, forces buffer, and index mappings)
 * @param params the simulation parameters
 * @return the simulation data
 */
SimData init_data(SimParams *params);

/**
 * Initializes the simulation parameters with the given values.
 * @param number_of_threads the number of threads to use
 * @param N the number of molecules
 * @param L the width of the square lattice
 * @param timestep the timestep
 * @param max_deviation the maximum deviation of the molecules' positions
 * @param max_v0 the maximum initial velocity of the molecules
 * @param periodic_boundaries whether to use periodic boundaries
 * @param dist_threshold the maximum distance at which forces get calculated
 * @param iterations unused param
 * @return a structure which contains all above parameters
 */
SimParams init_params(uint64_t number_of_threads, uint64_t N,
                      double L, double timestep, double max_deviation,
                      double max_v0, bool periodic_boundaries, double dist_threshold,
                      uint64_t iterations);

/**
  Frees all the allocated resources from initalize.
*/

/**
 * Frees all the resources that have been allocated at init_data
 * @param data the simulation data to delete
 * @param params the simulation parameters
 */
void finalize(SimData *data, SimParams *params);

/**
 * Do one simulation step:
 * 1) Calculate the forces and store them in the forces buffer
 * 2) Load the forces and apply them to the molecules
 * @param data the simulation data with molecules ans forces buffer
 * @param params the simulatiom params to use
 */
void step(SimData *data, SimParams *params);

/**
 * Enables a heat bath with the given arguments for the simulation
 * @param data the simulation data
 * @param target_temperature the temperature to adjust the system to
 * @param steps the number of steps to use
 * @param attack_coefficient The attack coefficient -1 < alpha < 1
 */
void enable_heat_bath(SimData *data, double target_temperature, uint64_t steps);


#endif // MDCORE_H
