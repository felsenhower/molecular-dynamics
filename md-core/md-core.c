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

#include "md-core.h"

void seed() {
   struct timeval time;
   gettimeofday(&time, NULL);
   srand48(((unsigned long long)time.tv_sec * 1000000) + time.tv_usec);
}

/**
 * The Lennard Jones potential between two molecules
 * @param r the distance between the two molecules
 * @return the resulting potential
 */
static double lennjo_potential(double r) {
  return 4 * (1/pow(r,12) - 1/pow(r,6));
}

/**
 * The force that gets induced by the Lennard jones potential
 * @param r the distance between the two molecules
 * @return The resulting force
 */
static double lennjo_force(double r) {
  return 24 * (2/pow(r,13) - 1/pow(r,7));
}

/**
 * Returns a vector from the given coordinates (that also contains the length
 * of the vector).
 * @param x the x coordinate
 * @param y the y coordinate
 * @return the resulting vector
 */
static Vector build_vector(double x, double y) {
  Vector result;
  result.coordinates.x = x;
  result.coordinates.y = y;
  result.length = sqrt(x*x + y*y);
  return result;
}

/**
 * Determines the distance between two molecules.
 * @param first the first position
 * @param second the second position
 * @param L the length of the square distance
 * @param periodic_boundaries whether to use periodic boundaries
 * @return the distance vector
 */
static Vector distance(const Coordinates *first, const Coordinates *second,
                       SimParams *params) {
  double dx = second->x - first->x;
  double dy = second->y - first->y;
  if (params->periodic_boundaries) {
    double L = params->L;
    if (fabs(dx) > L / 2) {
      dx = dx - copysign(1.0, dx) * L;
    }
    if (fabs(dy) > L / 2) {
      dy = dy - copysign(1.0, dy) * L;
    }
  }
  return build_vector(dx, dy);
}

/**
 * Gathers the kinetic energy, potential energy, total energy and temperature
 * of the system and writes that data back into their corresponding members
 * of data.
 * @param data The simulation data
 * @param params The simulation parameters
 */
void get_statistics(SimData *data, SimParams *params) {
  double Epot = 0.0;
  double Ekin = 0.0;
  for (uint64_t i=0; i < params->N; i++) {
    for (uint64_t j=i+1; j < params->N; j++) {
      Vector dist = distance(&(data->molecules[i].curr_pos),
                             &(data->molecules[j].curr_pos), params);
      if(params->dist_threshold == 0.0 || dist.length <= params->dist_threshold){
        Epot += lennjo_potential(dist.length);
      }
    }
    Ekin += (data->molecules[i].vel.x * data->molecules[i].vel.x
             + data->molecules[i].vel.y * data->molecules[i].vel.y) / 2;
  }
  data->kinetic_energy = Ekin;
  data->potential_energy = Epot;
  data->energy = Ekin + Epot;
  data->temperature = Ekin / params->N;
}

/**
 * Initializes the molecules array
 * @param params the simulation parameters to use
 * @return the molecules
 */
static Molecule *init_molecules(SimParams *params) {
  Molecule *molecules = malloc(params->N * sizeof(Molecule));
  ASSERT(molecules, "Could not initialize molecules");
  uint64_t grid_size = ceil(sqrt(params->N));
  double grid_dist = params->L / ((double)grid_size);

  for (uint64_t i = 0; i < params->N; i++) {
    double x_curr = ((((double)(i % grid_size)) + 0.5)
                    + params->max_deviation * (drand48() - 0.5) * sqrt(2)) * grid_dist;
    double y_curr = ((((double)(i / grid_size)) + 0.5)
                    + params->max_deviation * (drand48() - 0.5) * sqrt(2)) * grid_dist;
    double v_x = params->max_v0 * (drand48() - 0.5) * sqrt(2);
    double v_y = params->max_v0 * (drand48() - 0.5) * sqrt(2);
    double x_prev = x_curr - v_x * params->timestep;
    double y_prev = y_curr - v_y * params->timestep;
    Molecule molecule;
    molecule.curr_pos.x = x_curr;
    molecule.curr_pos.y = y_curr;
    molecule.prev_pos.x = x_prev;
    molecule.prev_pos.y = y_prev;
    molecule.vel.x = v_x;
    molecule.vel.y = v_y;
    molecules[i] = molecule;
  }
  return molecules;
}

void calculate_forces(SimData *data, SimParams *params){
  #pragma omp for simd
  for (uint64_t i = 0; i < params->N-1; i++) {
    for (uint64_t j = i + 1; j < params->N; j++) {
      Vector dist = distance(&(data->molecules[i].curr_pos),
                             &(data->molecules[j].curr_pos),
                             params);
  
      Coordinates *force_ij = &(data->forces_buffer[i][j]);
      Coordinates *force_ji = &(data->forces_buffer[j][i]);
  
      if(dist.length > 0.0 && (params->dist_threshold == 0.0
         || dist.length <= params->dist_threshold)){
        double force = lennjo_force(dist.length);
        double force_x = force * dist.coordinates.x / dist.length;
        double force_y = force * dist.coordinates.y / dist.length;
        force_ij->x = force_x;
        force_ij->y = force_y;
        force_ji->x = -force_x;
        force_ji->y = -force_y;
      } else {
        force_ij->x = 0.0;
        force_ij->y = 0.0;
        force_ji->x = 0.0;
        force_ji->y = 0.0;
      }
    }
  }
}

/**
 * Calculates the total force F_i that affects a particle i. This is just the sum
 * of all forces F_j,i that act on a particle.
 * @param data the simulation data
 * @param params the simulation params
 * @param i The index of the particle to get the force for.
 * @return the resulting force.
 */
static Coordinates get_total_force(SimData *data, SimParams *params, uint64_t i) {
  Coordinates result = {.x = 0.0, .y = 0.0};
  for (uint64_t j = 0; j < params->N; j++) {
    result.x += data->forces_buffer[j][i].x;
    result.y += data->forces_buffer[j][i].y;
  }
  return result;
}

/**
 * Corrects a position so that 0 <= position <= L.
 * @param curr the current position that is going to be positioned the same way
 *        relatively to next that it was before the correction was applies. This
 *        might include putting curr on the other side of next with the same
 *        distance if the movement of direction is being flipped during
 *        elastic reflection (if periodic_boundaries == false).
 * @param next the next position (the one to correct). If periodic_boundaries,
 *        next will be teleported to the other side of the lattice. Otherwise,
 *        it will be reflected elastically.
 * @param params the simulation params, including L and periodic_boundaries.
 */
static void fix_boundaries(double *curr, double *next, SimParams *params) {
  double L = params->L;
  if (*next > L || *next < 0) {
    // We'll need this to correct the current position later
    double curr_to_next = *next - *curr;
    if (params->periodic_boundaries) {
      if (*next > L) {
        // fmod should return a value between 0.0 and L here.
        *next = fmod(*next, L);
      } else {
        // fmod should return a value between -L and 0.0 here, so position should
        // mathematically be between 0.0 and and L, but numerical errors do happen.
        *next = L + fmod(*next, L);
        // The >= is there to catch numerical errors. The distance to L should be
        // at the order of magnitude of machine accuracy. Position can not be greater
        // than L mathematically at this point.
        if (*next >= L) {
          *next = 0.0;
        }
      }
    } else { /* Reflective walls */
      double dist_to_border = *next > L ? L - *curr : *curr;
      // The distance that remains after the first collision
      double rem_initial = fabs(curr_to_next) - dist_to_border;
      // The distance that remains after the last collision, if there is going
      // to be more than one, according to the velocity.
      double rem_final = fmod(rem_initial, L);
      // We'll need to turn the direction of movement around, if the particle
      // reflects an odd amount of times with any walls.
      bool keep_direction = (bool)((int)((rem_initial - rem_final) / L) % 2);
      // Determine at which border the particle is going to be placed
      if (keep_direction == (*next > L)) {
        *next = rem_final;
      } else {
        *next = L - rem_final;
      }
      // Determine on which side curr is going to be placed
      if (!keep_direction) {
        curr_to_next *= -1;
      }
    }
    // Correct curr.
    *curr = *next - curr_to_next;
  }
}


void apply_forces(SimData *data, SimParams *params){
  double squared_timestep = params->timestep * params->timestep;
  #pragma omp for simd
  for(uint64_t i = 0; i < params->N; i++) {

    Coordinates force = get_total_force(data, params, i);

    Coordinates *prev_pos = &(data->molecules[i].prev_pos);
    Coordinates *curr_pos = &(data->molecules[i].curr_pos);

    Coordinates next_pos = {.x = 2 * curr_pos->x - prev_pos->x
                            + force.x * squared_timestep,
                            .y = 2 * curr_pos->y - prev_pos->y
                            + force.y * squared_timestep};

    data->molecules[i].vel.x = (next_pos.x - prev_pos->x)
                                     / (2 * params->timestep);
    data->molecules[i].vel.y = (next_pos.y - prev_pos->y)
                                     / (2 * params->timestep);

    fix_boundaries(&(curr_pos->x), &(next_pos.x), params);
    fix_boundaries(&(curr_pos->y), &(next_pos.y), params);

    *prev_pos = *curr_pos;
    *curr_pos = next_pos;
  }
}

void enable_heat_bath(SimData *data, double target_temperature, uint64_t steps) {
  ASSERT(steps, "Steps must be greater than 0");

  double alpha_tilde = 3.0;
  HeatBath heat_bath;
  heat_bath.alpha_tilde = alpha_tilde;
  heat_bath.steps_to_the_alpha_tilde = pow(steps, alpha_tilde);
  heat_bath.steps = steps;
  heat_bath.current_step = 1;
  heat_bath.target_temperature = target_temperature;

  data->heat_bath = heat_bath;
  data->heat_bath_enabled = true;
}

void multiply_velocity(double *curr_pos, double *prev_pos, double *vel, double factor) {
  *prev_pos = *curr_pos - (*curr_pos - *prev_pos) * factor;
  *vel *= factor;
}

void apply_heat_bath(SimData *data, SimParams *params) {
  HeatBath *heat_bath = &(data->heat_bath);
  
  double foo1 = pow(heat_bath->current_step, heat_bath->alpha_tilde);
  double foo2 = foo1 / heat_bath->steps_to_the_alpha_tilde;
  double foo3 = heat_bath->target_temperature / data->temperature;
  double foo4 = pow(foo3, foo2);
  
  for (uint64_t i = 0; i < params->N; i++) {
    Coordinates *prev_pos = &(data->molecules[i].prev_pos);
    Coordinates *curr_pos = &(data->molecules[i].curr_pos);
    Coordinates *vel = &(data->molecules[i].vel);
    multiply_velocity(&(curr_pos->x), &(prev_pos->x), &(vel->x), foo4);
    multiply_velocity(&(curr_pos->y), &(prev_pos->y), &(vel->y), foo4);
  }  
  if (++(heat_bath->current_step) > heat_bath->steps) {
    data->heat_bath_enabled = false;
  }
}

void step(SimData *data, SimParams *params) {
  #pragma omp parallel
  calculate_forces(data, params);
  apply_forces(data, params);
  get_statistics(data, params);

  if (data->heat_bath_enabled) {
    apply_heat_bath(data, params);
  }
}

SimData init_data(SimParams *params) {
  SimData result;
  result.molecules = init_molecules(params);
  result.forces_buffer = malloc(sizeof(Coordinates*) * params->N);
  ASSERT(result.forces_buffer, "Could not allocate memory for forces buffer");

  for (uint64_t i = 0; i < params->N; i++) {
    result.forces_buffer[i] = malloc(sizeof(Coordinates) * params->N);
    ASSERT(result.forces_buffer[i], "Could not allocate memory for forces buffer");
    result.forces_buffer[i][i].x = 0.0;
    result.forces_buffer[i][i].y = 0.0;
  }
  
  result.heat_bath_enabled = false;
  
  return result;
}

SimParams init_params(uint64_t number_of_threads, uint64_t N,
                      double L, double timestep, double max_deviation,
                      double max_v0, bool periodic_boundaries, double dist_threshold,
                      uint64_t iterations) {
  SimParams result;
  result.number_of_threads = number_of_threads;
  result.N = N;
  result.L = L;
  result.timestep = timestep;
  result.max_deviation = max_deviation;
  result.max_v0 = max_v0;
  result.periodic_boundaries = periodic_boundaries;
  result.dist_threshold = dist_threshold;
  result.iterations = iterations;
  omp_set_num_threads(number_of_threads);
  return result;
}

void finalize(SimData *data, SimParams *params) {
  free(data->molecules);
  for (uint64_t i = 0; i < params->N; i++) {
    free(data->forces_buffer[i]);
  }
  free(data->forces_buffer);
}


