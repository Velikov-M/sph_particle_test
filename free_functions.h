#ifndef FREE_FUNCTIONS_H

#define FREE_FUNCTIONS_H

#include "Particle.h"
#include <vector>

double weight_function(double r, double b);

double ddr_weight_fun(double r, double b);

double calculate_distance(Particle& p1, Particle& p2);

double calculate_density(std::vector<Particle>& particles, int id_part, int num);

std::vector<std::vector<double>> calc_eulier_pres_tensor(std::vector<double>& densities, int id_part, int num);

std::vector <double> calc_grad_weight_fun(std::vector<Particle>& particles, int id_part, int j_part);

std::vector<std::vector<double>> calc_deform_tensor(std::vector<Particle>& particles, std::vector<double>& densities, int id_part, int num);

std::vector<std::vector<double>> calc_newton_pres_tensor(std::vector<Particle>& particles, std::vector<double>& densities, int id_part, int num);

std::vector<double> calculate_acceleration(std::vector<Particle>& particles, std::vector<std::vector<std::vector<double>>>& pres_tensors, std::vector<double>& densities, int id_part, int num, bool gravity_flag);

void apply_friction(std::vector<Particle>& particles, int id_part);

double calc_kinetic_nrg(std::vector<Particle>& particles, int start, int end);

double calc_momentum(std::vector<Particle>& particles, int start, int end);

double calc_norm_of_impulse(std::vector<Particle>& particles, int start, int end);

#endif