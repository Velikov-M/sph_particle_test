#ifndef CELL_FREE_FUNCTIONS_H

#define CELL_FREE_FUNCTIONS_H

// most of calculation functions must be specifically rewritten to apply cell calculations, as well as to apply periodic bcs.
//Though, It's a lot of copypaste in this function, but we still need separate functions

#include "Particle.h"
#include <vector>

double cell_calculate_density(std::vector<Particle>& particles, std::vector<std::vector<std::vector<int>>>& cells, int id_part, int i_cell, int j_cell, int i_end);

std::vector<std::vector<double>> cell_calc_newton_pt(std::vector<Particle>& particles, std::vector<double>& densities, std::vector<std::vector<std::vector<int>>>& cells, int id_part, int i_cell, int j_cell, int i_end);

std::vector<std::vector<double>> cell_calc_deform_tensor(std::vector<Particle>& particles, std::vector<double>& densities, std::vector<std::vector<std::vector<int>>>& cells, int id_part, int i_cell, int j_cell, int i_end);

std::vector<double> cell_calc_acc(std::vector<Particle>& particles, std::vector<std::vector<std::vector<double>>>& pres_tensors, std::vector<double>& densities, std::vector<std::vector<std::vector<int>>>& cells, int id_part, int i_cell, int j_cell, int i_end, bool gravity_flag);

std::vector<double> cell_calc_stone_accs(std::vector<Particle>& particles, std::vector<std::vector<std::vector<int>>>& cells, std::vector<double> stone_pos, int n_exterior, int n_stone, int i_end, bool gravity_flag);

#endif
