#ifndef XY_WORK_H

#define XY_WORK_H

#include "Particle.h"

void create_XY(std::vector<Particle>& particles, int num);

void add_ts_XY(std::vector<Particle>& particles, int num);

int read_initial_position(std::vector<Particle>& particles, int n);

#endif