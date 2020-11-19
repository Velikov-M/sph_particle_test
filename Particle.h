#ifndef PARTICLE_H
#define PARTICLE_H

#include <vector>

class Particle {
private:

	std::vector<double> position;
	std::vector<double> velosity;

public:
	Particle();
	void set_position(double x1, double x2);

	inline std::vector<double> get_position() { return position;};

	void print_position();

	void set_velosity(double v1, double v2);

	inline std::vector<double> get_velosity() { return velosity;};

	void print_velosity();

};

#endif
