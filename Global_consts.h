#ifndef GLOBAL_CONSTS_H

#define GLOBAL_CONSTS_H

const double mass = 1; // mass of one particle, for now i am assuming be equal
const double B = 6428;  // next 3 global vars are consts in equation of state for liquid (Cole form, see Monaghan, SPH, 2005) 6428.6
const double gamma = 7.0; //mb should think more about B and rho_0;
const double rho_0 = 1.0;
//const double rho_0 = 0.0105;
const double mu = 50.0; // dynamic viscosity
const double b = 1.5; // domain of interest
const double g = 9.81 * 300; // gravity const, increased for tests
const double alpha = 0.02; //linear dissipation coefficient
const double Ene_bond = 0.01; // next 2 global vars ara parametrs for Leonard-Jones potential (interaction water-stone)
const double a_bond = 1.5;
const double mass_stone = 40.0; // kg
const double moment_inertia_stone = 25.0; // kg * m ^ 2
const double length_of_area = 60; // important for BC!

#endif

