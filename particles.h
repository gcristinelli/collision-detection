#ifndef PARTICLES_H
#define PARTICLES_H

#include "utils.h"

struct Material {
    double density;
    double normal_stiffness;
    double normal_damping;
    double tangential_stiffness;
    double tangential_damping;
    double friction_coefficient;
};

struct Particle {
    Vec3 pos, vel, acc;
    double radius, mass;
};

#endif
