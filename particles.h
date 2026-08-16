#ifndef PARTICLES_H
#define PARTICLES_H

#include <string>
#include <vector>

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

std::vector<Particle> loadParticlesCsv(const std::string& filename,
                                       const Vec3& initial_acceleration);

#endif
